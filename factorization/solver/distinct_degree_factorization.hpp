// MIT License
//
// Copyright (c) 2026 Andrei Ishutin
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stack>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/utils.hpp>

/*! @file
 *  @brief Distinct-degree factorization (DDF).
 *
 *  DDF groups irreducible factors of a square-free polynomial by their
 *  degrees.
 *
 *  Input: square-free polynomial over a finite field.
 *  Output: products of irreducible factors grouped by degree.
 *
 *  This file contains several implementations:
 *    - naive: simple reference algorithm;
 *    - ntl_like: baby-step/giant-step algorithm close to NTL;
 *    - own_lazy: lazy giant-step variant;
 *    - own_tree: product-tree based variant.
 *
 *  To run a particular implementation, choose one of the namespaces listed
 *  above and call DistinctDegreeFactorize:
 *  @code
 *    auto factors =
 *        factorization::ddf::<implementation>::DistinctDegreeFactorize(poly);
 *  @endcode
 *
 *  DistinctDegreeFactorize functions accept a polynomial over a finite field,
 *  normalize it to monic form, and return an empty result for zero and constant
 *  polynomials.
 *
 *  DistinctDegreeFactorizer classes are lower-level helpers used by these
 *  wrappers. They expect the polynomial to already satisfy the DDF
 *  preconditions: monic and square-free.
 */

namespace factorization::ddf {

/*! @brief One output entry of DDF.
 *
 *  `factor` is the product of distinct irreducible polynomials of degree
 *  `degree`. Entries returned by DDF are pairwise coprime, and their product is
 *  the monic square-free input polynomial.
 */
template <concepts::Polynom Polynom>
struct DistinctDegreeFactor {
  //!< Product of irreducible factors of the same degree.
  Polynom factor;
  //!< Degree of every irreducible factor in factor.
  int degree;
};

namespace naive {

/*! @brief Straightforward reference implementation.
 *
 *  Assume q is field size. For each degree d it computes
 *    x^(q^d) (mod f)
 *  and extracts
 *    gcd(x^(q^d) - x, f).
 */
template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  using Element = typename Poly::Element;
  constexpr static auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

  std::vector<DistinctDegreeFactor<Poly>> result;

  // At the beginning of iteration `degree`,
  // h = x^(q^(degree - 1)) (mod poly).
  Poly h = x;
  size_t degree = 1;
  auto mod = poly.BuildModulus(2 * poly.Size());

  // At the beginning of iteration `degree`, poly contains only irreducible
  // factors of degree at least `degree`.
  while (2 * degree <= poly.Size() - 1) {
    h = polynomial::BinPowMod(std::move(h), kFieldSize, mod);
    Poly factor = poly.Gcd(h.Sub(x));
    // Extract all irreducible factors of degree `degree`.
    if (!factor.IsOne()) {
      result.emplace_back(factor, static_cast<int>(degree));
      poly = std::move(poly).Div(factor).MakeMonic();
      // factorization is complete
      if (poly.IsOne()) {
        break;
      }
      mod = poly.BuildModulus(2 * poly.Size());
      h = std::move(h).Rem(mod);
    }
    ++degree;
  }
  if (!poly.IsOne()) {
    const int factor_degree = static_cast<int>(poly.Size() - 1);
    result.emplace_back(std::move(poly), factor_degree);
  }
  return result;
}

}  // namespace naive

/*! @brief Strategy for choosing the baby-step block size.
 *
 *  Public DistinctDegreeFactorize wrappers currently use kSmallField:
 *  this is the better fit for the small finite fields.
 */
enum StepsMode {
  //!< Classical NTL choice: l ~= sqrt(n / 2).
  kExactNtl,
  //!< Larger baby-step blocks, useful when q-th powers (mod f) are cheap.
  kSmallField
};

template <concepts::Polynom Poly>
int Degree(const Poly& value) {
  return value.IsZero() ? 0 : static_cast<int>(value.Size()) - 1;
}

namespace ntl_like {

/*! @brief Practical baby-step/giant-step DDF in the style of NTL.
 *
 *  Notation used below:
 *    - n = deg(f), where f is the current square-free polynomial;
 *    - l = number of baby steps in one interval;
 *    - m = number of giant-step intervals;
 *    - h[i] = x^(q^i) (mod f) for 0 <= i <= l;
 *    - H[j] = x^(q^(j*l)) (mod f) for 1 <= j <= m;
 *    - F[j] stores factors with irreducible degrees in ((j - 1)l, jl].
 *
 *  This implementation follows the structure of NTL's `NewDDF`.
 *
 *  1. Generate baby steps:
 *         h[i] = x^(q^i) (mod f).
 *
 *  2. Generate giant steps:
 *         H[j] = x^(q^(j*l)) (mod f).
 *
 *  3. Giant refine:
 *         F[j] = gcd(f, product_{0 <= i < l}(H[j] - h[i]))
 *     for each interval 1 <= j <= m.
 *     At this point f no longer contains factors of degrees <= (j - 1)l, so
 *     F[j] collects factors with degrees in ((j - 1)l, jl].
 *
 *  4. Baby refine:
 *         gcd(F[j], H[j] - h[j*l - d]).
 *     This is the same kind of degree test as in the naive algorithm, but it
 *     is applied only inside a nontrivial interval F[j]. It splits F[j] into
 *     exact degree components by checking each degree d inside the interval.
 *
 *  Usage:
 *  @code
 *    DistinctDegreeFactorizer<Poly> factorizer(poly);
 *    auto factors = factorizer.Run();
 *  @endcode
 *
 *  @see https://www.shoup.net/papers/factorimpl.pdf
 *  @see https://libntl.org/doc/GF2EXFactoring.cpp.html
 */
template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  /*! @brief Batches interval products before computing gcd with poly.
   *
   *  Instead of computing gcd(poly, I[j]) for every interval separately, the
   *  buffer stores consecutive intervals and first computes
   *    gcd(poly, I[j] * I[j + 1] * ... * I[j + s - 1]).
   *  If the result is nontrivial, it is split back between the buffered
   *  intervals and written to F.
   */
  class GcdBuffer {
    // This value was obtained empirically.
    // NTL uses the same value.
    constexpr static int kBufferSize = 4;

   public:
    /*! @brief Creates an empty buffer for intervals of size interval_size. */
    explicit GcdBuffer(int interval_size)
        : interval_size_(interval_size) {
    }

    /*! @brief Returns true if there are no buffered interval products. */
    bool Empty() const {
      return size_ == 0;
    }

    /*! @brief Returns true if the buffer has reached its capacity. */
    bool Full() const {
      return size_ == kBufferSize;
    }

    /*! @brief Adds interval product for interval j.
     *
     *  Calls to Add must use consecutive interval numbers in increasing order.
     */
    void Add(int j, Poly factor) {
      assert(!Full());
      buf_[size_] = {j, std::move(factor)};
      ++size_;
    }

    /*! @brief Processes buffered products and clears the buffer.
     *
     *  Returns true if a nontrivial factor was extracted from poly.
     *  Extracted interval components are written to F.
     */
    bool Proceed(Poly& poly, std::vector<Poly>& F, const Modulus& mod) {
      if (Empty()) {
        return false;
      }
      // First test the product of all buffered intervals with one gcd.
      Poly product = buf_[0].second;
      for (int i = 1; i < size_; ++i) {
        product = std::move(product).Mul(buf_[i].second).Rem(mod);
      }
      product = std::move(product).Gcd(poly);
      if (product.IsOne()) {
        size_ = 0;
        return false;
      }

      poly = std::move(poly).Div(product).MakeMonic();
      int i = 0;
      int j = buf_[0].first;
      int min_degree = (j - 1) * interval_size_ + 1;
      // Split the nontrivial product back into interval slots while more than
      // one interval can still contribute.
      while (i < size_ - 1 && 2 * min_degree <= Degree(product)) {
        Poly factor = buf_[i].second.Gcd(product);
        if (!factor.IsOne()) {
          F[j] = std::move(factor);
          product = std::move(product).Div(F[j]).MakeMonic();
        }
        ++i;
        ++j;
        min_degree += interval_size_;
      }

      if (!product.IsOne()) {
        if (i != size_ - 1) {
          j = (Degree(product) + interval_size_ - 1) / interval_size_;
        }
        F[j] = std::move(product);
      }

      size_ = 0;
      return true;
    }

   private:
    int size_ = 0;
    int interval_size_;
    std::array<std::pair<int, Poly>, kBufferSize> buf_;
  };

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      // Classical NTL choice:
      // split degrees up to n / 2 into intervals of length
      //   l ~= sqrt(n / 2).
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      // Small-field heuristic: use larger baby-step blocks when raising to the
      // q-th power (mod f) is cheaper than modular composition.
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    // ceil(n / (2 * l)):
    m = (n + 2 * l - 1) / (2 * l);
  }

  /*! @brief Runs all stages and returns DDF components. */
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    auto mod = poly_.BuildModulus(2 * poly_.Size());
    GenerateBabySteps(mod);
    GenerateGiantSteps(mod);
    GiantRefine();
    BabyRefine();
    return std::move(result_);
  }

 private:
  /*! @brief Builds the baby-step table.
   *
   *  The table contains
   *    h[i] = x^(q^i) (mod f)
   *  for 0 <= i <= l.
   *
   *  All entries are computed modulo the same polynomial, so the precomputed
   *  Modulus can be reused.
   */
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    // The NTL choice advances the table by modular composition.
    // The small-field mode uses binary exponentiation to the q-th power.
    if constexpr (kMode == kExactNtl) {
      // Brent-Kung modular composition uses
      //   t ~= sqrt(n)
      // as its internal baby-step block size. This is unrelated to the
      // Kaltofen-Shoup baby-step table h.
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  /*! @brief Builds the giant-step table.
   *
   *  The table contains
   *    H[j] = x^(q^(j*l)) (mod f)
   *  for 1 <= j <= m.
   *
   *  All entries are computed modulo the same polynomial, so the precomputed
   *  Modulus can be reused.
   */
  void GenerateGiantSteps(const Modulus& mod) {
    // H[0] is never used.
    H.resize(m + 1);
    H[1] = h[l];

    // Brent-Kung modular composition uses
    //   t ~= sqrt(n)
    // as its internal baby-step block size. This is unrelated to the
    // Kaltofen-Shoup baby-step table h.
    int t = std::floor(std::sqrt(n));
    if (t == 0) {
      t = 1;
    }
    auto matrix = polynomial::BuildCompModMatrix(H[1], t, mod);
    for (int j = 2; j <= m; ++j) {
      // Advance by modular composition:
      //   H[j] = H[j - 1](H[1]) (mod f).
      H[j] = polynomial::CompMod(H[j - 1], matrix, mod);
    }
  }

  /*! @brief Fills F with nontrivial degree intervals.
   *
   *  After this stage, every nontrivial F[j] contains the product of all
   *  remaining irreducible factors with degrees in ((j - 1)l, jl].
   *  Extracted interval products are removed from poly_.
   *
   *  If a single large component remains after interval tests,
   *  it is appended directly to result_.
   *  After this stage Run no longer uses poly_.
   */
  void GiantRefine() {
    // F[0] is never used.
    F.assign(m + 1, Poly(Element::One()));
    // h was built modulo the original polynomial.
    // hh is reduced whenever the current polynomial poly_ shrinks.
    // h itself must stay unchanged because BabyRefine uses it later.
    std::vector<Poly> hh = h;

    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    GcdBuffer buffer(l);
    // At the beginning of interval j, poly_ contains only irreducible
    // factors with degrees greater than (j - 1)l.
    for (int j = 1; j <= m; ++j) {
      Poly HH = H[j].Rem(mod);  // NOLINT(readability-identifier-naming)
      // Interval product:
      //   I = product_{0 <= i < l}(H[j] - h[i]) (mod poly_).
      Poly I = HH.Sub(hh[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(HH.Sub(hh[i])).Rem(mod);
      }
      // The buffer batches several interval products and extracts their gcd
      // with poly_ in one call.
      // This happens at the very end or when buffer capacity is reached.
      buffer.Add(j, std::move(I));
      if (buffer.Full() && buffer.Proceed(poly_, F, mod)) {
        if (!poly_.IsOne() && Degree(poly_) >= 2 * j * l) {
          // poly_ has changed, so modular reductions must use a new modulus.
          mod = poly_.BuildModulus(2 * poly_.Size());
          for (int i = 1; i <= l; ++i) {
            hh[i] = std::move(hh[i]).Rem(mod);
          }
        }
      }
      // If
      //   deg(poly_) < 2jl,
      // all remaining factors of poly_ have degree > jl.
      // There can be at most one such factor, so interval tests are done.
      if (Degree(poly_) < 2 * j * l) {
        break;
      }
    }
    // Process the final, possibly non-full batch of interval products.
    buffer.Proceed(poly_, F, mod);
    // After all interval tests, any remaining nontrivial poly_ can only be one
    // component whose irreducible factors all have degree deg(poly_).
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
  }

  /*! @brief Splits all nontrivial intervals into exact degree components.
   *
   *  After this stage result_ also contains all components that were stored in
   *  F by GiantRefine.
   */
  void BabyRefine() {
    for (int j = 1; j < F.size(); ++j) {
      IntervalRefine(j, std::move(F[j]));
    }
  }

  /*! @brief Refines one interval product by testing individual degrees.
   *
   *  F_j contains factors with degrees in ((j - 1)l, jl]. This function tests
   *  those degrees one by one and appends exact DDF components to result_.
   */
  void IntervalRefine(int j,
                      Poly F_j) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F_j.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      // This loop visits degrees in increasing order:
      //   factor_degree = (j - 1)l + 1, ..., jl.
      const int factor_degree = j * l - i;
      if (Degree(F_j) < 2 * factor_degree) {
        // Any remaining nontrivial F_j can only be one component.
        if (!F_j.IsOne()) {
          const int d = Degree(F_j);
          result_.emplace_back(std::move(F_j), d);
        }
        break;
      }
      // Degree test:
      //   gcd(H[j] - h[i], F_j) = gcd((H[j] - h[i]) mod F_j, F_j).
      Poly factor = H[j].Sub(h[i]).Rem(mod).Gcd(F_j);
      if (!factor.IsOne()) {
        F_j = std::move(F_j).Div(factor).MakeMonic();
        mod = F_j.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> H;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> F;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace ntl_like

namespace own_lazy {

/*! @brief Lazy baby-step/giant-step DDF variant.
 *
 *  This variant keeps the baby-step table h, but does not materialize the
 *  giant-step table H. Instead, Run stores only the current giant step and
 *  advances it by modular composition after each interval.
 *
 *  Compared with ntl_like, this saves memory for H and can stop before
 *  computing future giant steps once the remaining polynomial becomes small.
 *  Nontrivial intervals are refined immediately.
 *
 *  Comments in this class focus on differences from ntl_like.
 *  This implementation is easier to read after the ntl_like variant.
 */
template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  /*! @brief Runs lazy DDF and returns exact degree components. */
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());

    GenerateBabySteps(mod);
    // Current giant step. Initially it is H[1].
    Poly H = h[l];  // NOLINT(readability-identifier-naming)

    // Brent-Kung modular composition uses
    //   t ~= sqrt(deg(poly_))
    // as its internal baby-step block size.
    // It is updated whenever poly_ shrinks.
    int t = std::floor(std::sqrt(n));
    t = std::max(t, 1);
    // matrix represents the fixed composition argument h[l] = x^(q^l).
    auto matrix = polynomial::BuildCompModMatrix(H, t, mod);

    // At the beginning of interval j:
    //   - H = H[j] (mod poly_),
    //   - h[i] = x^(q^i) (mod poly_) for 0 <= i <= l,
    //   - poly_ contains only irreducible factors with degrees > (j - 1)l,
    //   - t ~= sqrt(deg(poly_)),
    //   - matrix is the Brent-Kung table for h[l] modulo poly_ with block
    //     size t.
    //
    // Unlike ntl_like, the interval component F[j] is not stored.
    // It is computed, refined immediately, and removed from poly_.
    for (int j = 1; j <= m; ++j) {
      // Interval product:
      //   I = product_{0 <= i < l}(H - h[i]) (mod poly_).
      Poly I = H.Sub(h[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(H.Sub(h[i])).Rem(mod);
      }
      // Nontrivial F contains factors with degrees in ((j - 1)l, jl]:
      //   F = gcd(poly_, I).
      // This is the lazy counterpart of F[j] in the ntl_like variant.
      // NOLINTNEXTLINE(readability-identifier-naming)
      Poly F = std::move(I).Gcd(poly_);
      bool is_one = F.IsOne();
      if (!is_one) {
        poly_ = std::move(poly_).Div(F).MakeMonic();
        IntervalRefine(j, std::move(F), H);
      }
      // If
      //   deg(poly_) < 2jl,
      // all remaining factors of poly_ have degree > jl.
      // There can be at most one such factor, so interval tests are done.
      if (Degree(poly_) < 2 * j * l) {
        break;
      }
      if (!is_one) {
        // Restore the loop invariants after poly_ has changed.
        // Since deg(poly_) decreases, the modulus changes.
        // The Brent-Kung block size t is recomputed,
        // and the composition matrix must be updated accordingly.
        mod = poly_.BuildModulus(2 * poly_.Size());
        t = std::floor(std::sqrt(Degree(poly_)));
        t = std::max(t, 1);
        polynomial::UpdateCompModMatrix<Poly>(matrix, t, mod);
        H = std::move(H).Rem(mod);
        for (int i = 1; i <= l; ++i) {
          h[i] = std::move(h[i]).Rem(mod);
        }
      }
      // Advance to the next giant step:
      //   H[j + 1] = H[j](h[l]) (mod poly_).
      H = polynomial::CompMod(H, matrix, mod);
    }
    // Any remaining nontrivial poly_ is the single large component
    // left after interval tests.
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
    return std::move(result_);
  }

 private:
  /*! @brief Builds the same baby-step table as the NTL-like variant. */
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    if constexpr (kMode == kExactNtl) {
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  /*! @brief Refines one interval using the current giant step H.
   *
   *  Unlike ntl_like, H is passed as the current lazy value,
   *  not read from a precomputed giant-step table.
   */
  void IntervalRefine(int j, Poly F,    // NOLINT(readability-identifier-naming)
                      const Poly& H) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      const int factor_degree = j * l - i;
      if (Degree(F) < 2 * factor_degree) {
        if (!F.IsOne()) {
          const int d = Degree(F);
          result_.emplace_back(std::move(F), d);
        }
        break;
      }
      // Degree test:
      //   gcd(H - h[i], F) = gcd((H - h[i]) mod F, F).
      Poly factor = H.Sub(h[i]).Rem(mod).Gcd(F);
      if (!factor.IsOne()) {
        F = std::move(F).Div(factor).MakeMonic();
        mod = F.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace own_lazy

namespace own_tree {

/*! @brief Product-tree variant of baby-step/giant-step DDF.
 *
 *  This variant follows ntl_like except for GiantRefine. Instead of taking
 *  gcds interval by interval, it builds a product tree over interval products,
 *  computes one gcd at the root, and propagates the result back to leaves.
 *
 *  Comments in this class focus on product-tree specific logic. The shared
 *  baby-step, giant-step, and baby-refine stages are documented in ntl_like.
 */
template <concepts::Polynom Poly, StepsMode kMode = kSmallField>
class DistinctDegreeFactorizer {
  using Element = typename Poly::Element;
  using Modulus = typename Poly::Modulus;
  static constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  /*! @brief Product tree for interval products I[j].
   *
   *  The tree is one-shot: intervals are added, then
   *    - Build
   *    - Run
   *    - Extract
   *  are called once.
   *
   *  Run may move values out of non-leaf nodes while propagating backward
   *  values. Extract moves values out of leaves and finishes tree usage.
   *
   *  A node with range (from, to] represents consecutive intervals
   *  from + 1, ..., to.
   *
   *  Its forward value is
   *    I_{from + 1,to} = product_{from < j <= to} I[j] (mod f).
   *
   *  During the backward pass, the node receives
   *    F_{from + 1,to},
   *  the product of all irreducible factors whose degrees lie in
   *  (from*l, to*l].
   */
  class ComputationTree {
    struct Node {
      // I_{from + 1,to}: product of interval polynomials covered by this node.
      Poly forward;
      // F_{from + 1,to}: factor product assigned to this node.
      Poly backward;
      // Interval range represented by this node: (from, to].
      int from;
      int to;

      Node* left = nullptr;
      Node* right = nullptr;
    };

   public:
    /*! @brief Creates an empty one-shot tree for m intervals of size l. */
    ComputationTree(int interval_size, int interval_count, Modulus mod)
        : interval_size_(interval_size),
          nodes_(2 * interval_count),
          mod_(std::move(mod)) {
      for (auto& node : nodes_) {
        pool_.push(&node);
      }
    }

    // Adds a leaf with forward = I[interval].
    void Add(int interval, Poly poly) {
      Node* node = GetNode();
      leaves_.push_back(node);

      node->forward = std::move(poly);
      node->from = interval - 1;
      node->to = interval;
      node->left = nullptr;
      node->right = nullptr;
    }

    // Builds forward products bottom-up and returns I_{1,m}.
    Poly Build() {
      // Nodes on each level represent disjoint interval ranges.
      auto compare_nodes = [](Node* lhs, Node* rhs) {
        return lhs->from < rhs->from;
      };
      std::sort(leaves_.begin(), leaves_.end(), compare_nodes);
      std::vector<Node*> current_level = leaves_;
      while (current_level.size() > 1) {
        std::vector<Node*> next_level;
        int k = current_level.size() / 2;
        for (int i = 0; i < k; ++i) {
          Node* next = GetNode();
          Merge(current_level[2 * i], current_level[2 * i + 1], next);
          next_level.push_back(next);
        }
        if (current_level.size() % 2 != 0) {
          next_level.push_back(current_level.back());
        }
        current_level = std::move(next_level);
      }
      root_ = current_level[0];
      return std::move(root_->forward);
    }

    // Starts the backward pass with F_{1,m} stored at the root.
    // Non-leaf node values may be moved from during this pass.
    void Run(Poly root_val) {
      root_->backward = std::move(root_val);
      BackwardVisitNode(root_);
    }

    // Returns each interval j with the leaf value F[j].
    // Leaf values are moved out; the tree must not be used afterwards.
    std::vector<std::pair<int, Poly>> Extract() {
      std::vector<std::pair<int, Poly>> result;
      result.reserve(leaves_.size());
      for (Node* leave : leaves_) {
        result.emplace_back(leave->to, std::move(leave->backward));
      }
      return result;
    }

   private:
    void BackwardVisitNode(Node* node) {
      if (node->left == nullptr || node->right == nullptr) {
        return;
      }
      Split(node->left, node->right, node);
      BackwardVisitNode(node->left);
      BackwardVisitNode(node->right);
    }

    Node* GetNode() {
      Node* result = pool_.top();
      pool_.pop();
      return result;
    }

    void Merge(Node* left, Node* right, Node* root) {
      root->left = left;
      root->right = right;

      root->from = left->from;
      root->to = right->to;
      // right->forward is not needed during Split, so it can be moved.
      root->forward = std::move(right->forward).Mul(left->forward).Rem(mod_);
    }

    // Splits F_{root->from + 1,root->to} between child ranges using
    //   F_left = gcd(I_left, F_root),
    //   F_right = F_root / F_left.
    void Split(Node* left, Node* right, Node* root) {
      if (root->backward.IsOne()) {
        left->backward = root->backward;
        right->backward = root->backward;
        return;
      }
      int d = interval_size_ * root->from + 1;
      // If deg(F_root) is smaller than twice the minimal degree in this range,
      // F_root is one irreducible component and belongs to exactly one child.
      // root->backward is not needed after Split, so it can be moved.
      if (Degree(root->backward) < 2 * d) {
        int border = interval_size_ * right->from + 1;
        if (Degree(root->backward) < border) {
          left->backward = std::move(root->backward);
          right->backward = Poly(Element::One());
        } else {
          right->backward = std::move(root->backward);
          left->backward = Poly(Element::One());
        }
        return;
      }
      left->backward = root->backward.Gcd(std::move(left->forward));
      right->backward = root->backward.Div(left->backward);
    }

    const int interval_size_;
    std::vector<Node> nodes_;
    std::stack<Node*> pool_;
    std::vector<Node*> leaves_;

    const Modulus mod_;
    Node* root_;
  };

 public:
  explicit DistinctDegreeFactorizer(Poly poly)
      : poly_(std::move(poly)) {
    n = Degree(this->poly_);
    if constexpr (kMode == kExactNtl) {
      l = std::floor(std::sqrt(n / 2.0));
    } else {
      auto q = kFieldSize;
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(3.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  // One-shot
  std::vector<DistinctDegreeFactor<Poly>> Run() {
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    GenerateBabySteps(mod);
    GenerateGiantSteps(mod);
    GiantRefine();
    BabyRefine();
    return std::move(result_);
  }

 private:
  void GenerateBabySteps(const Modulus& mod) {
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

    h.resize(l + 1);
    h[0] = x;
    h[1] = polynomial::BinPowMod(x, kFieldSize, mod);

    if constexpr (kMode == kExactNtl) {
      int t = std::floor(std::sqrt(n));
      if (t == 0) {
        t = 1;
      }
      auto matrix = polynomial::BuildCompModMatrix(h[1], t, mod);
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::CompMod(h[i - 1], matrix, mod);
      }
    } else {
      for (size_t i = 2; i <= l; ++i) {
        h[i] = polynomial::BinPowMod(h[i - 1], kFieldSize, mod);
      }
    }
  }

  void GenerateGiantSteps(const Modulus& mod) {
    // H[0] is never used
    H.resize(m + 1);
    H[1] = h[l];

    int t = std::floor(std::sqrt(n));
    if (t == 0) {
      t = 1;
    }
    auto matrix = polynomial::BuildCompModMatrix(H[1], t, mod);
    for (int j = 2; j <= m; ++j) {
      // Advance by modular composition:
      //   H[j] = H[j - 1](H[1]) (mod f).
      H[j] = polynomial::CompMod(H[j - 1], matrix, mod);
    }
  }

  /*! @brief Fills F with nontrivial degree intervals using ComputationTree. */
  void GiantRefine() {
    // F[0] is never used.
    F.assign(m + 1, Poly(Element::One()));
    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    ComputationTree tree(l, m, mod);
    for (int j = 1; j <= m; ++j) {
      // Interval product:
      //   I[j] = product_{0 <= i < l}(H[j] - h[i]) (mod poly_).
      Poly I = H[j].Sub(h[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(H[j].Sub(h[i])).Rem(mod);
      }
      tree.Add(j, std::move(I));
    }
    Poly product = tree.Build();
    // product = I_{1,m}
    //   gcd(product, poly_)
    // gives F_{1,m}, the root value for the backward tree pass.
    product = std::move(product).Gcd(poly_);
    poly_ = std::move(poly_).Div(product);
    tree.Run(std::move(product));

    // Leaves now contain F[j].
    for (auto [j, factor] : tree.Extract()) {
      F[j] = std::move(factor);
    }
    // Any remaining nontrivial poly_ is the single large component
    // left after interval tests.
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
  }

  void BabyRefine() {
    for (int j = 1; j < F.size(); ++j) {
      IntervalRefine(j, std::move(F[j]));
    }
  }

  void IntervalRefine(int j,
                      Poly F_j) {  // NOLINT(readability-identifier-naming)
    Modulus mod = F_j.BuildModulus(n);
    for (int i = l; i-- > 0;) {
      const int factor_degree = j * l - i;
      if (Degree(F_j) < 2 * factor_degree) {
        if (!F_j.IsOne()) {
          const int d = Degree(F_j);
          result_.emplace_back(std::move(F_j), d);
        }
        break;
      }
      // factor = Gcd(H_j - h_i, F_j)
      // but since deg F_j is low we use rem first
      Poly factor = H[j].Sub(h[i]).Rem(mod).Gcd(F_j);
      if (!factor.IsOne()) {
        F_j = std::move(F_j).Div(factor).MakeMonic();
        mod = F_j.BuildModulus(n);
        result_.emplace_back(std::move(factor), factor_degree);
      }
    }
  }

  Poly poly_;
  int n = 0;            // NOLINT(readability-identifier-naming)
  int l = 1;            // NOLINT(readability-identifier-naming)
  int m = 1;            // NOLINT(readability-identifier-naming)
  std::vector<Poly> h;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> H;  // NOLINT(readability-identifier-naming)
  std::vector<Poly> F;  // NOLINT(readability-identifier-naming)
  std::vector<DistinctDegreeFactor<Poly>> result_;
};

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const size_t n = Degree(poly);
  if (n == 1) {
    return {{std::move(poly), 1}};
  }
  return DistinctDegreeFactorizer<Poly>(std::move(poly)).Run();
}

}  // namespace own_tree

}  // namespace factorization::ddf
