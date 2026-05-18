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

#include <cmath>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/utils.hpp>

namespace factorization::ddf {

template <concepts::Polynom Polynom>
struct DistinctDegreeFactor {
  Polynom factor;
  int degree;
};

namespace naive {

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

  Poly h = x;
  size_t degree = 1;
  auto mod = poly.BuildModulus(2 * poly.Size());

  while (2 * degree <= poly.Size() - 1) {
    h = polynomial::BinPowMod(std::move(h), kFieldSize, mod);
    Poly factor = poly.Gcd(h.Sub(x));
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

namespace ntl_like {

template <concepts::Polynom Poly>
int Degree(const Poly& value) {
  return value.IsZero() ? 0 : static_cast<int>(value.Size()) - 1;
}

enum Mode { kExactNtl, kSmallField };

template <concepts::Polynom Poly, Mode kMode = kSmallField>
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
      l = std::floor(std::pow(n, 0.75L) / std::sqrt(2.0L * std::log2(q)));
    }
    if (l == 0) {
      l = 1;
    }
    m = (n + 2 * l - 1) / (2 * l);
  }

  // One-shot
  std::vector<DistinctDegreeFactor<Poly>> Run() {

    return RunWithObserver([](const char*, bool) {});
  }

  template <typename Observer>
  std::vector<DistinctDegreeFactor<Poly>> RunWithObserver(Observer&& observer) {
    observer("BuildModulus", true);
    auto mod = poly_.BuildModulus(2 * poly_.Size());
    observer("BuildModulus", false);
    observer("GenerateBabySteps", true);
    GenerateBabySteps(mod);
    observer("GenerateBabySteps", false);
    observer("GenerateGiantSteps", true);
    GenerateGiantSteps(mod);
    observer("GenerateGiantSteps", false);
    observer("GiantRefine", true);
    GiantRefine();
    observer("GiantRefine", false);
    observer("BabyRefine", true);
    BabyRefine();
    observer("BabyRefine", false);
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
      // H_j = H_{j-1} (H_1) (mod f)
      H[j] = polynomial::CompMod(H[j - 1], matrix, mod);
    }
  }

  void GiantRefine() {
    // F[0] is never used btw
    F.assign(m + 1, Poly(Element::One()));
    std::vector<Poly> hh = h;

    Modulus mod = poly_.BuildModulus(2 * poly_.Size());
    for (int j = 1; j <= m; ++j) {
      Poly HH = H[j].Rem(mod);  // NOLINT(readability-identifier-naming)
      // I = (H_j - h_0) * ... * (H_j - h_{l-1})
      Poly I = HH.Sub(hh[0]);  // NOLINT(readability-identifier-naming)
      for (int i = 1; i < l; ++i) {
        I = std::move(I).Mul(HH.Sub(hh[i])).Rem(mod);
      }
      F[j] = std::move(I).Gcd(poly_);
      if (!F[j].IsOne()) {
        poly_ = std::move(poly_).Div(F[j]).MakeMonic();
      }
      if (Degree(poly_) < 2 * j * l) {
        break;
      }
      if (!F[j].IsOne()) {
        mod = poly_.BuildModulus(2 * poly_.Size());
        for (int i = 1; i <= l; ++i) {
          hh[i] = std::move(hh[i]).Rem(mod);
        }
      }
    }
    if (!poly_.IsOne()) {
      const int d = Degree(poly_);
      result_.emplace_back(std::move(poly_), d);
    }
  }

  void BabyRefine() {
    for (int j = 1; j < F.size(); ++j) {
      Modulus mod = F[j].BuildModulus(n);
      for (int i = l; i-- > 0;) {
        const int factor_degree = j * l - i;
        if (Degree(F[j]) < 2 * factor_degree) {
          if (!F[j].IsOne()) {
            const int d = Degree(F[j]);
            result_.emplace_back(std::move(F[j]), d);
          }
          break;
        }
        // factor = Gcd(H_j - h_i, F_j)
        // bu since deg F_j is low we use rem at first
        Poly factor = H[j].Sub(h[i]).Rem(mod).Gcd(F[j]);
        if (!factor.IsOne()) {
          F[j] = std::move(F[j]).Div(factor).MakeMonic();
          mod = F[j].BuildModulus(n);
          result_.emplace_back(std::move(factor), factor_degree);
        }
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

}  // namespace factorization::ddf
