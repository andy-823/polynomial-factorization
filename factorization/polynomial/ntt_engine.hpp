// SPDX-License-Identifier: LGPL-2.1-or-later
//
// Modifications Copyright (c) 2026 Andrei Ishutin
//
// Portions of the Half-GCD/GCD implementation are adapted from NTL's ZZ_pEX
// implementation. NTL is written and maintained by Victor Shoup and distributed
// under the GNU Lesser General Public License version 2.1 or later:
//   https://libntl.org

#pragma once

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>

namespace factorization::polynomial {

namespace detail {

// This implementation is taken from
//   https://codeforces.com/blog/entry/129600?locale=ru
template <uint64_t kMod, uint64_t kGenerator>
class IntegerNtt {
 public:
  [[nodiscard]]
  static std::vector<uint64_t> Convolve(std::vector<uint64_t> first,
                                        std::vector<uint64_t> second,
                                        size_t result_size) {
    const size_t ntt_size = NttSize(result_size);
    Montgomery red(kMod);

    Prepare(first, ntt_size, red);
    Prepare(second, ntt_size, red);

    const uint64_t inv_n = red.Transform(BinPow(ntt_size, kMod - 2, kMod));
    const uint64_t root =
        red.Transform(BinPow(kGenerator, (kMod - 1) / ntt_size, kMod));
    const uint64_t inv_root =
        red.Transform(BinPow(red.Reduce(root), kMod - 2, kMod));
    const auto& rev = BitSort(ntt_size);

    Ntt(first, rev, red, inv_n, root, inv_root, false);
    Ntt(second, rev, red, inv_n, root, inv_root, false);
    for (size_t i = 0; i < ntt_size; ++i) {
      first[i] = red.Multiply(first[i], second[i]);
    }
    Ntt(first, rev, red, inv_n, root, inv_root, true);

    first.resize(result_size);
    for (auto& value : first) {
      value = red.Reduce(value);
    }
    return first;
  }

 private:
  struct Montgomery {
    uint64_t n;
    uint64_t nr;

    constexpr explicit Montgomery(uint64_t value)
        : n(value),
          nr(1) {
      for (int i = 0; i < 6; ++i) {
        nr *= 2 - n * nr;
      }
    }

    [[nodiscard]]
    uint64_t Reduce(__uint128_t value) const {
      const uint64_t q = static_cast<uint64_t>(value) * nr;
      const uint64_t m = (static_cast<__uint128_t>(q) * n) >> 64;
      uint64_t result = (value >> 64) + n - m;
      if (result >= n) {
        result -= n;
      }
      return result;
    }

    [[nodiscard]]
    uint64_t Multiply(uint64_t first, uint64_t second) const {
      return Reduce(static_cast<__uint128_t>(first) * second);
    }

    [[nodiscard]]
    uint64_t Transform(uint64_t value) const {
      return (static_cast<__uint128_t>(value) << 64) % n;
    }
  };

  [[nodiscard]]
  static uint64_t BinPow(uint64_t base, uint64_t power, uint64_t mod) {
    uint64_t result = 1;
    while (power > 0) {
      if ((power & 1) != 0) {
        result = (static_cast<__uint128_t>(result) * base) % mod;
      }
      base = (static_cast<__uint128_t>(base) * base) % mod;
      power >>= 1;
    }
    return result;
  }

  [[nodiscard]]
  static size_t NttSize(size_t result_size) {
    size_t ntt_size = 1;
    while (ntt_size < result_size) {
      ntt_size <<= 1;
    }
    return ntt_size;
  }

  static void Prepare(std::vector<uint64_t>& values, size_t ntt_size,
                      const Montgomery& red) {
    const size_t old_size = values.size();
    values.resize(ntt_size);
    for (size_t i = 0; i < old_size; ++i) {
      values[i] = red.Transform(values[i]);
    }
  }

  [[nodiscard]]
  static const std::vector<int>& BitSort(size_t ntt_size) {
    static std::vector<int> rev;
    if (rev.size() == ntt_size) {
      return rev;
    }

    rev.assign(ntt_size, 0);
    int log = 0;
    while ((size_t{1} << log) < ntt_size) {
      ++log;
    }
    for (size_t i = 1; i < ntt_size; ++i) {
      rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (log - 1));
    }
    return rev;
  }

  static void Ntt(std::vector<uint64_t>& values, const std::vector<int>& rev,
                  const Montgomery& red, uint64_t inv_n, uint64_t root,
                  uint64_t inv_root, bool invert) {
    const auto size = static_cast<int>(values.size());

    for (int i = 0; i < size; ++i) {
      if (i < rev[i]) {
        std::swap(values[i], values[rev[i]]);
      }
    }

    const uint64_t step = invert ? inv_root : root;
    std::vector<uint64_t> roots(size >> 1);
    if (!roots.empty()) {
      roots[0] = red.Transform(1);
    }
    for (int i = 1; i < (size >> 1); ++i) {
      roots[i] = red.Multiply(roots[i - 1], step);
    }

    int log = 0;
    while ((1 << log) < size) {
      ++log;
    }
    for (int level = 0; level < log; ++level) {
      for (int i = 0; i < size; ++i) {
        if ((i & (1 << level)) != 0) {
          continue;
        }
        const int pair = i ^ (1 << level);
        const int root_index = (i & ((1 << level) - 1)) * (size >> (level + 1));
        const uint64_t value = red.Multiply(values[pair], roots[root_index]);
        values[pair] =
            values[i] >= value ? values[i] - value : values[i] + kMod - value;
        values[i] = values[i] + value < kMod ? values[i] + value
                                             : values[i] + value - kMod;
      }
    }

    if (invert) {
      for (auto& value : values) {
        value = red.Multiply(value, inv_n);
      }
    }
  }
};

}  // namespace detail

template <concepts::GaloisFieldElement Elem>
struct NttEngine {
  static_assert(Elem::FieldPower() == 1, "NttEngine requires prime field");

  // The Half-GCD routines in this engine follow the structure of NTL's ZZ_pEX
  // HalfGCD/GCD implementation, adapted to this polynomial representation.
  constexpr static size_t kHalfGCDThreshold = 512;

  struct HalfGCDMatrix {
    std::vector<Elem> data[2][2];
  };

  struct Modulus {
    std::vector<Elem> polynomial;
    std::vector<Elem> reversed_inverse;
    size_t max_quotient_size = 0;
  };

  [[nodiscard]]
  static std::vector<Elem> Mul(const std::vector<Elem>& a,
                               const std::vector<Elem>& b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    const size_t result_size = a.size() + b.size() - 1;
    return MulNtt(a, b, result_size);
  }

  // assume a.size() >= b.size()
  [[nodiscard]]
  static std::vector<Elem> Rem(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
    if (ShouldUsePlainDiv(a, b)) {
      return PlainRem(std::move(a), b);
    }
    auto quotient = Div(a, b);
    return Sub(std::move(a), Mul(std::move(quotient), b), b.size() - 1);
  }

  [[nodiscard]]
  static std::vector<Elem> Rem(std::vector<Elem> a, const Modulus& modulus) {
    if (a.size() < modulus.polynomial.size()) {
      return Trim(std::move(a));
    }
    if (ShouldUsePlainDiv(a, modulus.polynomial)) {
      return PlainRem(std::move(a), modulus.polynomial);
    }
    auto quotient = Div(a, modulus);
    return Sub(std::move(a), Mul(std::move(quotient), modulus.polynomial),
               modulus.polynomial.size() - 1);
  }

  // assume a.size() >= b.size()
  [[nodiscard]]
  static std::vector<Elem> Div(const std::vector<Elem>& a,
                               const std::vector<Elem>& b) {
    if (ShouldUsePlainDiv(a, b)) {
      return PlainDiv(a, b);
    }
    const size_t quotient_size = a.size() - b.size() + 1;
    std::vector<Elem> rev_a = ReverseTake(a, quotient_size);
    std::vector<Elem> rev_b = ReverseTake(b, quotient_size);
    auto inv = InverseMod(rev_b, quotient_size);
    auto quotient = Mul(std::move(rev_a), inv);
    quotient.resize(quotient_size);
    std::reverse(quotient.begin(), quotient.end());
    return quotient;
  }

  [[nodiscard]]
  static std::vector<Elem> Div(const std::vector<Elem>& a,
                               const Modulus& modulus) {
    if (a.size() < modulus.polynomial.size()) {
      return {};
    }
    if (ShouldUsePlainDiv(a, modulus.polynomial)) {
      return PlainDiv(a, modulus.polynomial);
    }
    const size_t quotient_size = a.size() - modulus.polynomial.size() + 1;
    if (quotient_size > modulus.max_quotient_size) {
      return Div(a, modulus.polynomial);
    }

    std::vector<Elem> rev_a = ReverseTake(a, quotient_size);
    std::vector<Elem> inv(
        modulus.reversed_inverse.begin(),
        modulus.reversed_inverse.begin() +
            std::min(modulus.reversed_inverse.size(), quotient_size));
    auto quotient = Mul(std::move(rev_a), inv);
    quotient.resize(quotient_size);
    std::reverse(quotient.begin(), quotient.end());
    return quotient;
  }

  [[nodiscard]]
  static std::pair<std::vector<Elem>, std::vector<Elem>> DivRem(
      std::vector<Elem> a, const std::vector<Elem>& b) {
    if (a.size() < b.size()) {
      return {{}, Trim(std::move(a))};
    }
    if (ShouldUsePlainDiv(a, b)) {
      return PlainDivRem(std::move(a), b);
    }
    auto quotient = Div(a, b);
    auto remainder = Sub(std::move(a), Mul(quotient, b), b.size() - 1);
    return {std::move(quotient), std::move(remainder)};
  }

  [[nodiscard]]
  static std::pair<std::vector<Elem>, std::vector<Elem>> DivRem(
      std::vector<Elem> a, const Modulus& modulus) {
    if (a.size() < modulus.polynomial.size()) {
      return {{}, Trim(std::move(a))};
    }
    if (ShouldUsePlainDiv(a, modulus.polynomial)) {
      return PlainDivRem(std::move(a), modulus.polynomial);
    }
    auto quotient = Div(a, modulus);
    auto remainder = Sub(std::move(a), Mul(quotient, modulus.polynomial),
                         modulus.polynomial.size() - 1);
    return {std::move(quotient), std::move(remainder)};
  }

  [[nodiscard]]
  static Modulus BuildModulus(const std::vector<Elem>& polynomial,
                              size_t max_dividend_size) {
    if (polynomial.empty()) {
      return {};
    }
    if (max_dividend_size < polynomial.size()) {
      return {polynomial, {}, 0};
    }
    const size_t quotient_size = max_dividend_size - polynomial.size() + 1;
    std::vector<Elem> rev_polynomial = ReverseTake(polynomial, quotient_size);
    return {
        polynomial,
        InverseMod(rev_polynomial, quotient_size),
        quotient_size,
    };
  }

  [[nodiscard]]
  static std::vector<Elem> Gcd(std::vector<Elem> a, std::vector<Elem> b) {
    return HalfGCDImpl(std::move(a), std::move(b));
    // return PlainGCD(std::move(a), std::move(b));
  }

 private:
  constexpr static uint64_t kNttMod = 2524775926340780033;
  constexpr static uint64_t kNttGenerator = 3;

  constexpr static size_t kPlainDivThreshold = 128;

  static void TrimInPlace(std::vector<Elem>& a) {
    while (!a.empty() && a.back() == Elem::Zero()) {
      a.pop_back();
    }
  }

  [[nodiscard]]
  static std::vector<Elem> Trim(std::vector<Elem> a) {
    TrimInPlace(a);
    return a;
  }

  [[nodiscard]]
  static std::vector<Elem> Sub(std::vector<Elem> a, const std::vector<Elem>& b,
                               size_t max_size = 0) {
    const size_t result_size = std::max(a.size(), b.size());
    if (a.size() < result_size) {
      a.resize(result_size, Elem::Zero());
    }
    for (size_t i = 0; i < b.size(); ++i) {
      a[i] -= b[i];
    }
    if (max_size != 0 && a.size() > max_size) {
      a.resize(max_size);
    }
    return Trim(std::move(a));
  }

  [[nodiscard]]
  static std::vector<Elem> ReverseTake(const std::vector<Elem>& a,
                                       size_t size) {
    std::vector<Elem> result;
    const size_t result_size = std::min(a.size(), size);
    result.reserve(result_size);
    for (size_t i = 0; i < result_size; ++i) {
      result.push_back(a[a.size() - 1 - i]);
    }
    return Trim(std::move(result));
  }

  [[nodiscard]]
  static bool ShouldUsePlainDiv(const std::vector<Elem>& a,
                                const std::vector<Elem>& b) {
    return b.size() <= kPlainDivThreshold ||
           a.size() - b.size() + 1 <= kPlainDivThreshold;
  }

  [[nodiscard]]
  static std::vector<Elem> PlainRem(std::vector<Elem> a,
                                    const std::vector<Elem>& b) {
    const size_t divisor_size = b.size();
    const size_t quotient_size = a.size() - divisor_size + 1;
    const Elem lead_inverse = b.back().Inverse();

    for (size_t i = quotient_size; i-- > 0;) {
      Elem coeff = a[i + divisor_size - 1] * lead_inverse;
      if (coeff == Elem::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j + 1 < divisor_size; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    a.resize(divisor_size - 1);
    return Trim(std::move(a));
  }

  [[nodiscard]]
  static std::vector<Elem> PlainDiv(std::vector<Elem> a,
                                    const std::vector<Elem>& b) {
    return PlainDivRem(std::move(a), b).first;
  }

  [[nodiscard]]
  static std::pair<std::vector<Elem>, std::vector<Elem>> PlainDivRem(
      std::vector<Elem> a, const std::vector<Elem>& b) {
    const size_t divisor_size = b.size();
    const size_t quotient_size = a.size() - divisor_size + 1;
    const Elem lead_inverse = b.back().Inverse();
    std::vector<Elem> quotient(quotient_size, Elem::Zero());

    for (size_t i = quotient_size; i-- > 0;) {
      Elem coeff = a[i + divisor_size - 1] * lead_inverse;
      quotient[i] = coeff;
      if (coeff == Elem::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j + 1 < divisor_size; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    a.resize(divisor_size - 1);
    return {Trim(std::move(quotient)), Trim(std::move(a))};
  }

  [[nodiscard]]
  static std::vector<Elem> MulTrunc(std::vector<Elem> a,
                                    const std::vector<Elem>& b, size_t size) {
    auto result = Mul(std::move(a), b);
    if (result.size() > size) {
      result.resize(size);
    }
    return result;
  }

  [[nodiscard]]
  static std::vector<Elem> InverseMod(const std::vector<Elem>& a, size_t size) {
    if (size == 1) {
      return {a[0].Inverse()};
    }

    const size_t k = (size + 1) / 2;
    std::vector<Elem> a1(a.begin(), a.begin() + std::min(a.size(), k));
    auto b1 = InverseMod(a1, k);

    auto product = MulTrunc(a, b1, size);
    std::vector<Elem> c(size - k, Elem::Zero());
    if (product.size() > k) {
      const size_t c_size = std::min(product.size() - k, c.size());
      std::copy(product.begin() + k, product.begin() + k + c_size, c.begin());
      TrimInPlace(c);
    }

    auto b2 = MulTrunc(b1, c, size - k);
    for (auto& value : b2) {
      value = -value;
    }

    std::vector<Elem> g(std::move(b1));
    g.resize(k, Elem::Zero());
    g.insert(g.end(), b2.begin(), b2.end());
    if (g.size() > size) {
      g.resize(size);
    }
    return Trim(std::move(g));
  }

  [[nodiscard]]
  static std::vector<Elem> PlainGCD(std::vector<Elem> a, std::vector<Elem> b) {
    if (a.size() < b.size()) {
      a.swap(b);
    }
    while (!b.empty()) {
      a = Rem(std::move(a), b);
      a.swap(b);
    }
    return a;
  }

  [[nodiscard]]
  static std::vector<Elem> MulNtt(const std::vector<Elem>& a,
                                  const std::vector<Elem>& b,
                                  size_t result_size) {
    std::vector<uint64_t> first;
    first.reserve(a.size());
    for (const auto& value : a) {
      first.push_back(value.Get()[0]);
    }

    std::vector<uint64_t> second;
    second.reserve(b.size());
    for (const auto& value : b) {
      second.push_back(value.Get()[0]);
    }

    auto convolution = detail::IntegerNtt<kNttMod, kNttGenerator>::Convolve(
        std::move(first), std::move(second), result_size);

    std::vector<Elem> result;
    result.reserve(result_size);
    const uint64_t field_base = Elem::FieldBase();
    for (const auto value : convolution) {
      result.emplace_back(
          Elem(static_cast<typename Elem::Coefficient>(value % field_base)));
    }
    while (!result.empty() && result.back() == Elem::Zero()) {
      result.pop_back();
    }
    return result;
  }

  [[nodiscard]]
  static std::vector<Elem> Add(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
    if (a.size() < b.size()) {
      a.resize(b.size(), Elem::Zero());
    }
    for (size_t i = 0; i < b.size(); ++i) {
      a[i] += b[i];
    }
    return Trim(std::move(a));
  }

  [[nodiscard]]
  static std::vector<Elem> HalfGCDImpl(std::vector<Elem> a,
                                       std::vector<Elem> b) {
    if (a.size() == b.size()) {
      if (a.empty()) {
        return {};
      }
      b = Rem(std::move(b), a);
    } else if (a.size() < b.size()) {
      a.swap(b);
    }

    while (a.size() >= kHalfGCDThreshold && !b.empty()) {
      HalfGCDReduce(a, b);
      if (!b.empty()) {
        a = Rem(std::move(a), b);
        a.swap(b);
      }
    }
    return PlainGCD(std::move(a), std::move(b));
  }

  static void HalfGCDReduce(std::vector<Elem>& u, std::vector<Elem>& v) {
    const size_t d_red = u.size() / 2;
    if (v.empty() || v.size() <= u.size() - d_red) {
      return;
    }

    const size_t du = u.size() - 1;
    size_t d1 = (d_red + 1) / 2;
    d1 = std::max<size_t>(d1, 1);
    if (d1 >= d_red) {
      d1 = d_red - 1;
    }

    HalfGCDMatrix m1;
    HalfGCD(m1, u, v, d1);
    ApplyMatrix(u, v, m1);

    const long d2 =
        Degree(v) - static_cast<long>(du) + static_cast<long>(d_red);
    if (v.empty() || d2 <= 0) {
      return;
    }

    auto [quotient, remainder] = DivRem(u, v);
    u = std::move(remainder);
    u.swap(v);

    HalfGCD(m1, u, v, static_cast<size_t>(d2));
    ApplyMatrix(u, v, m1);
  }

  static void HalfGCD(HalfGCDMatrix& result, const std::vector<Elem>& u,
                      const std::vector<Elem>& v, size_t d_red) {
    if (v.empty() || v.size() <= u.size() - d_red) {
      SetIdentity(result);
      return;
    }

    const long shift = Degree(u) - 2 * static_cast<long>(d_red) + 2;
    auto u1 = RightShift(u, shift > 0 ? static_cast<size_t>(shift) : 0);
    auto v1 = RightShift(v, shift > 0 ? static_cast<size_t>(shift) : 0);

    if (d_red <= kHalfGCDThreshold) {
      IterHalfGCD(result, u1, v1, d_red);
      return;
    }

    size_t d1 = (d_red + 1) / 2;
    d1 = std::max<size_t>(d1, 1);
    if (d1 >= d_red) {
      d1 = d_red - 1;
    }

    HalfGCDMatrix m1;
    HalfGCD(m1, u1, v1, d1);
    ApplyMatrix(u1, v1, m1);

    const long d2 = Degree(v1) - Degree(u) + shift + static_cast<long>(d_red);
    if (v1.empty() || d2 <= 0) {
      result = std::move(m1);
      return;
    }

    auto [quotient, remainder] = DivRem(u1, v1);
    u1 = std::move(remainder);
    u1.swap(v1);

    HalfGCDMatrix m2;
    HalfGCD(m2, u1, v1, static_cast<size_t>(d2));
    UpdateAfterDivision(m1, quotient);
    result = MulMatrix(std::move(m2), std::move(m1));
  }

  static void IterHalfGCD(HalfGCDMatrix& result, std::vector<Elem>& u,
                          std::vector<Elem>& v, size_t d_red) {
    SetIdentity(result);
    const size_t goal_size = u.size() - d_red;
    if (v.size() <= goal_size) {
      return;
    }

    while (v.size() > goal_size) {
      auto [quotient, remainder] = DivRem(u, v);
      u = std::move(remainder);
      u.swap(v);
      UpdateAfterDivision(result, quotient);
    }
  }

  static void SetIdentity(HalfGCDMatrix& matrix) {
    matrix.data[0][0] = {Elem::One()};
    matrix.data[0][1] = {};
    matrix.data[1][0] = {};
    matrix.data[1][1] = {Elem::One()};
  }

  static void ApplyMatrix(std::vector<Elem>& u, std::vector<Elem>& v,
                          const HalfGCDMatrix& matrix) {
    auto new_u = AddProducts(matrix.data[0][0], u, matrix.data[0][1], v);
    auto new_v = AddProducts(matrix.data[1][0], u, matrix.data[1][1], v);
    u = std::move(new_u);
    v = std::move(new_v);
  }

  [[nodiscard]]
  static HalfGCDMatrix MulMatrix(HalfGCDMatrix lhs, HalfGCDMatrix rhs) {
    return {{
        {AddProducts(lhs.data[0][0], rhs.data[0][0], lhs.data[0][1],
                     rhs.data[1][0]),
         AddProducts(lhs.data[0][0], rhs.data[0][1], lhs.data[0][1],
                     rhs.data[1][1])},
        {AddProducts(lhs.data[1][0], rhs.data[0][0], lhs.data[1][1],
                     rhs.data[1][0]),
         AddProducts(lhs.data[1][0], rhs.data[0][1], lhs.data[1][1],
                     rhs.data[1][1])},
    }};
  }

  [[nodiscard]]
  static std::vector<Elem> MulTerm(const std::vector<Elem>& lhs,
                                   const std::vector<Elem>& rhs) {
    if (lhs.empty() || rhs.empty()) {
      return {};
    }
    if (lhs.size() == 1 && lhs[0] == Elem::One()) {
      return rhs;
    }
    if (rhs.size() == 1 && rhs[0] == Elem::One()) {
      return lhs;
    }
    return Mul(lhs, rhs);
  }

  [[nodiscard]]
  static std::vector<Elem> AddProducts(const std::vector<Elem>& a,
                                       const std::vector<Elem>& b,
                                       const std::vector<Elem>& c,
                                       const std::vector<Elem>& d) {
    auto first = MulTerm(a, b);
    if (first.empty()) {
      return MulTerm(c, d);
    }
    return Add(std::move(first), MulTerm(c, d));
  }

  static void UpdateAfterDivision(HalfGCDMatrix& matrix,
                                  const std::vector<Elem>& quotient) {
    auto top_left = std::move(matrix.data[1][0]);
    auto bottom_left =
        Sub(std::move(matrix.data[0][0]), MulByQuotient(quotient, top_left));
    matrix.data[0][0] = std::move(top_left);
    matrix.data[1][0] = std::move(bottom_left);

    auto top_right = std::move(matrix.data[1][1]);
    auto bottom_right =
        Sub(std::move(matrix.data[0][1]), MulByQuotient(quotient, top_right));
    matrix.data[0][1] = std::move(top_right);
    matrix.data[1][1] = std::move(bottom_right);
  }

  [[nodiscard]]
  static std::vector<Elem> MulByQuotient(const std::vector<Elem>& quotient,
                                         const std::vector<Elem>& value) {
    if (quotient.size() > 4) {
      return MulTerm(quotient, value);
    }
    if (quotient.empty() || value.empty()) {
      return {};
    }
    if (quotient.size() == 1) {
      if (quotient[0] == Elem::Zero()) {
        return {};
      }
      if (quotient[0] == Elem::One()) {
        return value;
      }
      std::vector<Elem> result(value);
      for (auto& coeff : result) {
        coeff *= quotient[0];
      }
      return Trim(std::move(result));
    }

    std::vector<Elem> result(value.size() + quotient.size() - 1, Elem::Zero());
    for (size_t i = 0; i < quotient.size(); ++i) {
      if (quotient[i] == Elem::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < value.size(); ++j) {
        if (value[j] == Elem::Zero()) [[unlikely]] {
          continue;
        }
        result[i + j] += quotient[i] * value[j];
      }
    }
    return Trim(std::move(result));
  }

  [[nodiscard]]
  static std::vector<Elem> RightShift(const std::vector<Elem>& a,
                                      size_t shift) {
    if (shift >= a.size()) {
      return {};
    }
    return Trim(std::vector<Elem>(a.begin() + shift, a.end()));
  }

  [[nodiscard]]
  static long Degree(const std::vector<Elem>& a) {
    if (a.empty()) {
      return -1;
    }
    return static_cast<long>(a.size() - 1);
  }
};

}  // namespace factorization::polynomial
