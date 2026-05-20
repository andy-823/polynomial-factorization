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
    return PlainGCD(std::move(a), std::move(b));
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
};

}  // namespace factorization::polynomial
