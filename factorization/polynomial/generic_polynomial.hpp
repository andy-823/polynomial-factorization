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
#include <span>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>

#include "polynomial_base.hpp"

namespace factorization::polynomial {

template <concepts::GaloisFieldElement Elem>
struct PolynomialEngine {
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
    if (a.size() < b.size()) {
      return MulKaratsuba(b, a);
    }
    return MulKaratsuba(a, b);
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
  constexpr static size_t kKaratsubaThreshold = 128;
  constexpr static size_t kPlainDivThreshold = 128;
  constexpr static size_t kHalfGCDThreshold = 256;

  struct HalfGCDMatrix {
    std::vector<Elem> data[2][2];
  };

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
  static std::vector<Elem> PlainMul(std::span<const Elem> a,
                                    std::span<const Elem> b) {
    std::vector<Elem> result(a.size() + b.size() - 1, Elem::Zero());
    for (size_t i = 0; i < a.size(); ++i) {
      if (a[i] == Elem::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < b.size(); ++j) {
        result[i + j] += a[i] * b[j];
      }
    }
    return Trim(std::move(result));
  }

  [[nodiscard]]
  static std::vector<Elem> AddParts(std::span<const Elem> a,
                                    std::span<const Elem> b) {
    std::vector<Elem> result(a.begin(), a.end());
    if (result.size() < b.size()) {
      result.resize(b.size(), Elem::Zero());
    }
    for (size_t i = 0; i < b.size(); ++i) {
      result[i] += b[i];
    }
    return Trim(std::move(result));
  }

  [[nodiscard]]
  static std::vector<Elem> MulKaratsuba(std::span<const Elem> a,
                                        std::span<const Elem> b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    if (b.size() == 1) {
      if (b[0] == Elem::Zero()) {
        return {};
      }
      std::vector<Elem> result(a.begin(), a.end());
      if (b[0] != Elem::One()) {
        for (auto& value : result) {
          value *= b[0];
        }
      }
      return Trim(std::move(result));
    }
    if (std::min(a.size(), b.size()) <= kKaratsubaThreshold) {
      return PlainMul(a, b);
    }

    const size_t split = std::min(a.size(), b.size()) / 2;
    auto a_low = a.first(split);
    auto a_high = a.subspan(split);
    auto b_low = b.first(split);
    auto b_high = b.subspan(split);

    // (A1 + A2x)(B1 + B2x) =
    //   A1B1 + ((A1 + A2)(B1 + B2) - A1B1 - A2B2)x + A2B2x^2
    auto low = MulKaratsuba(a_low, b_low);
    auto high = MulKaratsuba(a_high, b_high);
    auto middle =
        Sub(Sub(MulKaratsuba(AddParts(a_low, a_high), AddParts(b_low, b_high)),
                low),
            high);

    std::vector<Elem> result(a.size() + b.size() - 1, Elem::Zero());
    AddShifted(result, low, 0);
    AddShifted(result, middle, split);
    AddShifted(result, high, 2 * split);
    return Trim(std::move(result));
  }

  static void AddShifted(std::vector<Elem>& target,
                         const std::vector<Elem>& value, size_t shift) {
    if (value.empty()) {
      return;
    }
    if (target.size() < value.size() + shift) {
      target.resize(value.size() + shift, Elem::Zero());
    }
    for (size_t i = 0; i < value.size(); ++i) {
      target[i + shift] += value[i];
    }
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

/*! \brief Polynomial adaptation for different multiplication algorithms
 *
 * Follows invariant that is doesn't have leading zeros
 */
template <concepts::GaloisFieldElement Elem,
          concepts::PolynomialEngine<Elem> Engine>
class GenericPolynomial
    : public PolynomialBase<Elem, GenericPolynomial<Elem, Engine>> {
  using Base = PolynomialBase<Elem, GenericPolynomial<Elem, Engine>>;
  friend Base;

 public:
  class Modulus {
    friend GenericPolynomial;

   public:
    Modulus() = default;

   private:
    explicit Modulus(typename Engine::Modulus data)
        : data_(std::move(data)) {
    }

    typename Engine::Modulus data_;
  };

  using Base::Base;
  using Base::Div;
  using Base::Mul;
  using Element = Base::Element;

 public:
  [[nodiscard]]
  GenericPolynomial Mul(const GenericPolynomial& rhs) const& {
    if (data_.empty() || rhs.data_.empty()) {
      return GenericPolynomial();
    }
    return GenericPolynomial(*this).MulInPlace(rhs);
  }

  [[nodiscard]]
  GenericPolynomial Mul(const GenericPolynomial& rhs) && {
    if (data_.empty() || rhs.data_.empty()) {
      data_.clear();
      return std::move(*this);
    }
    MulInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  GenericPolynomial Div(const GenericPolynomial& rhs) const& {
    if (data_.empty()) {
      return GenericPolynomial();
    }
    return GenericPolynomial(*this).DivInPlace(rhs);
  }

  [[nodiscard]]
  GenericPolynomial Div(const GenericPolynomial& rhs) && {
    if (data_.empty()) {
      return std::move(*this);
    }
    DivInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  GenericPolynomial Div(const Modulus& modulus) const& {
    if (data_.empty()) {
      return GenericPolynomial();
    }
    return GenericPolynomial(*this).DivInPlace(modulus);
  }

  [[nodiscard]]
  GenericPolynomial Div(const Modulus& modulus) && {
    if (data_.empty()) {
      return std::move(*this);
    }
    DivInPlace(modulus);
    return std::move(*this);
  }

  [[nodiscard]]
  GenericPolynomial Rem(const GenericPolynomial& rhs) const& {
    if (data_.empty()) {
      return GenericPolynomial();
    }
    return GenericPolynomial(*this).RemInPlace(rhs);
  }

  [[nodiscard]]
  GenericPolynomial Rem(const GenericPolynomial& rhs) && {
    if (data_.empty()) {
      return std::move(*this);
    }
    RemInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  GenericPolynomial Rem(const Modulus& modulus) const& {
    if (data_.empty()) {
      return GenericPolynomial();
    }
    return GenericPolynomial(*this).RemInPlace(modulus);
  }

  [[nodiscard]]
  GenericPolynomial Rem(const Modulus& modulus) && {
    if (data_.empty()) {
      return std::move(*this);
    }
    RemInPlace(modulus);
    return std::move(*this);
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRem(
      const GenericPolynomial& rhs) const& {
    if (data_.empty()) {
      return {GenericPolynomial(), GenericPolynomial()};
    }
    return GenericPolynomial(*this).DivRemInPlace(rhs);
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRem(
      const GenericPolynomial& rhs) && {
    if (data_.empty()) {
      return {GenericPolynomial(), std::move(*this)};
    }
    return std::move(*this).DivRemInPlace(rhs);
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRem(
      const Modulus& modulus) const& {
    if (data_.empty()) {
      return {GenericPolynomial(), GenericPolynomial()};
    }
    return GenericPolynomial(*this).DivRemInPlace(modulus);
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRem(
      const Modulus& modulus) && {
    if (data_.empty()) {
      return {GenericPolynomial(), std::move(*this)};
    }
    return std::move(*this).DivRemInPlace(modulus);
  }

  [[nodiscard]]
  Modulus BuildModulus(size_t max_dividend_size) const {
    return Modulus(Engine::BuildModulus(data_, max_dividend_size));
  }

  [[nodiscard]]
  Modulus BuildModulus() const {
    if (data_.empty()) {
      return BuildModulus(0);
    }
    return BuildModulus(2 * data_.size() - 1);
  }

  [[nodiscard]]
  GenericPolynomial Gcd(GenericPolynomial b) const& {
    return GenericPolynomial(*this).Gcd(std::move(b)).MakeMonic();
  }

  [[nodiscard]]
  GenericPolynomial Gcd(GenericPolynomial b) && {
    auto gcd = Engine::Gcd(std::move(data_), std::move(b.data_));
    return GenericPolynomial(std::move(gcd)).MakeMonic();
  }

 protected:
  using Base::data_;
  using Base::DivInPlace;
  using Base::MulInPlace;

  GenericPolynomial& MulInPlace(const GenericPolynomial& rhs) {
    if (rhs.data_.size() == 1) {
      return Base::MulInPlace(rhs.data_[0]);
    }
    data_ = Engine::Mul(std::move(data_), rhs.data_);
    return *this;
  }

  GenericPolynomial& DivInPlace(const GenericPolynomial& rhs) {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      data_.clear();
      return *this;
    }
    if (m == 1) {
      return Base::DivInPlace(rhs.data_[0]);
    }
    data_ = Engine::Div(std::move(data_), rhs.data_);
    return *this;
  }

  GenericPolynomial& DivInPlace(const Modulus& modulus) {
    const size_t n = data_.size();
    const size_t m = modulus.data_.polynomial.size();

    if (n < m) {
      data_.clear();
      return *this;
    }
    if (m == 1) {
      return Base::DivInPlace(modulus.data_.polynomial[0]);
    }
    data_ = Engine::Div(std::move(data_), modulus.data_);
    return *this;
  }

  GenericPolynomial& RemInPlace(const GenericPolynomial& rhs) {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      return *this;
    }
    if (m == 1) {
      data_.clear();
      return *this;
    }
    data_ = Engine::Rem(std::move(data_), rhs.data_);
    return *this;
  }

  GenericPolynomial& RemInPlace(const Modulus& modulus) {
    const size_t n = data_.size();
    const size_t m = modulus.data_.polynomial.size();

    if (n < m) {
      return *this;
    }
    if (m == 1) {
      data_.clear();
      return *this;
    }
    data_ = Engine::Rem(std::move(data_), modulus.data_);
    return *this;
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRemInPlace(
      const GenericPolynomial& rhs) && {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      return {GenericPolynomial(), std::move(*this)};
    }
    if (m == 1) {
      Base::DivInPlace(rhs.data_[0]);
      return {std::move(*this), GenericPolynomial()};
    }

    auto [quotient, remainder] = Engine::DivRem(std::move(data_), rhs.data_);
    return {GenericPolynomial(std::move(quotient)),
            GenericPolynomial(std::move(remainder))};
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRemInPlace(
      const Modulus& modulus) && {
    const size_t n = data_.size();
    const size_t m = modulus.data_.polynomial.size();

    if (n < m) {
      return {GenericPolynomial(), std::move(*this)};
    }
    if (m == 1) {
      Base::DivInPlace(modulus.data_.polynomial[0]);
      return {std::move(*this), GenericPolynomial()};
    }

    auto [quotient, remainder] =
        Engine::DivRem(std::move(data_), modulus.data_);
    return {GenericPolynomial(std::move(quotient)),
            GenericPolynomial(std::move(remainder))};
  }
};

}  // namespace factorization::polynomial
