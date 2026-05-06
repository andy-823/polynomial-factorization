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
  static std::vector<Elem> Mul(std::vector<Elem> a,
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
    auto quotient = Div(a, b);
    return Sub(std::move(a), Mul(std::move(quotient), b), b.size() - 1);
  }

  [[nodiscard]]
  static std::vector<Elem> Rem(std::vector<Elem> a, const Modulus& modulus) {
    if (a.size() < modulus.polynomial.size()) {
      return Trim(std::move(a));
    }
    auto quotient = Div(a, modulus);
    return Sub(std::move(a), Mul(std::move(quotient), modulus.polynomial),
               modulus.polynomial.size() - 1);
  }

  // assume a.size() >= b.size()
  [[nodiscard]]
  static std::vector<Elem> Div(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
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
  static std::vector<Elem> Div(std::vector<Elem> a, const Modulus& modulus) {
    if (a.size() < modulus.polynomial.size()) {
      return {};
    }
    const size_t quotient_size = a.size() - modulus.polynomial.size() + 1;
    if (quotient_size > modulus.max_quotient_size) {
      return Div(std::move(a), modulus.polynomial);
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
    auto quotient = Div(a, modulus);
    auto remainder = Sub(std::move(a), Mul(quotient, modulus.polynomial),
                         modulus.polynomial.size() - 1);
    return {std::move(quotient), std::move(remainder)};
  }

  [[nodiscard]]
  static Modulus BuildModulus(const std::vector<Elem>& polynomial,
                              size_t max_dividend_size) {
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
    return EuclidGcd(std::move(a), std::move(b));
  }

 private:
  constexpr static size_t kKaratsubaThreshold = 128;

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
  static std::vector<Elem> MulNaive(const std::vector<Elem>& a,
                                    const std::vector<Elem>& b) {
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
  static std::vector<Elem> MulKaratsuba(const std::vector<Elem>& a,
                                        const std::vector<Elem>& b) {
    if (a.empty() || b.empty()) {
      return {};
    }
    if (b.size() == 1) {
      if (b[0] == Elem::Zero()) {
        return {};
      }
      std::vector<Elem> result(a);
      if (b[0] != Elem::One()) {
        for (auto& value : result) {
          value *= b[0];
        }
      }
      return Trim(std::move(result));
    }
    if (std::min(a.size(), b.size()) <= kKaratsubaThreshold) {
      return MulNaive(a, b);
    }

    const size_t split = std::min(a.size(), b.size()) / 2;
    std::vector<Elem> a_low(a.begin(), a.begin() + split);
    std::vector<Elem> a_high(a.begin() + split, a.end());
    std::vector<Elem> b_low(b.begin(), b.begin() + split);
    std::vector<Elem> b_high(b.begin() + split, b.end());

    // (A1 + A2x)(B1 + B2x) =
    //   A1B1 + ((A1 + A2)(B1 + B2) - A1B1 - A2B2)x + A2B2x^2
    auto low = MulKaratsuba(a_low, b_low);
    auto high = MulKaratsuba(a_high, b_high);
    auto middle = Sub(Sub(MulKaratsuba(Add(std::move(a_low), a_high),
                                       Add(std::move(b_low), b_high)),
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
  static std::vector<Elem> EuclidGcd(std::vector<Elem> a, std::vector<Elem> b) {
    if (a.size() < b.size()) {
      a.swap(b);
    }
    while (!b.empty()) {
      a = Rem(std::move(a), b);
      a.swap(b);
    }
    return a;
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
    return GenericPolynomial(*this).Gcd(std::move(b));
  }

  [[nodiscard]]
  GenericPolynomial Gcd(GenericPolynomial b) && {
    auto gcd = Engine::Gcd(std::move(data_), std::move(b.data_));
    return GenericPolynomial(std::move(gcd));
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
