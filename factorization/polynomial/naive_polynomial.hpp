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

#include <vector>

#include <factorization/concepts.hpp>

#include "polynomial_base.hpp"

namespace factorization::polynomial {

/*! \brief Polynomial implementation which uses naive algorithms for operations
 *
 * Follows invariant that is doesn't have leading zeros
 */
template <concepts::GaloisFieldElement Elem>
class NaivePolynomial : public PolynomialBase<Elem, NaivePolynomial<Elem>> {
  using Base = PolynomialBase<Elem, NaivePolynomial<Elem>>;
  friend Base;

 public:
  using Base::Base;
  using Element = Base::Element;
  using Base::Div;
  using Base::Mul;
  using Modulus = NaivePolynomial;

 public:
  [[nodiscard]]
  NaivePolynomial Mul(const NaivePolynomial& rhs) const& {
    // multiply by zero
    if (data_.empty() || rhs.data_.empty()) {
      return NaivePolynomial();
    }
    return NaivePolynomial(*this).MulInPlace(rhs);
  }

  [[nodiscard]]
  NaivePolynomial Mul(const NaivePolynomial& rhs) && {
    // multiply by zero
    if (data_.empty() || rhs.data_.empty()) {
      data_.clear();
      return std::move(*this);
    }
    MulInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  NaivePolynomial Div(const NaivePolynomial& rhs) const& {
    // we are zero
    if (data_.empty()) {
      return NaivePolynomial();
    }
    return NaivePolynomial(*this).DivInPlace(rhs);
  }

  [[nodiscard]]
  NaivePolynomial Div(const NaivePolynomial& rhs) && {
    // we are zero
    if (data_.empty()) {
      return std::move(*this);
    }
    DivInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  NaivePolynomial Rem(const NaivePolynomial& rhs) const& {
    // we are zero
    if (data_.empty()) {
      return NaivePolynomial();
    }
    return NaivePolynomial(*this).RemInPlace(rhs);
  }

  [[nodiscard]]
  NaivePolynomial Rem(const NaivePolynomial& rhs) && {
    // we are zero
    if (data_.empty()) {
      return std::move(*this);
    }
    RemInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  std::pair<NaivePolynomial, NaivePolynomial> DivRem(
      const NaivePolynomial& rhs) const& {
    if (data_.empty()) {
      return {NaivePolynomial(), NaivePolynomial()};
    }
    return NaivePolynomial(*this).DivRemInPlace(rhs);
  }

  [[nodiscard]]
  std::pair<NaivePolynomial, NaivePolynomial> DivRem(
      const NaivePolynomial& rhs) && {
    if (data_.empty()) {
      return {NaivePolynomial(), std::move(*this)};
    }
    return DivRemInPlace(rhs);
  }

  [[nodiscard]]
  NaivePolynomial Gcd(NaivePolynomial b) const& {
    return NaivePolynomial(*this).Gcd(std::move(b));
  }

  [[nodiscard]]
  NaivePolynomial Gcd(NaivePolynomial b) && {
    while (!b.IsZero()) {
      RemInPlace(b);
      data_.swap(b.data_);
    }
    return std::move(*this).MakeMonic();
  }

  Modulus BuildModulus(size_t) const& {
    return *this;
  }

  Modulus BuildModulus(size_t) && {
    return std::move(*this);
  }

 protected:
  using Base::data_;
  using Base::DivInPlace;
  using Base::MulInPlace;
  using Base::RemoveLeadingZeros;

  // Both are nonzero
  NaivePolynomial& MulInPlace(const NaivePolynomial& rhs) {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    // multiplication by constant
    if (m == 1) {
      return Base::MulInPlace(rhs.data_[0]);
    }

    std::vector<Element> result(n + m - 1, Element::Zero());

    const auto* a = data_.data();
    const auto* b = rhs.data_.data();
    auto* res = result.data();

    for (size_t i = 0; i < n; ++i) {
      Element coeff = a[i];
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < m; ++j) {
        res[i + j] += coeff * b[j];
      }
    }
    data_ = std::move(result);
    return *this;
  }

  // assume division is not by zero
  NaivePolynomial& DivInPlace(const NaivePolynomial& rhs) {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      data_.clear();
      return *this;
    }
    if (m == 1) {
      return Base::DivInPlace(rhs.data_[0]);
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    const size_t quotient_size = n - m + 1;
    std::vector<Element> quotient(quotient_size);

    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    // will perform naive polinomial division
    // go from greatest power to lowest
    // we have smth like this every step
    //   a[0] + ... + a[k - 1] + a[k] + ... + a[n]
    // minus
    //                           b[0] + ... + b[n - k]
    for (size_t i = quotient_size; i-- > 0;) {
      Element coeff = a[i + m - 1] * inv_lead;
      quotient[i] = coeff;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_ = std::move(quotient);
    return *this;
  }

  // assume division is not by zero
  NaivePolynomial& RemInPlace(const NaivePolynomial& rhs) {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      return *this;
    }
    if (m == 1) {
      data_.clear();
      return *this;
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    // n >= m -> quotient_size >= 1
    const size_t quotient_size = n - m + 1;
    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    // will perform naive polinomial division
    // go from greatest power to lowest
    // we have smth like this every step
    //   a[0] + ... + a[k - 1] + a[k] + ... + a[n]
    // minus
    //                           b[0] + ... + b[n - k]
    for (size_t i = quotient_size; i-- > 0;) {
      Element coeff = a[i + m - 1] * inv_lead;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_.resize(m - 1);
    RemoveLeadingZeros();
    return *this;
  }

  // assume division is not by zero, <quotient, remainder>
  [[nodiscard]]
  std::pair<NaivePolynomial, NaivePolynomial> DivRemInPlace(
      const NaivePolynomial& rhs) && {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();
    if (n < m) {
      RemoveLeadingZeros();
      return {NaivePolynomial(), std::move(*this)};
    }
    if (m == 1) {
      Base::DivInPlace(rhs.data_[0]);
      RemoveLeadingZeros();
      return {std::move(*this), NaivePolynomial()};
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    const size_t quotient_size = n - m + 1;
    std::vector<Element> quotient(quotient_size);

    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    for (size_t i = quotient_size; i-- > 0;) {
      Element coeff = a[i + m - 1] * inv_lead;
      quotient[i] = coeff;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (size_t j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_.resize(m - 1);
    RemoveLeadingZeros();
    NaivePolynomial remainder(std::move(*this));
    return {NaivePolynomial(std::move(quotient)), std::move(remainder)};
  }
};

}  // namespace factorization::polynomial
