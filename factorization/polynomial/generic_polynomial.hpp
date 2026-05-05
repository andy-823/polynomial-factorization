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

#include <utility>
#include <vector>

#include <factorization/concepts.hpp>

#include "polynomial_base.hpp"

namespace factorization::polynomial {

template <concepts::GaloisFieldElement Elem>
struct PolynomialEngine {
  [[nodiscard]]
  static std::vector<Elem> Mul(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
  }

  [[nodiscard]]
  static std::vector<Elem> Rem(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
  }

  [[nodiscard]]
  static std::vector<Elem> Div(std::vector<Elem> a,
                               const std::vector<Elem>& b) {
  }

  [[nodiscard]]
  static std::pair<std::vector<Elem>, std::vector<Elem>> DivRem(
      std::vector<Elem> a, const std::vector<Elem>& b) {
  }

  [[nodiscard]]
  static std::vector<Elem> Gcd(std::vector<Elem> a, std::vector<Elem> b) {
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
  GenericPolynomial Gcd(GenericPolynomial b) const& {
    return GenericPolynomial(*this).Gcd(std::move(b));
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
  using Base::RemoveLeadingZeros;

  GenericPolynomial& MulInPlace(const GenericPolynomial& rhs) {
    if (rhs.data_.size() == 1) {
      return Base::MulInPlace(rhs.data_[0]);
    }
    data_ = Engine::Mul(std::move(data_), rhs.data_);
    return RemoveLeadingZeros();
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
    return RemoveLeadingZeros();
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
    return RemoveLeadingZeros();
  }

  [[nodiscard]]
  std::pair<GenericPolynomial, GenericPolynomial> DivRemInPlace(
      const GenericPolynomial& rhs) && {
    const size_t n = data_.size();
    const size_t m = rhs.data_.size();

    if (n < m) {
      RemoveLeadingZeros();
      return {GenericPolynomial(), std::move(*this)};
    }
    if (m == 1) {
      Base::DivInPlace(rhs.data_[0]);
      RemoveLeadingZeros();
      return {std::move(*this), GenericPolynomial()};
    }

    auto [quotient, remainder] = Engine::DivRem(std::move(data_), rhs.data_);
    return {GenericPolynomial(std::move(quotient)),
            GenericPolynomial(std::move(remainder))};
  }
};

}  // namespace factorization::polynomial
