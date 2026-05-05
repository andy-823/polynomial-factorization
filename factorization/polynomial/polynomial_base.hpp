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

namespace factorization::polynomial {

/*! \brief Implementation of basic polynomial operations
 * Implementation of this operation is common for all polynomials
 * Follows invariant that is doesn't have leading zeros
 */
template <concepts::GaloisFieldElement Elem, typename Derived>
class PolynomialBase {
 public:
  using Element = Elem;

 public:
  PolynomialBase() = default;

  // TODO: replace this constructor with better one for convenience
  template <typename T>
  explicit PolynomialBase(const std::vector<T>& elements) {
    static_assert(std::constructible_from<Element, T>);
    data_.reserve(elements.size());
    for (const auto& element : elements) {
      data_.emplace_back(element);
    }
    RemoveLeadingZeros();
  }

  explicit PolynomialBase(const std::vector<Element>& elements)
      : data_(elements) {
    RemoveLeadingZeros();
  }

  explicit PolynomialBase(std::vector<Element>&& elements)
      : data_(std::move(elements)) {
    RemoveLeadingZeros();
  }

  explicit PolynomialBase(const Element& element)
      : data_(1, element) {
    RemoveLeadingZeros();  // element may be equal to zero
  }

  PolynomialBase(const PolynomialBase&) = default;
  PolynomialBase(PolynomialBase&&) = default;

  PolynomialBase& operator=(const PolynomialBase&) = default;
  PolynomialBase& operator=(PolynomialBase&&) = default;

  ~PolynomialBase() = default;

  auto operator<=>(const PolynomialBase&) const = default;

  [[nodiscard]]
  Derived Neg() const& {
    return Derived(static_cast<const Derived&>(*this)).Neg();
  }

  [[nodiscard]]
  Derived Neg() && {
    for (auto& value : data_) {
      value = -value;
    }
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Add(const Derived& rhs) const& {
    if (data_.size() < rhs.data_.size()) {
      return Derived(rhs)
          .AddInPlace(static_cast<const Derived&>(*this))
          .RemoveLeadingZeros();
    }
    return Derived(static_cast<const Derived&>(*this))
        .AddInPlace(rhs)
        .RemoveLeadingZeros();
  }

  [[nodiscard]]
  Derived Add(const Derived& rhs) && {
    if (data_.size() < rhs.data_.size()) {
      data_.resize(rhs.data_.size(), Element::Zero());
    }
    AddInPlace(rhs);
    RemoveLeadingZeros();
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Add(const Element& element) const& {
    if (data_.empty() && element != Element::Zero()) {
      return Derived(element);
    }
    return Derived(static_cast<const Derived&>(*this)).Add(element);
  }

  [[nodiscard]]
  Derived Add(const Element& element) && {
    if (data_.empty()) {
      data_.emplace_back(element);
    } else {
      data_[0] += element;
    }
    RemoveLeadingZeros();
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Sub(const Derived& rhs) const& {
    if (data_.size() < rhs.data_.size()) {
      Derived result(static_cast<const Derived&>(*this));
      result.data_.resize(rhs.data_.size(), Element::Zero());
      return result.SubInPlace(rhs).RemoveLeadingZeros();
    }
    return Derived(static_cast<const Derived&>(*this))
        .SubInPlace(rhs)
        .RemoveLeadingZeros();
  }

  [[nodiscard]]
  Derived Sub(const Derived& rhs) && {
    if (data_.size() < rhs.data_.size()) {
      data_.resize(rhs.data_.size(), Element::Zero());
    }
    SubInPlace(rhs);
    RemoveLeadingZeros();
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Sub(const Element& element) const& {
    if (data_.empty() && element != Element::Zero()) {
      return Derived(-element);
    }
    return Derived(static_cast<const Derived&>(*this)).Sub(element);
  }

  [[nodiscard]]
  Derived Sub(const Element& element) && {
    if (data_.empty()) {
      data_.emplace_back(-element);
    } else {
      data_[0] -= element;
    }
    RemoveLeadingZeros();
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Mul(const Element& element) const& {
    if (data_.empty() || element == Element::Zero()) {
      return Derived();
    }
    return Derived(static_cast<const Derived&>(*this)).MulInPlace(element);
  }

  [[nodiscard]]
  Derived Mul(const Element& element) && {
    if (element == Element::Zero()) {
      data_.clear();
      return std::move(static_cast<Derived&>(*this));
    }
    MulInPlace(element);
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Div(const Element& element) const& {
    if (data_.empty()) {
      return Derived();
    }
    return Derived(static_cast<const Derived&>(*this)).DivInPlace(element);
  }

  [[nodiscard]]
  Derived Div(const Element& element) && {
    if (data_.empty()) {
      return std::move(static_cast<Derived&>(*this));
    }
    DivInPlace(element);
    return std::move(static_cast<Derived&>(*this));
  }

  // Makes polynomial monic
  [[nodiscard]]
  Derived MakeMonic() const& {
    if (data_.empty()) {
      return Derived();
    }
    // division by leading coefficient
    return Derived(static_cast<const Derived&>(*this)).DivInPlace(data_.back());
  }

  Derived MakeMonic() && {
    if (!data_.empty()) {
      DivInPlace(data_.back());
    }
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  Derived Derivative() const& {
    if (data_.size() <= 1) {
      return Derived();
    }
    std::vector<Element> result(data_.size() - 1);
    for (size_t i = 1; i < data_.size(); ++i) {
      result[i - 1] = Element(i) * data_[i];
    }
    return Derived(std::move(result));
  }

  [[nodiscard]]
  Derived Derivative() && {
    size_t n = data_.size();
    if (n <= 1) {
      data_.clear();
      return std::move(static_cast<Derived&>(*this));
    }
    for (size_t i = 1; i < n; ++i) {
      data_[i - 1] = Element(i) * data_[i];
    }
    data_.pop_back();
    RemoveLeadingZeros();
    return std::move(static_cast<Derived&>(*this));
  }

  [[nodiscard]]
  std::vector<Element> Get() const& {
    return data_;
  }

  [[nodiscard]]
  std::vector<Element>&& Get() && {
    return std::move(data_);
  }

  [[nodiscard]]
  size_t Size() const noexcept {
    return data_.size();
  }

  [[nodiscard]]
  bool IsZero() const noexcept {
    return data_.empty();
  }

  [[nodiscard]]
  bool IsOne() const noexcept {
    return data_.size() == 1 && data_[0] == Element::One();
  }

 protected:
  Derived& RemoveLeadingZeros() {
    while (!data_.empty() && data_.back() == Element::Zero()) {
      data_.pop_back();
    }
    return static_cast<Derived&>(*this);
  }

  // assume data_.size() >= rhs.data_.size()
  Derived& AddInPlace(const Derived& rhs) {
    for (size_t i = 0; i < rhs.data_.size(); ++i) {
      data_[i] += rhs.data_[i];
    }
    return static_cast<Derived&>(*this);
  }

  // assume data_.size() >= rhs.data_.size()
  Derived& SubInPlace(const Derived& rhs) {
    for (size_t i = 0; i < rhs.data_.size(); ++i) {
      data_[i] -= rhs.data_[i];
    }
    return static_cast<Derived&>(*this);
  }

  // assume element is not zero
  Derived& MulInPlace(const Element& element) {
    if (element != Element::One()) [[likely]] {
      for (auto& value : data_) {
        value *= element;
      }
    }
    return static_cast<Derived&>(*this);
  }

  // assume element is not zero
  Derived& DivInPlace(const Element& element) {
    if (element != Element::One()) {
      return MulInPlace(element.Inverse());
    }
    return static_cast<Derived&>(*this);
  }

  std::vector<Element> data_;
};

}  // namespace factorization::polynomial
