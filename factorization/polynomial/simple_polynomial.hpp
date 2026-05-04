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
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
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


/*! \brief Polynomial implementation which uses naive algorithms for operations
 *
 * Follows invariant that is doesn't have leading zeros
 */
template <concepts::GaloisFieldElement Elem>
class SimplePolynomial {
 public:
  using Element = Elem;

 public:
  SimplePolynomial() = default;

  // TODO: replace this constructor with better one for convenience
  template <typename T>
  SimplePolynomial(const std::vector<T>& elements) {
    static_assert(std::constructible_from<Element, T>);
    data_.reserve(elements.size());
    for (const auto& element : elements) {
      data_.emplace_back(element);
    }
    RemoveLeadingZeros();
  }

  SimplePolynomial(const std::vector<Element>& elements) 
      : data_(elements) {
    RemoveLeadingZeros();
  }

  SimplePolynomial(std::vector<Element>&& elements) 
      : data_(std::move(elements)) {
    RemoveLeadingZeros();
  }

  SimplePolynomial(const Element& element) : data_(1, element) {
    RemoveLeadingZeros();  // element may be equal to zero
  }

  SimplePolynomial(const SimplePolynomial&) = default;
  SimplePolynomial(SimplePolynomial&&) = default;

  SimplePolynomial& operator=(const SimplePolynomial&) = default;
  SimplePolynomial& operator=(SimplePolynomial&&) = default;

  ~SimplePolynomial() = default;

  auto operator<=>(const SimplePolynomial&) const = default;

  [[nodiscard]]
  SimplePolynomial Neg() const & {
    return SimplePolynomial(*this).Neg();
  }

  [[nodiscard]]
  SimplePolynomial Neg() && {
    for (auto& value : data_) {
      value = -value;
    }
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Add(const SimplePolynomial& rhs) const & {
    if (data_.size() < rhs.data_.size()) {
      return SimplePolynomial(rhs).AddInPlace(*this).RemoveLeadingZeros();
    }
    return SimplePolynomial(*this).AddInPlace(rhs).RemoveLeadingZeros();
  }

  [[nodiscard]]
  SimplePolynomial Add(const SimplePolynomial& rhs) && {
    if (data_.size() < rhs.data_.size()) {
      data_.resize(rhs.data_.size(), Element::Zero());
    }
    AddInPlace(rhs);
    RemoveLeadingZeros();
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Add(const Element& element) const & {
    if (data_.empty() && element != Element::Zero()) {
      return SimplePolynomial({element});
    }
    return SimplePolynomial(*this).Add(element);
  }

  [[nodiscard]]
  SimplePolynomial Add(const Element& element) && {
    if (data_.empty()) {
      data_.emplace_back(element);
    } else {
      data_[0] += element;
    }
    RemoveLeadingZeros();
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Sub(const SimplePolynomial& rhs) const & {
    if (data_.size() < rhs.data_.size()) {
      SimplePolynomial result(*this);
      result.data_.resize(rhs.data_.size(), Element::Zero());
      return result.SubInPlace(rhs).RemoveLeadingZeros();
    }
    return SimplePolynomial(*this).SubInPlace(rhs).RemoveLeadingZeros();
  }

  [[nodiscard]]
  SimplePolynomial Sub(const SimplePolynomial& rhs) && {
    if (data_.size() < rhs.data_.size()) {
      data_.resize(rhs.data_.size(), Element::Zero());
    }
    SubInPlace(rhs);
    RemoveLeadingZeros();
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Sub(const Element& element) const & {
    if (data_.empty() && element != Element::Zero()) {
      return SimplePolynomial({-element});
    }
    return SimplePolynomial(*this).Sub(element);
  }

  [[nodiscard]]
  SimplePolynomial Sub(const Element& element) && {
    if (data_.empty()) {
      data_.emplace_back(-element);
    } else {
      data_[0] -= element;
    }
    RemoveLeadingZeros();
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Mul(const SimplePolynomial& rhs) const & {
    // multiply by zero
    if (data_.empty() || rhs.data_.empty()) {
      return SimplePolynomial();
    }
    return SimplePolynomial(*this).MulInPlace(rhs);
  }

  [[nodiscard]]
  SimplePolynomial Mul(const SimplePolynomial& rhs) && {
    // multiply by zero
    if (data_.empty() || rhs.data_.empty()) {
      data_.clear();
      return std::move(*this);
    }
    MulInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Mul(const Element& element) const & {
    if (data_.empty() || element == Element::Zero()) {
      return SimplePolynomial();
    }
    return SimplePolynomial(*this).MulInPlace(element);
  }

  [[nodiscard]]
  SimplePolynomial Mul(const Element& element) && {
    if (element == Element::Zero()) {
      data_.clear();
      return std::move(*this);
    }
    MulInPlace(element);
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Div(const SimplePolynomial& rhs) const & {
    // we are zero
    if (data_.empty()) {
      return SimplePolynomial();
    }
    return SimplePolynomial(*this).DivInPlace(rhs);
  }

  [[nodiscard]]
  SimplePolynomial Div(const SimplePolynomial& rhs) && {
    // we are zero
    if (data_.empty()) {
      return std::move(*this);
    }
    DivInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Div(const Element& element) const & {
    if (data_.empty()) {
      return SimplePolynomial();
    }
    return SimplePolynomial(*this).DivInPlace(element);
  }

  [[nodiscard]]
  SimplePolynomial Div(const Element& element) && {
    if (data_.empty()) {
      return std::move(*this);
    }
    DivInPlace(element);
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Rem(const SimplePolynomial& rhs) const & {
    // we are zero
    if (data_.empty()) {
      return SimplePolynomial();
    }
    return SimplePolynomial(*this).RemInPlace(rhs);
  }

  [[nodiscard]]
  SimplePolynomial Rem(const SimplePolynomial& rhs) && {
    // we are zero
    if (data_.empty()) {
      return std::move(*this);
    }
    RemInPlace(rhs);
    return std::move(*this);
  }

  [[nodiscard]]
  std::pair<SimplePolynomial, SimplePolynomial>
  DivRem(const SimplePolynomial& rhs) const & {
    if (data_.empty()) {
      return {SimplePolynomial(), SimplePolynomial()};
    }
    return SimplePolynomial(*this).DivRemInPlace(rhs);
  }

  [[nodiscard]]
  std::pair<SimplePolynomial, SimplePolynomial>
  DivRem(const SimplePolynomial& rhs) && {
    if (data_.empty()) {
      return {SimplePolynomial(), std::move(*this)};
    }
    return DivRemInPlace(rhs);
  }

  [[nodiscard]]
  SimplePolynomial Gcd(SimplePolynomial b) const & {
    return SimplePolynomial(*this).Gcd(std::move(b));
  }

  [[nodiscard]]
  SimplePolynomial Gcd(SimplePolynomial b) && {
    while (!b.IsZero()) {
      RemInPlace(b);
      data_.swap(b.data_);
    }
    return std::move(*this).MakeMonic();
  }

  // Makes polynomial monic
  [[nodiscard]]
  SimplePolynomial MakeMonic() const & {
    if (data_.empty()) {
      return SimplePolynomial();
    }
    // division by leading coefficient
    return SimplePolynomial(*this).DivInPlace(data_.back());
  }

  SimplePolynomial MakeMonic() && {
    if (!data_.empty()) {
      DivInPlace(data_.back());
    }
    return std::move(*this);
  }

  [[nodiscard]]
  SimplePolynomial Derivative() const & {
    if (data_.size() <= 1) {
      return SimplePolynomial();
    }
    std::vector<Element> result(data_.size() - 1);
    for (size_t i = 1; i < data_.size(); ++i) {
      result[i - 1] = Element::AsPolyConstant(i) * data_[i];
    }
    return SimplePolynomial(std::move(result));
  }

  [[nodiscard]]
  SimplePolynomial Derivative() && {
    int n = static_cast<int>(data_.size());
    if (n <= 1) {
      data_.clear();
      return std::move(*this);
    }
    for (int i = 1; i < n; ++i) {
      data_[i - 1] = Element::AsPolyConstant(i) * data_[i];
    }
    data_.pop_back();
    RemoveLeadingZeros();
    return std::move(*this);
  }

  [[nodiscard]]
  std::vector<Element> Get() const & {
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
    return data_.size() == 0;
  }

  [[nodiscard]]
  bool IsOne() const noexcept {
    return data_.size() == 1 && data_[0] == Element::One();
  }


 private:
  SimplePolynomial& RemoveLeadingZeros() {
    while (!data_.empty() && data_.back() == Element::Zero()) {
      data_.pop_back();
    }
    return *this;
  }

  // assume data_.size() >= rhs.data_.size()
  SimplePolynomial& AddInPlace(const SimplePolynomial& rhs) {
    for (size_t i = 0; i < rhs.data_.size(); ++i) {
      data_[i] += rhs.data_[i];
    }
    return *this;
  }

  // assume data_.size() >= rhs.data_.size()
  SimplePolynomial& SubInPlace(const SimplePolynomial& rhs) {
    for (size_t i = 0; i < rhs.data_.size(); ++i) {
      data_[i] -= rhs.data_[i];
    }
    return *this;
  }

  // Both are nonzero
  SimplePolynomial& MulInPlace(const SimplePolynomial& rhs) {
    const int n = static_cast<int>(data_.size());
    const int m = static_cast<int>(rhs.data_.size());

    // multiplication by constant
    if (m == 1) {
      return MulInPlace(rhs.data_[0]);
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

  // assume element is not zero
  SimplePolynomial& MulInPlace(const Element& element) {
    if (element != Element::One()) [[likely]] {
      for (auto& value : data_) {
        value *= element;
      }
    }
    return *this;
  }

  // assume division is not by zero
  SimplePolynomial& DivInPlace(const SimplePolynomial& rhs) {
    const int n = static_cast<int>(data_.size());
    const int m = static_cast<int>(rhs.data_.size());

    if (n < m) {
      data_.clear();
      return *this;
    }
    if (m == 1) {
      return DivInPlace(rhs.data_[0]);
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    const int quotient_size = n - m + 1;
    std::vector<Element> quotient(quotient_size);

    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    // will perform naive polinomial division
    // go from greatest power to lowest
    // we have smth like this every step
    //   a[0] + ... + a[k - 1] + a[k] + ... + a[n]
    // minus
    //                           b[0] + ... + b[n - k]
    for (int i = quotient_size - 1; i >= 0; --i) {
      Element coeff = a[i + m - 1] * inv_lead;
      quotient[i] = coeff;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (int j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_ = std::move(quotient);
    return *this;
  }

  // assume element is not zero
  SimplePolynomial& DivInPlace(const Element& element) {
    if (element != Element::One()) {
      return MulInPlace(element.Inverse());
    }
    return *this;
  }

  // assume division is not by zero
  SimplePolynomial& RemInPlace(const SimplePolynomial& rhs) {
    const int n = static_cast<int>(data_.size());
    const int m = static_cast<int>(rhs.data_.size());

    if (n < m) {
      return *this;
    }
    if (m == 1) {
      data_.clear();
      return *this;
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    const int quotient_size = n - m + 1;
    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    // will perform naive polinomial division
    // go from greatest power to lowest
    // we have smth like this every step
    //   a[0] + ... + a[k - 1] + a[k] + ... + a[n]
    // minus
    //                           b[0] + ... + b[n - k]
    for (int i = quotient_size - 1; i >= 0; --i) {
      Element coeff = a[i + m - 1] * inv_lead;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (int j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_.resize(m - 1);
    RemoveLeadingZeros();
    return *this;
  }

  // assume division is not by zero, <quotient, remainder>
  [[nodiscard]]
  std::pair<SimplePolynomial, SimplePolynomial> 
  DivRemInPlace(const SimplePolynomial& rhs) && {
    const int n = static_cast<int>(data_.size());
    const int m = static_cast<int>(rhs.data_.size());
    if (n < m) {
      RemoveLeadingZeros();
      return {SimplePolynomial(), std::move(*this)};
    }
    if (m == 1) {
      DivInPlace(rhs.data_[0]);
      RemoveLeadingZeros();
      return {std::move(*this), SimplePolynomial()};
    }
    const Element inv_lead = rhs.data_.back().Inverse();
    const int quotient_size = n - m + 1;
    std::vector<Element> quotient(quotient_size);

    Element* a = data_.data();
    const Element* b = rhs.data_.data();
    for (int i = quotient_size - 1; i >= 0; --i) {
      Element coeff = a[i + m - 1] * inv_lead;
      quotient[i] = coeff;
      if (coeff == Element::Zero()) [[unlikely]] {
        continue;
      }
      for (int j = 0; j < m - 1; ++j) {
        a[i + j] -= coeff * b[j];
      }
    }
    data_.resize(m - 1);
    RemoveLeadingZeros();
    SimplePolynomial remainder(std::move(*this));
    return {SimplePolynomial(std::move(quotient)), std::move(remainder)};
  }

  std::vector<Element> data_;
};

}  // namespace factorization::polynomial