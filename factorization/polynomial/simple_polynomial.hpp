// MIT License
//
// Copyright (c) 2025 Andrei Ishutin
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
  using Value = Element::Value;

 public:
  SimplePolynomial() = default;

  SimplePolynomial(const std::vector<Value>& elements) {
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

  bool operator==(const SimplePolynomial&) const = default;

  bool operator<(const SimplePolynomial& other) const {
    if (data_.size() != other.data_.size()) {
      return data_.size() < other.data_.size();
    }
    for (size_t i = 0; i < data_.size(); ++i) {
      if (data_[i].Get() < other.data_[i].Get()) {
        return true;
      }
    }
    return false;  // equal
  }

  SimplePolynomial& operator+=(const SimplePolynomial& other) {
    if (data_.size() < other.data_.size()) {
      data_.resize(other.data_.size(), Element::Zero());
    }
    for (size_t i = 0; i < other.data_.size(); ++i) {
      data_[i] += other.data_[i];
    }
    RemoveLeadingZeros();
    return *this;
  }

  SimplePolynomial& operator+=(const Element& element) {
    if (data_.empty()) {
      data_.emplace_back(element);
    } else {
      data_[0] += element;
    }
    RemoveLeadingZeros();
    return *this;
  }

  SimplePolynomial& operator-=(const SimplePolynomial& other) {
    if (data_.size() < other.data_.size()) {
      data_.resize(other.data_.size(), Element::Zero());
    }
    for (size_t i = 0; i < other.data_.size(); ++i) {
      data_[i] -= other.data_[i];
    }
    RemoveLeadingZeros();
    return *this;
  }

  SimplePolynomial& operator-=(const Element& element) {
    if (data_.empty()) {
      data_.emplace_back(element);
    }
    else {
      data_[0] -= element;
    }
    RemoveLeadingZeros();
    return *this;
  }

  SimplePolynomial& operator*=(const SimplePolynomial& other) {
    // check if result is zero
    if (data_.empty() || other.data_.empty()) {
      data_.clear();
      return *this;
    }
    // other is actually constant
    if (other.data_.size() == 1) {
      for (auto& value : data_) {
        value *= other.data_[0];
      }
      return *this;
    }
    size_t result_power = data_.size() + other.data_.size() - 1;
    std::vector<Element> result(result_power, Element::Zero());
    for (size_t power = 0; power < result_power; ++power) {
      for (size_t i = 0; i <= power; ++i) {
        result[power] += data_[i] * other.data_[power - i];
      }
    }
    data_ = std::move(result);
    return *this;
  }

  SimplePolynomial& operator*=(const Element& element) {
    // check if result is zero
    if (element == Element::Zero()) {
      data_.clear();
    } else {
      for (auto& value : data_) {
        value *= element;
      }
    }
    return *this;
  }

  SimplePolynomial& operator/=(const SimplePolynomial& other) {
    // result is zero
    if (data_.size() < other.data_.size()) {
      data_.clear();
      return *this;
    }
    // other is actually constant
    if (other.data_.size() == 1) {
      return *this /= other.data_[0];
    }
    size_t result_size = data_.size() - other.data_.size() + 1;
    std::vector<Element> result(result_size);
    for (int power = static_cast<int>(result_size - 1); power >= 0; --power) {
      // will perform naive polinomial division
      // go from greatest power to lowest
      // we have smth like this every step
      //   a[0] + ... + a[k - 1] + a[k] + ... + a[n]
      // minus
      //                           b[0] + ... + b[n - k]
      // probably the best idea how to make it look fine
      // is to use iterators
      auto divident = data_.rbegin() + result_size - 1 - power;
      auto divisor = other.data_.rbegin();

      Element coefficient = *divident / *divisor;
      result[power] = coefficient;

      if (coefficient == Element::Zero()) {
        continue;
      }
      while (divisor != other.data_.rend()) {
        *divident -= *divisor * coefficient;
        --divident;
        --divisor;
      }
    }
    data_ = std::move(result);
    return *this;
  }

  SimplePolynomial& operator/=(const Element& element) {
    // Divide by zero is UB
    auto inverse = element.Inverse();
    for (auto& value : data_) {
      value *= inverse;
    }
    return *this;
  }

  SimplePolynomial& operator%=(const SimplePolynomial& other) {
    // nothing to do
    if (data_.size() < other.data_.size()) {
      return *this;
    }
    // Remainder by zero is UB
    size_t num_steps = data_.size() - other.data_.size() + 1;
    for (size_t step = 0; step < num_steps; ++step) {
      // act almost the same way as division
      auto divident = data_.rbegin() + step;
      auto divisor = other.data_.rbegin();

      Element coefficient = *divident / *divisor;
      if (coefficient == Element::Zero()) {
        continue;
      }
      while (divisor != other.data_.rend()) {
        *divident -= *divisor * coefficient;
        --divident;
        --divisor;
      }
    }
    RemoveLeadingZeros();
    return *this;
  }

  SimplePolynomial operator-() const {
    std::vector<Element> result = data_;
    for (auto& value : result) {
      value = -value;
    }
    return SimplePolynomial(std::move(result));
  }

  std::vector<Element> GetElements() const {
    return data_;
  }

  size_t Size() const {
    return data_.size();
  }

  SimplePolynomial Derivative() const {
    if (data_.size() <= 1) {
      return SimplePolynomial();
    }
    std::vector<Element> result(data_.size() - 1);
    for (size_t i = 1; i < data_.size(); ++i) {
      result[i - 1] = Element::AsPolyConstant(i) * data_[i];
    }
    return SimplePolynomial(std::move(result));
  }

  // Makes polynomial monic
  void MakeMonic() {
    if (data_.empty()) {
      return;
    }
    auto leading = data_.back();
    if (leading == Element::One()) {  // already monic
      return;
    }
    for (auto& element : data_) {
      element /= leading;
    }
  }

  bool IsZero() const {
    return data_.size() == 0;
  }

  bool IsOne() const {
    return data_.size() == 1 && data_[0] == Element::One();
  }

 private:
  void RemoveLeadingZeros() {
    size_t new_size = data_.size();
    while (new_size > 0 && data_[new_size - 1] == Element::Zero()) {
      --new_size;
    }
    data_.resize(new_size);
  }

  std::vector<Element> data_;
};

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator+(SimplePolynomial<Element> first,
                                    const SimplePolynomial<Element>& second) {
  return first += second;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator+(SimplePolynomial<Element> poly,
                                    Element value) {
  return poly += value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator+(Element value,
                                    SimplePolynomial<Element> poly) {
  return poly += value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator-(SimplePolynomial<Element> first,
                                    const SimplePolynomial<Element>& second) {
  return first -= second;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator-(SimplePolynomial<Element> poly,
                                    Element value) {
  return poly -= value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator-(Element value,
                                    const SimplePolynomial<Element>& poly) {
  return -poly += value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator*(SimplePolynomial<Element> first,
                                    const SimplePolynomial<Element>& second) {
  return first *= second;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator*(SimplePolynomial<Element> poly,
                                    Element value) {
  return poly *= value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator*(Element value,
                                    SimplePolynomial<Element> poly) {
  return poly *= value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator/(SimplePolynomial<Element> first,
                                    const SimplePolynomial<Element>& second) {
  return first /= second;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator/(SimplePolynomial<Element> poly,
                                    Element value) {
  return poly /= value;
}

template <concepts::GaloisFieldElement Element>
SimplePolynomial<Element> operator%(SimplePolynomial<Element> first,
                                    const SimplePolynomial<Element>& second) {
  return first %= second;
}

}  // namespace factorization::polynomial