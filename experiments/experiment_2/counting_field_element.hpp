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

#include <cstdint>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

namespace factorization::galois_field {

/*! \brief Wrapper of GaloisField API
 *
 * Implemented only for experiment
 * Do not use it elsewhere please
 */
template <concepts::GaloisField Field>
class CountingFieldElement {
  inline thread_local static uint64_t actions{0};

  static void Action() {
    ++actions;
  }

 public:
  using Value = typename Field::Value;
  constexpr static bool kCounting = true;

 public:
  constexpr CountingFieldElement() = default;

  template <typename T>
  constexpr CountingFieldElement(T value)
      : value_(std::forward<T>(value)) {
    static_assert(std::constructible_from<Value, T>);
  }

  constexpr CountingFieldElement(const CountingFieldElement&) = default;
  constexpr CountingFieldElement(CountingFieldElement&&) = default;

  constexpr CountingFieldElement& operator=(const CountingFieldElement&) = default;
  constexpr CountingFieldElement& operator=(CountingFieldElement&&) = default;

 ~CountingFieldElement() = default;

  constexpr auto operator<=>(const CountingFieldElement&) const = default;

  constexpr static CountingFieldElement Zero() {
    return CountingFieldElement(kField.Zero());
  }

  constexpr static CountingFieldElement One() {
    return CountingFieldElement(kField.One());
  }

  constexpr static CountingFieldElement AsPolyConstant(Value value) {
    return kField.FieldValueFromConstant(value);
  }

  constexpr Value Get() const {
    return value_;
  }

  static void ResetActions() {
    actions = 0;
  }

  static uint64_t GetActions() {
    return actions;
  }

  constexpr CountingFieldElement& operator+=(
      const CountingFieldElement& other) {
    Action();
    value_ = kField.Add(value_, other.value_);
    return *this;
  }

  constexpr CountingFieldElement& operator-=(
      const CountingFieldElement& other) {
    Action();
    value_ = kField.Sub(value_, other.value_);
    return *this;
  }

  constexpr CountingFieldElement operator-() const {
    Action();
    return CountingFieldElement(kField.Negative(value_));
  }

  constexpr CountingFieldElement& operator*=(
      const CountingFieldElement& other) {
    Action();
    value_ = kField.Multiply(value_, other.value_);
    return *this;
  }

  constexpr CountingFieldElement& operator/=(
      const CountingFieldElement& other) {
    Action();
    value_ = kField.Divide(value_, other.value_);
    return *this;
  }

  constexpr CountingFieldElement Inverse() const {
    Action();
    return CountingFieldElement(kField.Inverse(value_));
  }

  template <typename Power>
  constexpr CountingFieldElement Pow(Power power) const {
    Action();
    return CountingFieldElement(kField.Pow(value_, power));
  }

  constexpr static uint32_t FieldBase() {
    return Field::FieldBase();
  }

  constexpr static uint32_t FieldPower() {
    return Field::FieldPower();
  }

  // std vector is temporary option
  constexpr static std::vector<CountingFieldElement> AllFieldElements() {
    Value current = kField.FirstFieldValue();
    std::vector<CountingFieldElement> result;
    result.emplace_back(current);
    while (current != kField.LastFieldValue()) {
      current = kField.NextFieldValue(current);
      result.emplace_back(current);
    }
    return result;
  }

 private:
  constexpr static Field kField{};
  Value value_;
};

template <concepts::GaloisField Field>
constexpr inline CountingFieldElement<Field> operator+(
    CountingFieldElement<Field> first, CountingFieldElement<Field> second) {
  return first += second;
}

template <concepts::GaloisField Field>
constexpr inline CountingFieldElement<Field> operator-(
    CountingFieldElement<Field> first, CountingFieldElement<Field> second) {
  return first -= second;
}

template <concepts::GaloisField Field>
constexpr inline CountingFieldElement<Field> operator*(
    CountingFieldElement<Field> first, CountingFieldElement<Field> second) {
  return first *= second;
}

template <concepts::GaloisField Field>
constexpr inline CountingFieldElement<Field> operator/(
    CountingFieldElement<Field> first, CountingFieldElement<Field> second) {
  return first /= second;
}

}  // namespace factorization::galois_field