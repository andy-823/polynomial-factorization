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

#include <cstdint>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

namespace factorization::galois_field {

/*! \brief Wrapper of GaloisField API
 *
 * Expects all Field methods to be constexpr
 */
template <concepts::GaloisField Field>
class FieldElementWrapper {
 public:
  using Value = typename Field::Value;

 public:
  constexpr inline FieldElementWrapper() = default;

  template <typename T>
  constexpr inline FieldElementWrapper(T value)
      : value_(std::forward<T>(value)) {
    static_assert(std::constructible_from<Value, T>);
  }

  constexpr inline FieldElementWrapper(const FieldElementWrapper&) = default;
  constexpr inline FieldElementWrapper(FieldElementWrapper&&) = default;

  constexpr inline FieldElementWrapper& operator=(const FieldElementWrapper&) = default;
  constexpr inline FieldElementWrapper& operator=(FieldElementWrapper&&) = default;

  inline ~FieldElementWrapper() = default;

  constexpr inline bool operator==(const FieldElementWrapper&) const = default;

  constexpr static inline FieldElementWrapper Zero() {
    return FieldElementWrapper(kField.Zero());
  }

  constexpr static inline FieldElementWrapper One() {
    return FieldElementWrapper(kField.One());
  }

  constexpr static inline FieldElementWrapper AsPolyConstant(Value value) {
    return kField.FieldValueFromConstant(value);
  }

  constexpr inline Value Get() const {
    return value_;
  }

  constexpr inline FieldElementWrapper& operator+=(
      const FieldElementWrapper& other) {
    value_ = kField.Add(value_, other.value_);
    return *this;
  }

  constexpr inline FieldElementWrapper& operator-=(
      const FieldElementWrapper& other) {
    value_ = kField.Sub(value_, other.value_);
    return *this;
  }

  constexpr inline FieldElementWrapper operator-() const {
    return FieldElementWrapper(kField.Negative(value_));
  }

  constexpr inline FieldElementWrapper& operator*=(
      const FieldElementWrapper& other) {
    value_ = kField.Multiply(value_, other.value_);
    return *this;
  }

  constexpr inline FieldElementWrapper& operator/=(
      const FieldElementWrapper& other) {
    value_ = kField.Divide(value_, other.value_);
    return *this;
  }

  constexpr inline FieldElementWrapper Inverse() const {
    return FieldElementWrapper(kField.Inverse(value_));
  }

  template <typename Power>
  constexpr inline FieldElementWrapper Pow(Power power) const {
    return FieldElementWrapper(kField.Pow(value_, power));
  }

  constexpr static inline uint32_t FieldBase() {
    return Field::FieldBase();
  }

  constexpr static inline uint32_t FieldPower() {
    return Field::FieldPower();
  }

  // std vector is temporary option
  constexpr static inline std::vector<FieldElementWrapper> AllFieldElements() {
    Value current = kField.FirstFieldValue();
    std::vector<FieldElementWrapper> result;
    result.emplace_back(current);
    while (current != kField.LastFieldValue()) {
      current = kField.NextFieldValue(current);
      result.emplace_back(current);
    }
    return result;
  }

 private:
  constexpr static inline Field kField{};
  Value value_;
};

template <concepts::GaloisField Field>
constexpr inline FieldElementWrapper<Field> operator+(
    FieldElementWrapper<Field> first, FieldElementWrapper<Field> second) {
  return first += second;
}

template <concepts::GaloisField Field>
constexpr inline FieldElementWrapper<Field> operator-(
    FieldElementWrapper<Field> first, FieldElementWrapper<Field> second) {
  return first -= second;
}

template <concepts::GaloisField Field>
constexpr inline FieldElementWrapper<Field> operator*(
    FieldElementWrapper<Field> first, FieldElementWrapper<Field> second) {
  return first *= second;
}

template <concepts::GaloisField Field>
constexpr inline FieldElementWrapper<Field> operator/(
    FieldElementWrapper<Field> first, FieldElementWrapper<Field> second) {
  return first /= second;
}

}  // namespace factorization::galois_field