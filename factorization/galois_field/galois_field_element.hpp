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

#include <array>
#include <concepts>

#include <factorization/concepts.hpp>

namespace factorization::galois_field {

template <std::integral Int>
bool IsZero(Int value) {
  return value == 0;
}

template <typename T>
concept HasIsZero = requires(T value) {
  { value.IsZero() } -> std::same_as<bool>;
};

template <HasIsZero Value>
bool IsZero(const Value& value) {
  return value.IsZero();
}

template <typename Value, std::size_t kSize>
bool IsZero(const std::array<Value, kSize>& value) {
  bool is_zero = true;
  for (const auto& v : value) {
    is_zero = is_zero && IsZero(value);
  }
  return is_zero;
}


/*! \brief Wrapper of GaloisField API
 *
 * all methods are constexpt
 */
template <concepts::GaloisField Field>
class GaloisFieldElement {
 public:
  using Value = typename Field::Value;
  using GaloisField = Field;

 public:
  constexpr GaloisFieldElement(Value value = Value{}) : value_(value) {
  }

  constexpr GaloisFieldElement(const GaloisFieldElement&) = default;
  constexpr GaloisFieldElement(GaloisFieldElement&&) = default;

  constexpr GaloisFieldElement& operator=(const GaloisFieldElement&) = default;
  constexpr GaloisFieldElement& operator=(GaloisFieldElement&&) = default;

  ~GaloisFieldElement() = default;

  constexpr bool operator==(const GaloisFieldElement&) const = default;

  constexpr Value IsZero() const {
    return IsZero(value_);
  }

  constexpr Value Get() const {
    return value_;
  }

  constexpr GaloisFieldElement& operator+=(const GaloisFieldElement& other) {
    value_ = kField.Add(value_, other.value_);
    return *this;
  }

  constexpr GaloisFieldElement& operator-=(const GaloisFieldElement& other) {
    value_ = kField.Add(value_, kField.Negative(other.value_));
    return *this;
  }

  constexpr GaloisFieldElement operator-() const {
    return GaloisFieldElement(kField.Negative(value_));
  }

  constexpr GaloisFieldElement& operator*=(const GaloisFieldElement& other) {
    value_ = kField.Multiply(value_, other.value_);
    return *this;
  }

  constexpr GaloisFieldElement& operator/=(const GaloisFieldElement& other) {
    value_ = kField.Multiply(value_, kField.Inverse(other.value_));
    return *this;
  }

 private:
  constexpr static Field kField;
  Value value_;
};

template <concepts::GaloisField Field>
constexpr typename Field::Value operator+(GaloisFieldElement<Field> first,
                                          GaloisFieldElement<Field> second) {
  return first += second;
}

template <concepts::GaloisField Field>
constexpr typename Field::Value operator-(GaloisFieldElement<Field> first,
                                          GaloisFieldElement<Field> second) {
  return first -= second;
}

template <concepts::GaloisField Field>
constexpr typename Field::Value operator*(GaloisFieldElement<Field> first,
                                          GaloisFieldElement<Field> second) {
  return first *= second;
}

template <concepts::GaloisField Field>
constexpr typename Field::Value operator/(GaloisFieldElement<Field> first,
                                          GaloisFieldElement<Field> second) {
  return first /= second;
}

}  // namespace factorization::galois_field