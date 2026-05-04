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
  constexpr static uint32_t kFieldBase = Field::FieldBase();
  constexpr static uint32_t kFieldPower = Field::FieldPower();

 public:
  using Coefficient = typename Field::Coefficient;

 public:
  constexpr CountingFieldElement() = default;

  template <std::convertible_to<Coefficient> T = Coefficient>
  constexpr CountingFieldElement(T value)
      : value_(kField.Encode(value)) {
  }

  template <std::convertible_to<Coefficient> T = Coefficient>
  constexpr CountingFieldElement(const std::array<T, kFieldPower>& value) 
      : value_(kField.Encode(Convert(value))) {
  }

  constexpr std::array<Coefficient, kFieldPower> Get() const {
    return kField.Decode(value_);
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
    return Construct(kField.Negative(value_));
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
    return Construct(kField.Inverse(value_));
  }

  template <typename Power>
  constexpr CountingFieldElement Pow(Power power) const {
    Action();
    return Construct(kField.Pow(value_, power));
  }

  constexpr static uint32_t FieldBase() {
    return Field::FieldBase();
  }

  constexpr static uint32_t FieldPower() {
    return Field::FieldPower();
  }

  // std vector is temporary option
  constexpr static std::vector<CountingFieldElement> AllFieldElements() {
    constexpr uint64_t kFieldSize = utils::BinPow(kFieldBase, kFieldPower);

    std::vector<CountingFieldElement> result;
    std::array<Coefficient, kFieldPower> coeffs{};
    for (uint64_t i = 0; i < kFieldSize; ++i) {
      uint64_t val = i;
      for (size_t i = 0; i < kFieldPower; ++i) {
        coeffs[i] = static_cast<Coefficient>(val % kFieldBase);
        val /= kFieldBase;
      }
      result.emplace_back(CountingFieldElement(coeffs));
    }
    return result;
  }

 private:
  template <std::convertible_to<Coefficient> T>
  constexpr std::array<Coefficient, kFieldPower>
  Convert(const std::array<T, kFieldPower>& value) {
    std::array<Coefficient, kFieldPower> result;
    for (size_t i = 0; i < kFieldPower; ++i) {
      result[i] = value[i];
    }
    return result;
  }

  static CountingFieldElement Construct(Value value) {
    CountingFieldElement result;
    result.value_ = value;
    return result;
  }

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