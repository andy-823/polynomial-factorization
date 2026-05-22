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

#include <array>
#include <cstdint>
#include <iterator>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

namespace factorization::galois_field {

/*! \brief Wrapper of GaloisField API
 *
 * Expects all Field methods to be constexpr
 */
template <concepts::GaloisField Field>
class FieldElementWrapper {
  using Value = typename Field::Value;
  constexpr static auto kFieldBase = Field::FieldBase();
  constexpr static size_t kFieldPower = Field::FieldPower();

 public:
  using Coefficient = typename Field::Coefficient;

 private:
  struct ElementsRange {
    constexpr static uint64_t kFieldSize =
        utils::BinPow(kFieldBase, kFieldPower);

    struct Iterator {
      using iterator_concept = std::input_iterator_tag;  // NOLINT
      using value_type = FieldElementWrapper;            // NOLINT
      using difference_type = std::ptrdiff_t;            // NOLINT

      constexpr value_type operator*() const {
        std::array<Coefficient, kFieldPower> coeffs{};
        uint64_t val = value;
        for (size_t i = 0; i < kFieldPower; ++i) {
          coeffs[i] = static_cast<Coefficient>(val % kFieldBase);
          val /= kFieldBase;
        }
        return FieldElementWrapper(coeffs);
      }

      constexpr Iterator& operator++() {
        ++value;
        return *this;
      }

      constexpr void operator++(int) {
        ++(*this);
      }

      constexpr bool operator==(const Iterator&) const = default;

      uint64_t value;
    };

    constexpr Iterator begin() const {
      return {0};
    }  // NOLINT
    constexpr Iterator end() const {
      return {kFieldSize};
    }  // NOLINT
  };

 public:
  constexpr FieldElementWrapper() = default;

  template <std::convertible_to<Coefficient> T = Coefficient>
  explicit constexpr FieldElementWrapper(T value)
      : value_(kField.Encode(value)) {
  }

  template <std::convertible_to<Coefficient> T = Coefficient>
  explicit constexpr FieldElementWrapper(
      const std::array<T, kFieldPower>& value)
      : value_(kField.Encode(Convert(value))) {
  }

  constexpr std::array<Coefficient, kFieldPower> Get() const {
    return kField.Decode(value_);
  }

  constexpr FieldElementWrapper(const FieldElementWrapper&) = default;
  constexpr FieldElementWrapper(FieldElementWrapper&&) = default;

  constexpr FieldElementWrapper& operator=(const FieldElementWrapper&) =
      default;
  constexpr FieldElementWrapper& operator=(FieldElementWrapper&&) = default;

  ~FieldElementWrapper() = default;

  constexpr auto operator<=>(const FieldElementWrapper&) const = default;

  constexpr static FieldElementWrapper Zero() {
    return Construct(kField.Zero());
  }

  constexpr static FieldElementWrapper One() {
    return Construct(kField.One());
  }

  constexpr FieldElementWrapper& operator+=(const FieldElementWrapper& other) {
    value_ = kField.Add(value_, other.value_);
    return *this;
  }

  constexpr FieldElementWrapper& operator-=(const FieldElementWrapper& other) {
    value_ = kField.Sub(value_, other.value_);
    return *this;
  }

  constexpr FieldElementWrapper operator-() const {
    return Construct(kField.Negative(value_));
  }

  constexpr FieldElementWrapper& operator*=(const FieldElementWrapper& other) {
    value_ = kField.Multiply(value_, other.value_);
    return *this;
  }

  constexpr FieldElementWrapper& operator/=(const FieldElementWrapper& other) {
    value_ = kField.Divide(value_, other.value_);
    return *this;
  }

  constexpr FieldElementWrapper Inverse() const {
    return Construct(kField.Inverse(value_));
  }

  template <typename Power>
  constexpr FieldElementWrapper Pow(Power power) const {
    return Construct(kField.Pow(value_, power));
  }

  [[nodiscard]]
  constexpr static auto FieldBase() {
    return Field::FieldBase();
  }

  [[nodiscard]]
  constexpr static auto FieldPower() {
    return Field::FieldPower();
  }

  [[nodiscard]]
  constexpr static ElementsRange AllFieldElements() {
    return {};
  }

 private:
  template <std::convertible_to<Coefficient> T>
  constexpr std::array<Coefficient, kFieldPower> Convert(
      const std::array<T, kFieldPower>& value) {
    std::array<Coefficient, kFieldPower> result;
    for (size_t i = 0; i < kFieldPower; ++i) {
      result[i] = value[i];
    }
    return result;
  }

  constexpr static FieldElementWrapper Construct(Value value) {
    FieldElementWrapper result;
    result.value_ = value;
    return result;
  }

 private:
  constexpr static Field kField{};
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
