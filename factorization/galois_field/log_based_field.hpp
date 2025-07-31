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
#include <climits>
#include <cstdint>

#include <factorization/utils.hpp>

namespace factorization::galois_field {

/*! \brief Field implementation based on using logarithm tables.
 *
 *  @tparam kFieldBase Field characteristic
 *  @tparam kFieldPower Field power
 *  @tparam kFieldGenerator Primitive polynomial from lower power to greater
 *  @tparam Int Type used inside, uint32_t by default
 *
 *  Operates over elements written in polynomial form.
 *  All data is located on stack, can be constexpr.
 *  For fields with base 2 consumes O(n) space where n is field size
 *  For fields with other base consumes
 *    0(n * (2c)^k)
 *  where
 *    n is field size
 *    k is field power
 *    c is ratio of lowest power of 2 greater than field base to field base
 *  For example with field base 3 and power 2, it consumes
 *    O(n * (2 * 4 / 3)^2)
 *  Still better than use of add and multiply tables
 *  Operations have constant complexity.
 */
template <uint32_t kFieldBase,
          uint32_t kFieldPower,
          std::array<uint32_t, kFieldPower + 1> kFieldGenerator,
          std::integral Int = uint32_t>
class LogBasedField {
 public:
  using Value = Int;

 public:
  constexpr inline LogBasedField() {
    for (Int value = 0; value < (1 << kBitsPerSymbol * kFieldPower); ++value) {
      constexpr int kBitMask = (1 << kBitsPerSymbol) - 1;
      Int good_value = 0;
      for (uint32_t i = 0; i < kFieldPower; ++i) {
        // extract i-th symbol
        Int digit = (value >> (i * kBitsPerSymbol)) & kBitMask;
        // that means we got incorrect value - need to skip it 
        if (digit >= 2 * kFieldBase) {
          break;
        }
        digit = digit >= kFieldBase ? digit - kFieldBase : digit;
        good_value |= digit << (kBitsPerSymbol * i);
      }
      to_good_view_[value] = good_value;
    }
    // we have
    //   a[0] * alpha^0 + ... + a[k-1] * alpha^{k-1} + alpha^k = 0
    // thus
    //   alpha^k = -a[0] * alpha^0 - ... - a[k-1] * alpha^{k-1}
    Int polynom = 1;
    for (Int power = 0; power + 1 < kFieldSize; ++power) {
      log_to_poly_[power] = polynom;
      poly_to_log_[polynom] = power;
      // Now we want to multiply polynom by alpha
      // Consider
      //   polynom = c * alpha^{k - 1} + p(alpha), deg(p(alpha)) < k - 1
      // Then
      //   polynom * alpha = c * alpha^{k} + alpha * p(alpha)
      constexpr int kBitsShift = kBitsPerSymbol * (kFieldPower  - 1);
      Int overflow = polynom >> kBitsShift;
      // Remove c * alpha^{k - 1} at the beginning
      // and multiply remaining part by alpha
      polynom &= (1 << kBitsShift) - 1;
      polynom <<= kBitsPerSymbol;
      if (overflow > 0) {
        // now we want to add
        //   c * alpha^{k}, here c != 0
        // which is equal to
        //   c * (-a[0] * alpha^0 + ... + -a[k - 1] * alpha^{k - 1})
        // or simply
        //  -c * (a[0] * alpha^0 + ... + a[k - 1] * alpha^{k - 1})
        overflow = Negative(overflow);
        Int adder = 0;
        for (uint32_t i = kFieldPower - 1; i > 0; --i) {
          adder |= kFieldGenerator[i] * overflow % kFieldBase;
          adder <<= kBitsPerSymbol;
        }
        adder |= kFieldGenerator[0] * overflow % kFieldBase;
        polynom = Add(polynom, adder);
      }
    }
    for (uint32_t i = 0; i < kFieldSize - 1; ++i) {
      log_to_poly_[i + kFieldSize - 1] = log_to_poly_[i];
    }
  }

  constexpr inline Int Zero() const {
    return Int{0};
  }

  constexpr inline Int One() const {
    return Int{1};
  }

  //! Returns field element which is sum of 2 given.
  constexpr inline Int Add(Int first, Int second) const {
    return to_good_view_[first + second];
  }

  //! Returns field element that equal first - second
  //
  constexpr inline Int Sub(Int first, Int second) const {
    return Add(first, Negative(second));
  }

  //! returns element -value that -value + value = 0.
  constexpr inline Int Negative(Int value) const {
    constexpr Int kOtherZero = CalcMask() * kFieldBase;
    return to_good_view_[kOtherZero - value];
  }

  //! Returns field element which is product of 2 given.
  constexpr inline Int Multiply(Int first, Int second) const {
    if (first == 0 || second == 0) {
      return 0;
    }
    return log_to_poly_[poly_to_log_[first] + poly_to_log_[second]];
  }

  //! Returns field element that is result of multiply first by second**-1
  constexpr inline Int Divide(Int first, Int second) const {
    if (first == 0) {
      return 0;
    }
    return log_to_poly_[kFieldSize - 1 -
                        poly_to_log_[second] + poly_to_log_[first]];
  }

  //! Returns field element that equal given**power
  template <typename Power>
  constexpr inline Int Pow(Int base, Power power) const {
    if (base == 0) {
      return 0;
    }
    base = poly_to_log_[base];
    return log_to_poly_[power * base % (kFieldSize - 1)];
  }

  //! Returns field element b that ab = 1
  constexpr inline Int Inverse(Int value) const {
    return log_to_poly_[kFieldSize - 1 - poly_to_log_[value]];
  }

  //! Returns power that alpha^power = value
  constexpr inline Int Log(Int value) const {
    return poly_to_log_[value];
  }

  //! Returns field characteristic
  constexpr static inline uint32_t FieldBase() {
    return kFieldBase;
  }

  //! Returns field dimension  
  constexpr static inline uint32_t FieldPower() {
    return kFieldPower;
  }

  //! Return constant (if to use polynomial form)
  constexpr inline Int FieldValueFromConstant(Int value) const {
    return value % kFieldBase;
  }

  //! First to iterate over field values
  constexpr inline Int FirstFieldValue() const {
    return 0;
  }

  //! Next to iterate over field value
  // if value is last one behaviour is undefined
  // O(1*) - armotized constant - seems enough
  constexpr inline Int NextFieldValue(Int value) const {
    ++value;
    // here I use property
    // that correct field value will not change in this table 
    while (value != to_good_view_[value]) {
      ++value;
    }
    return value;
  }

  //! Last to iterate over field values
  // I hope it will be calculated during compilation
  constexpr inline Int LastFieldValue() const {
    return CalcMask() * (kFieldBase - 1);
  }

 private:
  constexpr static inline int GetBitesPerSymbol() {
    int bits_per_symbol = 1;
    while ((1 << bits_per_symbol) < kFieldBase) {
      ++bits_per_symbol;
    }
    return bits_per_symbol + 1;
  }

  constexpr static inline Int CalcMask() {
    Int result = 1;
    for (uint32_t i = 1; i < kFieldPower; ++i) {
      result = (result << kBitsPerSymbol) + 1;
    }
    return result;
  }

  constexpr static inline int kFieldSize = utils::BinPow(kFieldBase, kFieldPower);
  constexpr static inline int kBitsPerSymbol = GetBitesPerSymbol();

  static_assert(kBitsPerSymbol * kFieldPower <= sizeof(Int) * CHAR_BIT,
                "Integer type too small for this field");

  std::array<Int, 2 * kFieldSize> log_to_poly_{};
  // many values of this array will be inconsistent
  std::array<Int, 1 << (kBitsPerSymbol * kFieldPower - 1)> poly_to_log_{};
  // values in polynomial form are presented as
  //   a(x) = a_0 + a_1 x + ...
  // usually x = kFieldBase is used to store a, but it make addition difficult
  // here i use x = 2^kBitsPerSymbol, addition becomes like this
  //  (a_0 + b_0) + (a_1 + b_1) x + ...
  // Note: a_0 + b_0 MUST be less than x
  // then we know that a_i + b_i goes to (a_i + b_i) (mod kFieldBase)
  // so I precalculate outcome for each possible a + b
  // Addition becomes much faster
  // some values here can be inconsistent
  std::array<Int, 1 << (kBitsPerSymbol * kFieldPower)> to_good_view_{};
};

/** Specialization for base equal to 2.
 *  Addition can be replaced with XOR which much faster.
 */
template <uint32_t kFieldPower,
          std::array<uint32_t, kFieldPower + 1> kFieldGenerator,
          std::integral Int>
class LogBasedField<2, kFieldPower, kFieldGenerator, Int> {
 public:
  using Value = Int;

 public:
  constexpr inline LogBasedField() {
    // Simple version of above
    Int generator = kFieldGenerator[0];
    Int alpha = 1;
    for (uint32_t i = 1; i < kFieldPower; ++i) {
      alpha <<= 1;
      generator ^= alpha * kFieldGenerator[i];
    }

    Int polynom = 1;
    for (Int power = 0; power + 1 < kFieldSize; ++power) {
      log_to_poly_[power] = polynom;
      poly_to_log_[polynom] = power;

      polynom = polynom >= alpha 
                  ? Add((polynom - alpha) << 1, generator)
                  : polynom << 1;
    }
    for (uint32_t i = 0; i < kFieldSize; ++i) {
      log_to_poly_[i + kFieldSize - 1] = log_to_poly_[i];
    }
  }

  constexpr inline Int Zero() const {
    return Int{0};
  }

  constexpr inline Int One() const {
    return Int{1};
  }

  constexpr inline Int Add(Int first, Int second) const {
    return first ^ second;
  }

  constexpr inline Int Sub(Int first, Int second) const {
    return first ^ second;
  }

  constexpr inline Int Negative(Int value) const {
    return value;
  }

  constexpr inline Int Multiply(Int first, Int second) const {
    if (first == 0 || second == 0) {
      return 0;
    }
    return log_to_poly_[poly_to_log_[first] + poly_to_log_[second]];
  }

  constexpr inline Int Divide(Int first, Int second) const {
    if (first == 0) {
      return 0;
    }
    return log_to_poly_[kFieldSize - 1 -
                        poly_to_log_[second] + poly_to_log_[first]];
  }

  template <typename Power>
  constexpr inline Int Pow(Int base, Power power) const {
    if (base == 0) {
      return 0;
    }
    base = poly_to_log_[base];
    return log_to_poly_[power * base % (kFieldSize - 1)];
  }

  constexpr inline Int Inverse(Int value) const {
    if (value == 1) {
      return value;
    }
    return log_to_poly_[kFieldSize - 1 - poly_to_log_[value]];
  }

  constexpr inline Int Log(Int value) const {
    return poly_to_log_[value];
  }

  constexpr static inline uint32_t FieldBase() {
    return 2;
  }

  constexpr static inline uint32_t FieldPower() {
    return kFieldPower;
  }

  constexpr inline Int FieldValueFromConstant(Int value) const {
    return value & 1;
  }

  constexpr inline Int FirstFieldValue() const {
    return 0;
  }

  constexpr inline Int NextFieldValue(Int value) const {
    return ++value;
  }

  constexpr inline Int LastFieldValue() const {
    return kFieldSize - 1;
  }

 private:
  constexpr static inline uint32_t kFieldSize = 1u << kFieldPower;
  // can it be done constexpr since we have constexpr constructor?
  std::array<Int, 2 * kFieldSize> log_to_poly_{};
  std::array<Int, kFieldSize> poly_to_log_{};
};

}  // namespace factorization::galois_field