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
 *  Has linear construction complexity, takes linear amount of space. 
 *  Most of operations have constant complexity.
 */
template <uint32_t kFieldBase,
          uint32_t kFieldPower,
          std::array<uint32_t, kFieldPower + 1> kFieldGenerator,
          std::integral Int = uint32_t>
class LogBasedField {
 public:
  using Value = Int;

 public:
  constexpr LogBasedField() {
    // these 2 forms are impossible
    poly_to_log_[0] = kFieldSize - 1;
    log_to_poly_[kFieldSize - 1] = 0;
    // we have
    //   a[0] * alpha^0 + ... + a[k-1] * alpha^{k-1} + alpha^k = 0
    // thus
    //   alpha^k = -a[0] * alpha^0 - ... - a[k-1] * alpha^{k-1}
    Int alpha = utils::BinPow(kFieldBase, kFieldPower - 1);
    Int polynom = 1;
    for (Int power = 0; power + 1 < kFieldSize; ++power) {
      log_to_poly_[power] = polynom;
      poly_to_log_[polynom] = power;
      // Now we want to multiply polynom by alpha
      // Consider
      //   polynom = c * alpha^{k - 1} + p(alpha), deg(p(alpha)) < k - 1
      // Then
      //   polynom * alpha = c * alpha^{k} + alpha * p(alpha)
      Int overflow = polynom / alpha;
      // Remove c * alpha^{k - 1} at the beginning
      // and multiply remaining part by alpha
      polynom = (polynom - overflow * alpha) * kFieldBase;
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
          adder += kFieldGenerator[i] * overflow % kFieldBase;
          adder *= kFieldBase;
        }
        adder += kFieldGenerator[0] * overflow % kFieldBase;
        polynom = Add(polynom, adder);
      }
    }
  }

  //! Returns field element which is sum of 2 given.
  constexpr Int Add(Int first, Int second) const {
    // using polynomial form:
    //   first  = a[0] * alpha^0 + ... + a[k] * alpha^k
    //   second = b[0] * alpha^0 + ... + b[k] * alpha^k
    // we add these 2 polynomials element by element
    Int result{(first % kFieldBase + second % kFieldBase) % kFieldBase};
    // done a little dirty in terms of appearance
    Int alpha = kFieldBase;
    first /= kFieldBase;
    second /= kFieldBase;
    while (first != 0 || second != 0) {
      result += (first + second) % kFieldBase * alpha;
      alpha *= kFieldBase;
      first /= kFieldBase;
      second /= kFieldBase;
    }
    return result;
  }

  //! returns element -value that -value + value = 0.
  constexpr Int Negative(Int value) {
    Int result = 0;
    Int alpha = 1;
    while (value != 0) {
      auto a = value % kFieldBase;
      result += a != 0 ? (kFieldBase - a) * alpha : 0;
      value /= kFieldBase;
      alpha *= kFieldBase;
    }
    return result;
  }

  //! Returns field element which is product of 2 given.
  constexpr Int Multiply(Int first, Int second) const {
    if (first == 0 || second == 0) {
      return 0;
    }
    first = poly_to_log_[first];
    second = poly_to_log_[second];
    return log_to_poly_[(first + second) % (kFieldSize - 1)];
  }

  //! Returns field element that equal given**power
  template <typename Power>
  constexpr Int Pow(Int base, Power power) const {
    if (base == 0) {
      return 0;
    }
    base = poly_to_log_[base];
    return log_to_poly_[power * base % (kFieldSize - 1)];
  }

  //! Returns field element b that ab = 1
  constexpr Int Inverse(Int value) const {
    if (value == 1) {
      return value;
    }
    value = poly_to_log_[value];
    return log_to_poly_[kFieldSize - 1 - value];
  }

  //! Returns power that alpha^power = value
  constexpr Int Log(Int value) const {
    return poly_to_log_[value];
  }

 private:
  constexpr static uint32_t kFieldSize = utils::BinPow(kFieldBase, kFieldPower);
  // can it be done constexpr since we have constexpr constructor?
  std::array<Int, kFieldSize> log_to_poly_;
  std::array<Int, kFieldSize> poly_to_log_;
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
  constexpr LogBasedField() {
    // these 2 forms are impossible
    poly_to_log_[0] = kFieldSize - 1;
    log_to_poly_[kFieldSize - 1] = 0;
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
  }

  constexpr Int Add(Int first, Int second) const {
    return first ^ second;
  }

  constexpr Int Negative(Int value) {
    return value;
  }

  constexpr Int Multiply(Int first, Int second) const {
    if (first == 0 || second == 0) {
      return 0;
    }
    first = poly_to_log_[first];
    second = poly_to_log_[second];
    return log_to_poly_[(first + second) % (kFieldSize - 1)];
  }

  template <typename Power>
  constexpr Int Pow(Int base, Power power) const {
    if (base == 0) {
      return 0;
    }
    base = poly_to_log_[base];
    return log_to_poly_[power * base % (kFieldSize - 1)];
  }

  constexpr Int Inverse(Int value) const {
    if (value == 1) {
      return value;
    }
    value = poly_to_log_[value];
    return log_to_poly_[kFieldSize - 1 - value];
  }

  constexpr Int Log(Int value) const {
    return poly_to_log_[value];
  }

 private:
  constexpr static uint32_t kFieldSize = 1u << kFieldPower;
  // can it be done constexpr since we have constexpr constructor?
  std::array<Int, kFieldSize> log_to_poly_;
  std::array<Int, kFieldSize> poly_to_log_;
};

}  // namespace factorization::galois_field