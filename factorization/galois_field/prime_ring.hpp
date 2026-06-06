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
#include <concepts>
#include <cstdint>

namespace factorization::galois_field {

/*! \brief Field implementation of Z_p
 */
template <uint64_t kFieldBase, std::integral Int = uint32_t,
          std::integral DoubleInt = uint64_t>
class PrimeRing {
 public:
  using Value = Int;
  using Coefficient = Int;

 public:
  constexpr PrimeRing() {
  }

  constexpr Int Encode(const std::array<Coefficient, 1>& arr) const {
    return arr[0] % kFieldBase;
  }

  constexpr Int Encode(Coefficient value) const {
    return value % kFieldBase;
  }

  constexpr std::array<Coefficient, 1> Decode(Int value) const {
    return {value};
  }

  constexpr Int Zero() const {
    return Int{0};
  }

  constexpr Int One() const {
    return Int{1};
  }

  //! Returns field element which is sum of 2 given.
  constexpr Int Add(Int first, Int second) const {
    Int result = first + second;
    return result >= kFieldBase ? result - kFieldBase : result;
  }

  //! Returns field element that equal first - second
  constexpr Int Sub(Int first, Int second) const {
    if (first >= second) {
      return first - second;
    }
    return kFieldBase - second + first;
  }

  //! returns element -value that -value + value = 0.
  constexpr Int Negative(Int value) const {
    return value != 0 ? kFieldBase - value : value;
  }

  //! Returns field element which is product of 2 given.
  constexpr Int Multiply(Int first, Int second) const {
    return static_cast<DoubleInt>(first) * second % kFieldBase;
  }

  //! Returns field element that is result of multiply first by second**-1
  constexpr Int Divide(Int first, Int second) const {
    if (first == 0) {
      return 0;
    }
    return Multiply(first, Inverse(second));
  }

  //! Returns field element that equal given**power
  template <typename Power>
  constexpr Int Pow(Int base, Power power) const {
    Int result = One();
    while (power > 0) {
      if (power % 2 != 0) {
        result = Multiply(result, base);
      }
      base = Multiply(base, base);
      power /= 2;
    }
    return result;
  }

  //! Returns field element b that ab = 1
  constexpr Int Inverse(Int value) const {
    return Pow(value, kFieldBase - 2);
  }

  //! Returns field characteristic
  constexpr static Value FieldBase() {
    return kFieldBase;
  }

  //! Returns field dimension
  constexpr static uint32_t FieldPower() {
    return 1;
  }
};

}  // namespace factorization::galois_field
