#pragma once

#include <array>
#include <cstdint>

#include <factorization/utils.hpp>

namespace factorization::internal {

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
    Int generator = 0;
    Int alpha = 1;
    for (uint32_t i = 0; i < kFieldPower; ++i) {
      generator += alpha * kFieldGenerator[i];
      alpha *= kFieldBase;
    }
    generator = Negative(generator);

    Int polynom = 1;
    for (Int power = 0; power < kFieldSize; ++power) {
      log_to_poly_[power] = polynom;
      poly_to_log_[polynom] = power;

      polynom *= kFieldBase;
      if (polynom >= alpha) {
        Int overflow = polynom / alpha;
        polynom = Add(polynom, overflow * generator);
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
  constexpr static uint32_t kFieldSize = BinPow(kFieldBase, kFieldPower);
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
  constexpr LogBasedField() {
    // these 2 forms are impossible
    poly_to_log_[0] = kFieldSize - 1;
    log_to_poly_[kFieldSize - 1] = 0;
    // Simple version of above
    Int generator = 0;
    Int alpha = 1;
    for (uint32_t i = 0; i < kFieldPower; ++i) {
      generator ^= alpha * kFieldGenerator[i];
      alpha <<= 1;
    }

    Int polynom = 1;
    for (Int power = 0; power < kFieldSize; ++power) {
      log_to_poly_[power] = polynom;
      poly_to_log_[polynom] = power;

      polynom <<= 1;
      if (polynom >= alpha) {
        polynom = Add(polynom, generator);
      }
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

}  // namespace factorization::internal