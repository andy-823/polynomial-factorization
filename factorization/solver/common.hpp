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

namespace factorization::solver {

template <concepts::Polynom Polynom>
struct Factor {
  Polynom factor;
  int power;

  bool operator==(const Factor&) const = default;
};

template <concepts::Polynom Polynom>
constexpr inline Polynom Gcd(Polynom first, Polynom second) {
  while (!second.IsZero()) {
    first %= second;
    std::swap(first, second);
  }
  first.MakeMonic();
  return first;
}

// assume polynom is
//   f(x) = g(x)^p
// this function returns such g(x) by given f(x)
template <concepts::Polynom Polynom>
inline Polynom FieldBaseRoot(const Polynom& polynom) {
  using Element = typename Polynom::Element;
  // use auto because in some cases field can be big
  constexpr auto kFieldBase = Element::FieldBase();
  constexpr auto kFieldPower = Element::FieldPower();

  std::vector<Element> elements(polynom.GetElements());
  for (size_t i = 0; i < elements.size(); i += kFieldBase) {
    // want to get y that y^p = x
    // consider p is field base, q = p^k is field size
    // then y = x^{q/p} = x^p^{k - 1}
    constexpr auto kPower = utils::BinPow(kFieldBase,
                                          kFieldPower - 1);
    elements[i / kFieldBase] = elements[i].Pow(kPower);
  }
  elements.resize((elements.size() + kFieldBase - 1) / kFieldBase);
  return Polynom(std::move(elements));
}

template <concepts::Polynom Polynom>
constexpr inline std::vector<Factor<Polynom>> SquareFreeFactorize(
    Polynom polynom) {
  std::vector<Factor<Polynom>> result;
  // f  = g_1 g_2^2 ... g_m^m * h^p
  // f' = g_2 g_m^{m-1} * h^p * (sum i g'_i * g / g_i)
  Polynom derivative = polynom.Derivative();
  // c = g_2^2 ... g_m^{m-1} h^p
  Polynom c = Gcd(polynom, derivative);
  // w = g_1 ... g_m
  Polynom square_free_factors = polynom / c;
  // d = sum (i - 1) * g'_i * g / g_i
  Polynom d = derivative / c - square_free_factors.Derivative();
  int j = 1;
  // on each iteration we extract one factor
  // from square_free_factors
  while (!square_free_factors.IsOne()) {
    polynom /= square_free_factors;
    auto factor = Gcd(square_free_factors, d);
    result.emplace_back(factor, j);
    square_free_factors /= factor;
    d = d / factor - square_free_factors.Derivative();
    ++j;
  }
  // that means polynom is p-th power
  if (!polynom.IsOne()) {
    auto factors = SquareFreeFactorize(FieldBaseRoot(polynom));
    for (const auto& [factor, power] : factors) {
      int factor_power = power * Polynom::Element::FieldBase();
      result.emplace_back(factor, factor_power);
    }
  }
  return result;
}

}  // namespace factorization::solver