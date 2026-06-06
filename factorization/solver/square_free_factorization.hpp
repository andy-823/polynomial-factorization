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

#include <cstdint>
#include <vector>

#include <factorization/concepts.hpp>

namespace factorization::sff {

// assume polynom is
//   f(x) = g(x)^p
// this function returns such g(x) by given f(x)
template <concepts::Polynom Polynom>
inline Polynom FieldBaseRoot(const Polynom& polynom) {
  using Element = typename Polynom::Element;
  // use auto because in some cases field can be big
  constexpr auto kFieldBase = Element::FieldBase();
  constexpr auto kFieldPower = Element::FieldPower();

  std::vector<Element> elements(polynom.Get());
  for (size_t i = 0; i < elements.size(); i += kFieldBase) {
    // want to get y that y^p = x
    // consider p is field base, q = p^k is field size
    // then y = x^{q/p} = x^p^{k - 1}
    constexpr auto kPower = utils::BinPow(kFieldBase, kFieldPower - 1);
    elements[i / kFieldBase] = elements[i].Pow(kPower);
  }
  elements.resize((elements.size() + kFieldBase - 1) / kFieldBase);
  return Polynom(std::move(elements));
}

template <concepts::Polynom Polynom>
constexpr inline std::vector<solver::Factor<Polynom>> SquareFreeFactorize(
    Polynom polynom) {
  std::vector<solver::Factor<Polynom>> result;
  // f  = g_1 g_2^2 ... g_m^m * h^p
  // f' = g_2 g_m^{m-1} * h^p * (sum i g'_i * g / g_i)
  // c = g_2^1 ... g_m^{m-1} h^p
  Polynom c = polynom.Gcd(polynom.Derivative());
  // factors = g_1 ... g_m
  Polynom factors = polynom.Div(c);
  int j = 1;
  // on each iteration we extract one factor
  // from square_free_factors
  while (!factors.IsOne()) {
    // next_factors = g_2 ... g_m
    Polynom next_factors = factors.Gcd(c);
    if (factors != next_factors) {
      auto factor = std::move(factors).Div(next_factors);
      result.emplace_back(factor, j);
    }
    factors = std::move(next_factors);
    c = std::move(c).Div(factors);
    ++j;
  }
  // that means polynom is p-th power
  if (!c.IsOne()) {
    auto factors = SquareFreeFactorize(FieldBaseRoot(c));
    for (const auto& [factor, power] : factors) {
      auto factor_power = Polynom::Element::FieldBase() * power;
      result.emplace_back(factor, factor_power);
    }
  }
  return result;
}

}  // namespace factorization::sff