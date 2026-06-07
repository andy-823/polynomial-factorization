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

/*! @brief Extracts the p-th root of a polynomial known to be a p-th power.
 *
 *  Let p be the field characteristic. If
 *    f(x) = g(x)^p,
 *  then only powers x^{pi} can have nonzero coefficients in f.
 *  That means
 *    f(x) = f_0 + f_p x^p + f_{2p} x^{2p} + ...
 *  and
 *    g(x) = f_0^{1/p} + f_p^{1/p} x + f_{2p}^{1/p} x^2 + ... .
 *
 *  @pre polynom is a p-th power.
 */
template <concepts::Polynom Polynom>
inline Polynom FieldBaseRoot(const Polynom& polynom) {
  using Element = typename Polynom::Element;
  // Use auto because field sizes may exceed a fixed small integer type.
  constexpr auto kFieldBase = Element::FieldBase();
  constexpr auto kFieldPower = Element::FieldPower();

  std::vector<Element> elements(polynom.Get());
  for (size_t i = 0; i < elements.size(); i += kFieldBase) {
    // If q = p^k and y^p = x, then
    //   y = x^{p^{k-1}}.
    constexpr auto kPower = utils::BinPow(kFieldBase, kFieldPower - 1);
    elements[i / kFieldBase] = elements[i].Pow(kPower);
  }
  elements.resize((elements.size() + kFieldBase - 1) / kFieldBase);
  return Polynom(std::move(elements));
}

/*! @brief Decomposes a polynomial into square-free factors.
 *
 *  The result contains pairs (g_i, i), where each g_i is square-free and is the
 *  product of all irreducible factors that occur in the input with multiplicity
 *  exactly i.
 *
 *  @pre polynom is nonzero.
 */
template <concepts::Polynom Polynom>
constexpr inline std::vector<solver::Factor<Polynom>> SquareFreeFactorize(
    Polynom polynom) {
  std::vector<solver::Factor<Polynom>> result;
  // If f = g_1 g_2^2 ... g_m^m h^p, then
  //   gcd(f, f') = g_2 g_3^2 ... g_m^{m-1} h^p.
  Polynom c = polynom.Gcd(polynom.Derivative());
  // factors = g_1 g_2 ... g_m.
  Polynom factors = polynom.Div(c);
  int j = 1;
  while (!factors.IsOne()) {
    // At step j:
    //   factors = g_j g_{j+1} ...
    //   c contains g_{j+1} g_{j+2}^2 ... g_m^{m - j} h^p.
    // So next_factors is
    //   g_{j+1} ... g_m,
    // and g_j = factors / next_factors.
    Polynom next_factors = factors.Gcd(c);
    if (factors != next_factors) {
      auto factor = std::move(factors).Div(next_factors);
      result.emplace_back(factor, j);
    }
    factors = std::move(next_factors);
    c = std::move(c).Div(factors);
    ++j;
  }
  // Only h^p is left.
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
