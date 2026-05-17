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

#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/utils.hpp>

namespace factorization::ddf {

template <concepts::Polynom Polynom>
struct DistinctDegreeFactor {
  Polynom factor;
  int degree;
};

namespace naive {

template <concepts::Polynom Poly>
std::vector<DistinctDegreeFactor<Poly>> DistinctDegreeFactorize(Poly poly) {
  using Element = typename Poly::Element;
  constexpr auto kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());

  poly = std::move(poly).MakeMonic();
  if (poly.IsZero() || poly.IsOne()) {
    return {};
  }
  const Poly x(std::vector<Element>{Element::Zero(), Element::One()});

  std::vector<DistinctDegreeFactor<Poly>> result;

  Poly h = x;
  size_t degree = 1;
  auto mod = poly.BuildModulus(2 * poly.Size());

  while (2 * degree <= poly.Size() - 1) {
    h = polynomial::BinPowMod(std::move(h), kFieldSize, mod);
    Poly factor = poly.Gcd(h.Sub(x));
    if (!factor.IsOne()) {
      result.emplace_back(factor, static_cast<int>(degree));
      poly = std::move(poly).Div(factor).MakeMonic();
      // factorization is complete
      if (poly.IsOne()) {
        break;
      }
      mod = poly.BuildModulus(2 * poly.Size());
      h = std::move(h).Rem(mod);
    }
    ++degree;
  }
  if (!poly.IsOne()) {
    result.emplace_back(std::move(poly), static_cast<int>(poly.Size() - 1));
  }
  return result;
}

}  // namespace naive

}  // namespace factorization::ddf
