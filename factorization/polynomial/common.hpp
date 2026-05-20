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

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>

namespace factorization::polynomial {

/* Performs binary powering
 * Assume that poly.Rem(mod) = poly
 */
template <concepts::Polynom Poly, std::integral Int>
Poly BinPowMod(Poly poly, Int pow, const typename Poly::Modulus& mod) {
  Poly result(Poly::Element::One());
  while (pow > 0) {
    if (pow % 2 != 0) {
      result = std::move(result).Mul(poly).Rem(mod);
    }
    pow /= 2;
    if (pow > 0) {
      poly = poly.Mul(poly).Rem(mod);
    }
  }
  return result;
}

template <concepts::Polynom Poly>
std::vector<std::vector<typename Poly::Element>> BuildCompModMatrix(
    const Poly& poly, size_t t, const typename Poly::Modulus& mod) {
  assert(t > 0);

  using Element = typename Poly::Element;
  std::vector<std::vector<Element>> matrix;
  matrix.reserve(t + 1);

  matrix.emplace_back(Poly(Element::One()).Get());
  for (size_t i = 1; i <= t; ++i) {
    Poly current = Poly(matrix.back()).Mul(poly).Rem(mod);
    matrix.emplace_back(std::move(current).Get());
  }
  return matrix;
}

template <concepts::Polynom Poly>
std::vector<std::vector<typename Poly::Element>> UpdateCompModMatrix(
    std::vector<std::vector<typename Poly::Element>>& matrix, size_t t,
    const typename Poly::Modulus& mod) {
  assert(t > 0);
  assert(matrix.size() >= t + 1);

  matrix.resize(t + 1);
  for (auto& row : matrix) {
    row = Poly(std::move(row)).Rem(mod).Get();
  }
  return matrix;
}

template <concepts::Polynom Poly>
Poly CompMod(const Poly& poly,
             const std::vector<std::vector<typename Poly::Element>>& matrix,
             const typename Poly::Modulus& mod) {
  assert(matrix.size() >= 2);
  using Element = typename Poly::Element;

  const auto coefficients = poly.Get();
  if (coefficients.empty()) {
    return Poly();
  }
  const size_t t = matrix.size() - 1;
  const size_t n = coefficients.size();
  const size_t block_count = (n + t - 1) / t;

  std::vector<Poly> blocks(block_count);
  for (size_t block_id = 0; block_id < block_count; ++block_id) {
    std::vector<Element> block;
    const size_t from = block_id * t;
    for (size_t i = 0; i < t && from + i < n; ++i) {
      const auto& c = coefficients[from + i];
      if (c == Element::Zero()) {
        continue;
      }
      const auto& row = matrix[i];
      if (block.size() < row.size()) {
        block.resize(row.size(), Element::Zero());
      }
      for (size_t j = 0; j < row.size(); ++j) {
        block[j] += c * row[j];
      }
    }
    blocks[block_id] = Poly(std::move(block));
  }

  const Poly giant_step(matrix.back());
  Poly result;
  for (size_t block = block_count; block-- > 0;) {
    if (!result.IsZero()) {
      result = std::move(result).Mul(giant_step).Rem(mod);
    }
    result = std::move(result).Add(blocks[block]);
  }
  return result;
}

}  // namespace factorization::polynomial
