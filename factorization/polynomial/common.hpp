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

/*! @brief Computes poly^pow modulo mod by binary exponentiation.
 *
 *  @pre poly is already reduced modulo mod.
 *  @pre pow is non-negative.
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

/*! @brief Builds a power table for modular composition.
 *
 *  For a composition argument h, the returned matrix stores the coefficient
 *  vectors of
 *    1, h, h^2, ..., h^t (mod f),
 *  where f is represented by mod. The last row, h^t (mod f),
 *  is later used as the step in Horner's rule.
 */
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

/*! @brief Adapts a composition table to a smaller modulus.
 *
 *  Suppose the table contains
 *    1, g, g^2, ..., g^s (mod f).
 *
 *  If the current modulus changes from f to its divisor f*, the same powers of
 *  g can be reused: each row only has to be reduced modulo f*. This function
 *  keeps the first t + 1 rows and updates them in place.
 *
 *  @pre matrix was produced by BuildCompModMatrix for the same g and an old
 *       modulus divisible by the new one.
 *  @pre matrix.size() >= t + 1.
 */
template <concepts::Polynom Poly>
void UpdateCompModMatrix(
    std::vector<std::vector<typename Poly::Element>>& matrix, size_t t,
    const typename Poly::Modulus& mod) {
  assert(t > 0);
  assert(matrix.size() >= t + 1);

  matrix.resize(t + 1);
  for (auto& row : matrix) {
    row = Poly(std::move(row)).Rem(mod).Get();
  }
}

/*! @brief Computes poly(h) (mod f) using a precomputed power table.
 *
 *  The matrix must be produced by BuildCompModMatrix(h, t, mod) or updated by
 *  UpdateCompModMatrix for the same h.
 *
 *  The algorithm writes
 *    poly(x) = g_0(x) + g_1(x) x^t + ...,
 *  where every block g_j has degree less than t.
 *
 *  First, the rows
 *    1, h, ..., h^{t-1}
 *  are used to compute all values
 *    g_j(h) (mod f).
 *  Then these values are combined by Horner's rule with
 *    h^t (mod f),
 *  which is stored as the last row of the table.
 */
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

  // Write
  //   poly(x) = c_0 x^0 + ... + c_{n-1} x^{n-1}
  //           = g_0(x) + g_1(x) x^t + ... .
  // Values g_i(h) (mod f) are computed by the matrix product
  //   [ g_0                 ]   [ h^0     (mod f) ]
  //   [ ...                 ] * [ ...             ] = G * H.
  //   [ g_{block_count - 1} ]   [ h^{t-1} (mod f) ]
  // Polynomials in the matrices above are coefficient row vectors.
  // G is not computed directly, H is stored in the precomputed table.
  std::vector<Poly> blocks(block_count);
  for (size_t block_id = 0; block_id < block_count; ++block_id) {
    std::vector<Element> block;
    const size_t from = block_id * t;
    // The block row is
    //   g_{block_id}(x) = c_{from} + ... + c_{from + t - 1} x^{t-1}.
    // The last block may contain fewer than t coefficients.
    // This iteration computes
    //   g_{block_id}(h) (mod f) = g_{block_id} * H.
    // Result is stored in variable `block`.
    // The loop order is chosen to be more cache-friendly.
    for (size_t i = 0; i < t && from + i < n; ++i) {
      const auto& c = coefficients[from + i];
      if (c == Element::Zero()) {
        continue;
      }
      // row is h^i (mod f).
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

  // Starting from the highest block, maintain
  //   result <- result * h^t + g_block(h) (mod f).
  // In the end we have
  //   result = g_0(h) + g_1(h) h^t + ... (mod f).
  const Poly block_step(matrix.back());
  Poly result;
  for (size_t block = block_count; block-- > 0;) {
    if (!result.IsZero()) {
      result = std::move(result).Mul(block_step).Rem(mod);
    }
    result = std::move(result).Add(blocks[block]);
  }
  return result;
}

}  // namespace factorization::polynomial
