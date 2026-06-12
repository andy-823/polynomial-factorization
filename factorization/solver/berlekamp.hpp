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

#include <cstddef>
#include <utility>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

#include "common.hpp"
#include "square_free_factorization.hpp"

namespace factorization::solver {

/*! @brief Berlekamp factorization over a finite field.
 *
 *  The public entry point first removes repeated factors with square-free
 *  factorization. Each square-free component is then split into irreducible
 *  factors using the Berlekamp subalgebra.
 */
template <concepts::Polynom Polynom>
class Berlekamp {
  using Element = typename Polynom::Element;

 public:
  /*! @brief Factorizes a polynomial into irreducible factors with powers. */
  inline std::vector<Factor<Polynom>> Factorize(Polynom polynom) const {
    std::vector<Factor<Polynom>> result;
    polynom = std::move(polynom).MakeMonic();
    if (polynom.IsZero() || polynom.IsOne()) {
      return {};
    }
    for (const auto& [sf_factor, power] : sff::SquareFreeFactorize(polynom)) {
      for (const auto& factor : FactorizeImpl(sf_factor)) {
        result.emplace_back(factor, power);
      }
    }
    return result;
  }

 private:
  /*! @brief Splits a monic square-free polynomial into irreducible factors.
   *
   *  Input has the form
   *    f = f_1 ... f_k,
   *  where all f_i are distinct irreducible factors.
   *
   *  An element b of the Berlekamp subalgebra is constant modulo every
   *  f_i. Therefore, for each field element c, the polynomial
   *    gcd(f, b - c)
   *  collects exactly those irreducible factors f_i for which
   *    b = c (mod f_i).
   *  Applying this splitting to enough basis elements separates all f_i.
   */
  inline std::vector<Polynom> FactorizeImpl(Polynom polynom) const {
    std::vector<Polynom> basis = FindFactorizingBasis(polynom);
    // The basis size equals the number of irreducible factors. If it is one,
    // then f is irreducible.
    if (basis.size() == 1) {
      return {polynom};
    }

    const auto field_elements = Element::AllFieldElements();
    // factors is the current partition of f. Each nonconstant basis element
    // refines every part and writes the refined partition to new_factors.
    std::vector<Polynom> factors = {polynom};
    std::vector<Polynom> new_factors;
    new_factors.reserve(basis.size());

    for (const auto& factorizing : basis) {
      if (factorizing.Size() == 1) {  // constant
        continue;
      }
      for (const auto& factor : factors) {
        // For any current divisor of f and any Berlekamp subalgebra element b,
        //   factor = gcd(factor, b - c_1) * ... * gcd(factor, b - c_q),
        // where c_1, ..., c_q are all field elements. Thus the loop over c
        // splits factor without losing any irreducible divisor.
        for (const auto& c : field_elements) {
          Polynom new_factor = factor.Gcd(factorizing.Sub(c));
          if (!new_factor.IsOne()) {
            new_factors.emplace_back(std::move(new_factor));
          }
          // All irreducible factors have been found.
          if (new_factors.size() == basis.size()) {
            return new_factors;
          }
        }
      }
      factors.swap(new_factors);
      new_factors.clear();
    }
    return factors;
  }

  /*! @brief Finds a basis of the Berlekamp subalgebra.
   *
   *  For an input polynomial f over GF(q), the subalgebra consists of
   *  polynomials g such that
   *    g^q = g (mod f).
   *  Since the map g -> g^q (mod f) is linear over GF(q), the basis is found as
   *  the kernel of A - I, where A is the matrix of this linear map.
   */
  inline std::vector<Polynom> FindFactorizingBasis(
      const Polynom& polynom) const {
    // We want to solve Ax = 0.
    std::vector<Polynom> result;
    // First, build the matrix described above.
    auto matrix = BuildMatrix(polynom);
    // Put it into reduced row echelon form.
    matrix = PerformGaussElimination(std::move(matrix));
    // Then extract the kernel basis vectors.
    size_t rank = matrix.size();
    size_t n = polynom.Size() - 1;
    // Columns that correspond to free variables in the linear system.
    std::vector<size_t> free_coefficient_positions;
    // data_position[row] is the pivot column in this row.
    std::vector<size_t> data_position;
    free_coefficient_positions.reserve(n - rank);
    data_position.reserve(rank);
    result.reserve(n - rank);

    {
      size_t column = 0;
      for (size_t row = 0; row < rank; ++row) {
        while (column < n && matrix[row][column] == Element::Zero()) {
          free_coefficient_positions.emplace_back(column);
          ++column;
        }
        data_position.emplace_back(column);
        ++column;
      }
      // Columns after the last pivot are also free variables.
      while (column < n) {
        free_coefficient_positions.emplace_back(column);
        ++column;
      }
    }

    // The idea used here is the following. After elimination, we know pivot
    // columns and free coefficient positions. For example, for matrix
    //   0 1 1 0
    //   0 0 0 1
    // free positions are 0 and 2, while pivot columns are 1 and 3.
    //
    // To obtain a basis of solutions to Ax = 0, we process free positions one
    // by one. For the current free position, set it to one and set all other
    // free positions to zero; then restore pivot coefficients from the row
    // equations. In the example above, setting x_2 = 1 gives
    //   x_1 + x_2 = 0,
    //   x_3 = 0,
    // hence the vector
    //   (0, -1, 1, 0).
    for (const auto& column : free_coefficient_positions) {
      std::vector<Element> current(n, Element::Zero());
      // Choose the current free coefficient. All other free coefficients remain
      // zero because current was initialized with zeros.
      current[column] = Element::One();
      for (size_t row = 0; row < rank; ++row) {
        current[data_position[row]] = -matrix[row][column];
      }
      result.emplace_back(std::move(current));
    }
    return result;
  }

  /*! @brief Builds (A - I)^T for the map y -> y^q (mod f).
   *
   *  Coefficients are treated as row vectors. The matrix A is defined by
   *    y A = y^q (mod f).
   *  To solve y(A - I) = 0 with ordinary column-vector Gaussian elimination, the
   *  function returns (A - I)^T.
   */
  inline std::vector<std::vector<Element>> BuildMatrix(
      const Polynom& factorizing) const {
    constexpr int kFieldSize =
        utils::BinPow(Element::FieldBase(), Element::FieldPower());
    size_t n = factorizing.Size() - 1;
    std::vector<std::vector<Element>> result(n, std::vector<Element>(n));
    // Rows of A are images of basis monomials:
    //   A_0     = x^{0 * q} (mod f)
    //   A_1     = x^{1 * q} (mod f)
    //   ...
    //   A_{n-1} = x^{(n - 1) * q} (mod f)
    {
      // base = x^q (mod f).
      Polynom base;
      Polynom current(Element::One());
      {
        // Construct x^q directly and reduce it modulo f.
        std::vector<Element> tmp(kFieldSize + 1);
        tmp.back() = Element::One();  // the only nonzero element is last
        base = Polynom(std::move(tmp)).Rem(factorizing);
      }
      for (size_t power = 0; power < n; ++power) {
        auto elems = current.Get();
        for (size_t i = 0; i < elems.size(); ++i) {
          result[power][i] = elems[i];
        }
        current = std::move(current).Mul(base).Rem(factorizing);
      }
    }
    // Convert A to (A - I)^T in place.
    for (size_t i = 0; i < n; ++i) {
      result[i][i] -= Element::One();
      for (size_t j = i + 1; j < n; ++j) {
        std::swap(result[i][j], result[j][i]);
      }
    }
    return result;
  }

  /*! @brief Reduces a square matrix to reduced row echelon form.
   *
   *  Zero rows are removed from the returned matrix.
   */
  inline std::vector<std::vector<Element>> PerformGaussElimination(
      std::vector<std::vector<Element>> matrix) const {
    size_t n = matrix.size();
    size_t row = 0;
    for (size_t column = 0; column < n; ++column) {
      // Find a pivot row with a nonzero value in the current column.
      size_t next_row = row;
      while (next_row < n && matrix[next_row][column] == Element::Zero()) {
        ++next_row;
      }
      if (next_row != n) {
        std::swap(matrix[next_row], matrix[row]);
        // Normalize the pivot row so that the pivot value becomes one.
        Element coefficient = matrix[row][column].Inverse();
        for (size_t i = column; i < n; ++i) {
          matrix[row][i] *= coefficient;
        }
        // Remove this column from every other row. This gives reduced row
        // echelon form, not only an upper triangular form.
        for (size_t other_row = 0; other_row < n; ++other_row) {
          if (row == other_row ||
              matrix[other_row][column] == Element::Zero()) {
            continue;
          }
          coefficient = matrix[other_row][column];
          matrix[other_row][column] = Element::Zero();
          for (size_t i = column + 1; i < n; ++i) {
            matrix[other_row][i] -= matrix[row][i] * coefficient;
          }
        }
        ++row;
      }
    }
    matrix.resize(row);
    return matrix;
  }
};

}  // namespace factorization::solver
