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

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <map>
#include <vector>

#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>
#include <factorization/solver/common.hpp>

namespace factorization::solver {

template <concepts::Polynom Polynom>
class BerlekampExperiment {
  using Element = typename Polynom::Element;
  static_assert(Element::kCounting == true, "Only counting element is allowed here");

 public:
  inline std::vector<Factor<Polynom>> Factorize(Polynom polynom) {
    std::vector<Factor<Polynom>> result;
    polynom.MakeMonic();
    if (polynom.IsZero() || polynom.IsOne()) {
      return {};
    }
    std::map<Polynom, int> factors = FactorizeImpl(std::move(polynom));
    result.reserve(factors.size());
    for (const auto& [factor, power] : factors) {
      result.emplace_back(factor, power);
    }
    return result;
  }

  uint64_t GetGaussActions() const {
    return gauss_actions_;
  }

  uint64_t GetDivisionsActions() const {
    return divisions_actions_;
  }

  uint64_t GetTotalActions() const {
    return total_actions_;
  }

 private:
  inline std::map<Polynom, int> FactorizeImpl(Polynom polynom) {
    Element::ResetActions();
    gauss_actions_ = 0;
    divisions_actions_ = 0;
    total_actions_ = 0;

    std::map<Polynom, int> result;
    while (!polynom.IsOne()) {
      Polynom derivative = polynom.Derivative();
      if (derivative.IsZero()) {
        polynom = FieldBaseRoot(polynom);  // this polynom will be monic
        // very unlikely case so using recursion inside seems ok
        for (const auto& [factor, power] : FactorizeImpl(std::move(polynom))) {
          result[factor] += power * Polynom::Element::FieldBase();
        }
        break;  // factorization is complete!
      }
      // gcd is already monic
      Polynom gcd = Gcd(polynom, derivative);
      // polynom / gcd has no repeating factors
      for (const auto& factor : SquareFreeFactorize(polynom / gcd)) {
        ++result[factor];
      }
      polynom = gcd;
    }
    total_actions_ = Element::GetActions();
    return result;
  }

  // assume our polynom is
  //   f(x) = g(x)^p
  // this function returns such g(x) by given f(x)
  inline Polynom FieldBaseRoot(const Polynom& polynom) {
    constexpr int kFieldBase = Element::FieldBase();
    constexpr int kFieldPower = Element::FieldPower();

    std::vector<Element> elements(polynom.GetElements());
    for (size_t i = 0; i < elements.size(); i += kFieldBase) {
      // want to get y that y^p = x
      // consider p is field base, q = p^k is field size
      // then y = x^{q/p} = x^p^{k - 1}
      // Note:
      //   Berlekamp algorithm works only with relatively small fields
      //   caring about overflow is meaningless
      //   if your field is so big use another algo
      constexpr int64_t kPower = utils::BinPow(kFieldBase,
                                               kFieldPower - 1);
      elements[i / kFieldBase] = elements[i].Pow(kPower);
    }
    elements.resize((elements.size() + kFieldBase - 1) / kFieldBase);
    return Polynom(std::move(elements));
  }

  inline Polynom DoModuloPower(Polynom polynom, const auto& matrix) {
    auto elements = polynom.GetElements();
    std::vector<Element> result(matrix.size(), Element::Zero());
    for (size_t i = 0; i < elements.size(); ++i) {
      for (size_t j = 0; j < matrix.size(); ++j) {
        result[j] += elements[i] * matrix[i][j];
      }
    }
    // for (size_t _ = 0; _ < Element::FieldPower(); ++_) {
    //   size_t power = Element::FieldBase();

    //   Polynom squared = std::move(polynom);
    //   polynom = Polynom({Element::One()});
    //   while (power > 0) {
    //     if (power % 2 != 0) {
    //       polynom = polynom * squared % modulo;
    //     }
    //     squared = squared * squared % modulo;
    //     power /= 2;
    //   }
    // }
    return Polynom(result);
  }

  inline std::vector<Polynom> SquareFreeFactorize(Polynom polynom) {
    size_t n = polynom.Size() - 1;
    std::vector<std::vector<Element>> matrix(n, std::vector<Element>(n));
    {
      constexpr int kFieldSize = utils::BinPow(Element::FieldBase(),
                                               Element::FieldPower());
      Polynom base;
      Polynom current(Element::One());
      {
        // here we get that x, done a little weird
        std::vector<Element> tmp(kFieldSize + 1);
        tmp.back() = Element::One();  // the only nonzero element is last
        base = Polynom(std::move(tmp)) % polynom;
      }
      for (size_t power = 0; power < n; ++power) {
        auto elems = current.GetElements();
        for (size_t i = 0; i < elems.size(); ++i) {
          matrix[power][i] = elems[i];
        }
        current = current * base % polynom;
      }
    }
    
    std::vector<Polynom> result;

    Polynom factorizing = polynom;
    Polynom current({Element::Zero(), Element::One()});
    Polynom x = current;
    size_t power = 1;
    while (2 * power < factorizing.Size()) {
      current = DoModuloPower(std::move(current), matrix);
      Polynom gcd = Gcd(factorizing, current - x);
      if (gcd.Size() > 1) {
        factorizing /= gcd;
        if (gcd.Size() == power + 1) {
          result.emplace_back(std::move(gcd));
        } else {
          auto divisors = DistinctDegreeFactorize(gcd, power);
          for (auto& divisor : divisors) {
            result.emplace_back(std::move(divisor));
          }
        }
      }
      ++power;
    }
    if (factorizing.Size() > 1) {
      result.emplace_back(std::move(factorizing));
    }
    return result;
  }

  // input is monic f(x) = f_1(x) ... f_k(x)
  // where f_1 ... f_k are irreducible
  // return vector because of no repeating factors
  inline std::vector<Polynom> DistinctDegreeFactorize(Polynom polynom, size_t power) {
    
    uint64_t before_gauss = Element::GetActions();
    std::vector<Polynom> basis = FindFactorizingBasis(polynom);
    gauss_actions_ = Element::GetActions() - before_gauss;

    // that means polynom is irreducible
    if (basis.size() == 1) {
      return {polynom};
    }
    // this is supposed to be range
    // but inside it can be everything - better to get it here
    const auto field_elements = Element::AllFieldElements();
    std::vector<Polynom> factors = {polynom};
    // put outside of cycle to remove amount of dynamic allocations
    // i understand that most of them are done inside polynoms arithmetic
    // but if we can use less - why not?
    std::vector<Polynom> new_factors;
    new_factors.reserve(basis.size());

    uint64_t before_divisions = Element::GetActions();
    for (const auto& factorizing : basis) {
      if (factorizing.Size() == 1) {
        continue;
      }
      for (const auto& factor : factors) {
        for (const auto& c : field_elements) {
          Polynom new_factor = Gcd(factor, factorizing - c);
          // new factor is non trivial
          if (!new_factor.IsOne()) {
            new_factors.emplace_back(std::move(new_factor));
          }
          // it means that we have already found all necessary factors
          if (new_factors.size() == basis.size()) {
            divisions_actions_ = Element::GetActions() - before_divisions;
            return new_factors;
          }
        }
      }
      // again: less dynamic allocations
      factors.swap(new_factors);
      new_factors.clear();
    }
    divisions_actions_ = Element::GetActions() - before_divisions;
    return factors;
  }

  // find basis of Berlekamp subalgebra
  // it consists of polynomials g which
  //   g^q = g (mod f)
  // q is field size
  inline std::vector<Polynom> FindFactorizingBasis(const Polynom& polynom) {
    std::vector<Polynom> result;
    // since powering to q-th power is linear, it can be done with matrix
    // we want not to power but to find specific polynomials
    // this matrix has view (A - E)^T
    //   where yA = y^q
    auto matrix = BuildMatrix(polynom);
    // now we want to find basis of solutions Ax = 0
    // at first we want to perform gauss elimination
    matrix = PerformGaussElimination(std::move(matrix));
    // Idea i use here is the following
    // We have positions of free coefficients
    // Example:
    //   0 1 1 0
    //   0 0 0 1
    // Positions of free coefficients are [0, 3]
    // For each row we remember which column is "main"
    // Here we have main columns:
    //   row 0: column 1
    //   row 1: column 3
    // Then i apply each free coefficient to one, other zero
    // And retrieve necessary vector via matrix
    // For example above:
    //   first free coefficient is one: we get [1, 0, 0, 0]
    //   second free coefficient is one: we get [0, -1, 1, 0]
    // As free coefficient is at only one column we this equations
    // We have equations containing only 2 values
    // Example for second
    //   x_0 = 0, here c_0 = 0
    //   x_1 + 1 * x_2 = 0
    //   x_2 = 1, here c_1 = 1
    //   x_3 + 0 * x_2 = 0
    size_t rank = matrix.size();
    size_t n = polynom.Size() - 1;
    std::vector<size_t> free_coefficient_positions;
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
      // free coefficient after the last data position
      while (column < n) {
        free_coefficient_positions.emplace_back(column);
        ++column;
      }
    }

    for (const auto& column : free_coefficient_positions) {
      std::vector<Element> current(n, Element::Zero());
      current[column] = Element::One();
      for (size_t row = 0; row < rank; ++row) {
        current[data_position[row]] = -matrix[row][column];
      }
      result.emplace_back(std::move(current));
    }
    return result;
  }

  // returns (A - E)^T
  // where A is equivalent to powering to q-th power
  inline std::vector<std::vector<Element>> BuildMatrix(const Polynom& factorizing) {
    constexpr int kFieldSize = utils::BinPow(Element::FieldBase(),
                                             Element::FieldPower());
    size_t n = factorizing.Size() - 1;
    std::vector<std::vector<Element>> result(n, std::vector<Element>(n));
    // At start we want to build matrix A, such that
    //   y A = y^q (mod f)
    // where y is polynom
    // so we want to fill its rows this way
    //   A_0     = x^{0 * q} (mod f)
    //   A_1     = x^{1 * q} (mod f)
    //   ...
    //   A_{n-1} = x^{(n - 1) * q} (mod f)
    // here we have q is field size
    {
      // base is equal to x^q % f
      Polynom base;
      Polynom current(Element::One());
      {
        // here we get that x, done a little weird
        std::vector<Element> tmp(kFieldSize + 1);
        tmp.back() = Element::One();  // the only nonzero element is last
        base = Polynom(std::move(tmp)) % factorizing;
      }
      for (size_t power = 0; power < n; ++power) {
        auto elems = current.GetElements();
        for (size_t i = 0; i < elems.size(); ++i) {
          result[power][i] = elems[i];
        }
        current = current * base % factorizing;
      }
    }
    // yA = y
    // y(A - E) = 0
    // (A - E)^T y^T = 0
    // so we want to get (A - E)^T
    for (size_t i = 0; i < n; ++i) {
      result[i][i] -= Element::One();
      for (size_t j = i + 1; j < n; ++j) {
        std::swap(result[i][j], result[j][i]);
      }
    }
    return result;
  }

  inline std::vector<std::vector<Element>> PerformGaussElimination(
    std::vector<std::vector<Element>> matrix) {

    size_t n = matrix.size();
    size_t row = 0;
    for (size_t column = 0; column < n; ++column) {
      // at first we want to find row which differs from zero at column
      size_t next_row = row;
      while (next_row < n && matrix[next_row][column] == Element::Zero()) {
        ++next_row;
      }
      if (next_row != n) {
        // we have found this row
        std::swap(matrix[next_row], matrix[row]);
        // we want to set coefficient at matrix[row][column] to one
        Element coefficient = matrix[row][column].Inverse();
        for (size_t i = column; i < n; ++i) {
          matrix[row][i] *= coefficient;
        }
        for (size_t other_row = 0; other_row < n; ++other_row) {
          if (row == other_row ||
              matrix[other_row][column] == Element::Zero()) {
            continue;
          }
          // we want to perform A_i = A_i - A_j * a
          // where a = matrix[other_row][column]
          coefficient = matrix[other_row][column];
          matrix[other_row][column] = Element::Zero();
          for (size_t i = column + 1; i < n; ++i) {
            matrix[other_row][i] -= matrix[row][i] * coefficient;
          }
        }
        // now we can go to other row
        ++row;
      }
    }
    matrix.resize(row);
    return matrix;
  }

 private:
  uint64_t gauss_actions_ = 0;
  uint64_t divisions_actions_ = 0;
  uint64_t total_actions_ = 0;
};

}  // namespace factorization::solver