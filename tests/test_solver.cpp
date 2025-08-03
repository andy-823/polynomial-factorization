#include <cstdint>
#include <map>
#include <random>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/simple_polynomial.hpp>
#include <factorization/solver/berlekamp.hpp>

#include "generator.hpp"

using namespace factorization;  // NOLINT

TEST_CASE("Berlekamp") {
  std::mt19937 random_gen;

  SECTION("RandomPowers1") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    constexpr int kBudget = 100;
    constexpr int kTestCount = 1000;

    std::vector<Poly> polynoms = {
      Poly({0, 1}),
      Poly({1, 1 }),
      Poly({1, 1, 1}),  // 1 + x + x^2
      Poly({1, 0, 1, 1}),  // 1 + x^2 + x^3
      Poly({1, 1, 0, 1}),  // 1 + x + x^3
      Poly({1, 1, 0, 0, 1}),  // 1 + x + x^4
      Poly({1, 0, 0, 1, 1}),  // 1 + x^3 + x^4
      Poly({1, 1, 1, 1, 1})
    };

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestCount; ++test) {
      std::map<Poly, int> expected;
      int budget = kBudget;
      Poly factorizing(Element::One());

      static_assert(utils::BinPow(2, 6) == 64);
      for (const auto& poly : polynoms) {
        int power = budget != 0 ? random_gen() % budget : 0;
        budget -= power;
        if (power != 0) {
          expected[poly] = power;
          auto tmp = utils::BinPow(poly, power);
          factorizing *= utils::BinPow(poly, power);
        }
      }

      auto factors = solver.Factorize(factorizing);
      for (const auto& factor : factors) {
        REQUIRE(factor.power == expected[factor.factor]);
      }
    }
  }

  SECTION("RandomPowers2") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    constexpr int kBudget = 100;
    constexpr int kTestCount = 1000;

    std::vector<Poly> polynoms = {
      Poly({0, 1}),
      Poly({1, 1}),
      Poly({2, 1}),
      Poly({3, 1}),
      Poly({4, 1}),
      Poly({5, 1}),
      Poly({6, 1}),
      Poly({7, 1}),
    };

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestCount; ++test) {
      std::map<Poly, int> expected;
      int budget = kBudget;
      Poly factorizing(Element::One());

      static_assert(utils::BinPow(2, 6) == 64);
      for (const auto& poly : polynoms) {
        int power = std::min(uint64_t{10}, random_gen() % (budget + 1));
        budget -= power;
        if (power != 0) {
          expected[poly] = power;
          auto tmp = utils::BinPow(poly, power);
          factorizing *= utils::BinPow(poly, power);
        }
      }

      auto factors = solver.Factorize(factorizing);
      for (const auto& factor : factors) {
        REQUIRE(factor.power == expected[factor.factor]);
      }
    }
  }

  SECTION("Stress") {
    using GaloisField = galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
    // using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    // using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    constexpr int kTestsCount = 1000;

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestsCount; ++test) {
      auto poly = GenPoly<Poly>(random_gen);
      poly.MakeMonic();

      auto check = Poly(Element::One());

      for (const auto& factor : solver.Factorize(poly)) {
        auto got = solver.Factorize(factor.factor);
        std::vector<solver::Factor<Poly>> expected(1, solver::Factor<Poly>(factor.factor, 1));

        REQUIRE(got == expected);
  
        check *= utils::BinPow(factor.factor, factor.power);
      }
      REQUIRE(poly == check);
    }
  }
}