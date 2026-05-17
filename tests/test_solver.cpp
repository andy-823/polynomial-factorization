#include <cstdint>
#include <map>
#include <random>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/polynomial/naive_polynomial.hpp>
#include <factorization/solver/berlekamp.hpp>
#include <factorization/solver/distinct_degree_factorization.hpp>
#include <factorization/solver/square_free_factorization.hpp>

#include "generator.hpp"

#include <iostream>

using namespace factorization;  // NOLINT

template <concepts::Polynom Poly>
Poly BinPow(Poly base, int power) {
  using Element = typename Poly::Element;
  Poly result(Element(1));
  while (power > 0) {
    if (power % 2 != 0) {
      result = std::move(result).Mul(base);
    }
    base = base.Mul(base);
    power /= 2;
  }
  return result;
}

TEST_CASE("Berlekamp") {
  std::mt19937 random_gen;

  SECTION("RandomPowers1") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::NaivePolynomial<Element>;

    constexpr int kBudget = 100;
    constexpr int kTestCount = 1000;

    std::vector<Poly> polynoms = {
        Poly(std::vector<int>{0, 1}),          Poly(std::vector<int>{1, 1}),
        Poly(std::vector<int>{1, 1, 1}),        // 1 + x + x^2
        Poly(std::vector<int>{1, 0, 1, 1}),     // 1 + x^2 + x^3
        Poly(std::vector<int>{1, 1, 0, 1}),     // 1 + x + x^3
        Poly(std::vector<int>{1, 1, 0, 0, 1}),  // 1 + x + x^4
        Poly(std::vector<int>{1, 0, 0, 1, 1}),  // 1 + x^3 + x^4
        Poly(std::vector<int>{1, 1, 1, 1, 1})};

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestCount; ++test) {
      std::map<Poly, int> expected;
      int budget = kBudget;
      Poly factorizing{Element::One()};

      static_assert(utils::BinPow(2, 6) == 64);
      for (const auto& poly : polynoms) {
        int power = budget != 0 ? random_gen() % budget : 0;
        budget -= power;
        if (power != 0) {
          expected[poly] = power;
          auto tmp = BinPow(poly, power);
          factorizing = std::move(factorizing).Mul(BinPow(poly, power));
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
    using Poly = polynomial::NaivePolynomial<Element>;

    constexpr int kBudget = 100;
    constexpr int kTestCount = 1000;

    std::vector<Poly> polynoms = {
        Poly(std::vector<Element>{Element({0, 0, 0}), Element(1)}),
        Poly(std::vector<Element>{Element({1, 0, 0}), Element(1)}),
        Poly(std::vector<Element>{Element({0, 1, 0}), Element(1)}),
        Poly(std::vector<Element>{Element({1, 1, 0}), Element(1)}),
        Poly(std::vector<Element>{Element({0, 0, 1}), Element(1)}),
        Poly(std::vector<Element>{Element({1, 0, 1}), Element(1)}),
        Poly(std::vector<Element>{Element({0, 1, 1}), Element(1)}),
        Poly(std::vector<Element>{Element({1, 1, 1}), Element(1)}),
    };

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestCount; ++test) {
      std::map<Poly, int> expected;
      int budget = kBudget;
      Poly factorizing{Element::One()};

      static_assert(utils::BinPow(2, 6) == 64);
      for (const auto& poly : polynoms) {
        int power =
            std::min(uint64_t{10}, uint64_t(random_gen() % (budget + 1)));
        budget -= power;
        if (power != 0) {
          expected[poly] = power;
          auto tmp = BinPow(poly, power);
          factorizing = std::move(factorizing).Mul(BinPow(poly, power));
        }
      }

      auto factors = solver.Factorize(factorizing);
      REQUIRE(factors.size() == expected.size());
      for (const auto& factor : factors) {
        REQUIRE(factor.power == expected[factor.factor]);
      }
    }
  }

  SECTION("Stress") {
    using GaloisField =
        galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
    // using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    // using GaloisField = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::NaivePolynomial<Element>;

    constexpr int kTestsCount = 1000;

    solver::Berlekamp<Poly> solver;

    for (int test = 0; test < kTestsCount; ++test) {
      auto poly = GenPoly<Poly>(random_gen).MakeMonic();

      auto check = Poly(Element::One());

      for (const auto& factor : solver.Factorize(poly)) {
        auto got = solver.Factorize(factor.factor);
        std::vector<solver::Factor<Poly>> expected(
            1, solver::Factor<Poly>(factor.factor, 1));

        REQUIRE(got == expected);

        check = std::move(check).Mul(BinPow(factor.factor, factor.power));
      }
      REQUIRE(poly == check);
    }
  }
}

TEST_CASE("DistinctDegreeFactorization") {
  std::mt19937 random_gen;

  SECTION("KnownDegrees") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    const Poly degree_1_a(std::vector<int>{0, 1});
    const Poly degree_1_b(std::vector<int>{1, 1});
    const Poly degree_2(std::vector<int>{1, 1, 1});
    const Poly degree_3_a(std::vector<int>{1, 0, 1, 1});
    const Poly degree_3_b(std::vector<int>{1, 1, 0, 1});

    const Poly degree_1_part = degree_1_a.Mul(degree_1_b);
    const Poly degree_3_part = degree_3_a.Mul(degree_3_b);
    const Poly factorizing = degree_1_part.Mul(degree_2).Mul(degree_3_part);

    const auto factors = ddf::naive::DistinctDegreeFactorize(factorizing);

    REQUIRE(factors.size() == 3);
    REQUIRE(factors[0].degree == 1);
    REQUIRE(factors[0].factor == degree_1_part);
    REQUIRE(factors[1].degree == 2);
    REQUIRE(factors[1].factor == degree_2);
    REQUIRE(factors[2].degree == 3);
    REQUIRE(factors[2].factor == degree_3_part);
  }

  SECTION("RandomSquareFreeProducts") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    const std::vector<std::pair<Poly, int>> irreducibles = {
        {Poly(std::vector<int>{0, 1}), 1},
        {Poly(std::vector<int>{1, 1}), 1},
        {Poly(std::vector<int>{1, 1, 1}), 2},
        {Poly(std::vector<int>{1, 0, 1, 1}), 3},
        {Poly(std::vector<int>{1, 1, 0, 1}), 3},
        {Poly(std::vector<int>{1, 1, 0, 0, 1}), 4},
        {Poly(std::vector<int>{1, 0, 0, 1, 1}), 4},
    };

    constexpr int kTestsCount = 100;
    for (int test = 0; test < kTestsCount; ++test) {
      std::map<int, Poly> expected;
      Poly factorizing(Element::One());

      for (const auto& [factor, degree] : irreducibles) {
        if (random_gen() % 2 == 0) {
          continue;
        }
        factorizing = std::move(factorizing).Mul(factor);
        auto it = expected.find(degree);
        if (it == expected.end()) {
          expected.emplace(degree, factor);
        } else {
          it->second = std::move(it->second).Mul(factor);
        }
      }
      if (expected.empty()) {
        continue;
      }
      const auto factors = ddf::naive::DistinctDegreeFactorize(factorizing);

      REQUIRE(factors.size() == expected.size());
      for (const auto& [factor, degree] : factors) {
        REQUIRE(expected.contains(degree));
        REQUIRE(factor == expected[degree]);
      }
    }
  }

  SECTION("StressAgainstBerlekamp") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    using Key = std::pair<int, int>;

    constexpr int kTestsCount = 1000;
    solver::Berlekamp<Poly> berlekamp;

    for (int test = 0; test < kTestsCount; ++test) {
      Poly poly = GenPoly<Poly, 500>(random_gen).MakeMonic();
      if (poly.Size() <= 1) {
        continue;
      }

      std::map<Key, Poly> expected;
      for (const auto& [factor, power] : berlekamp.Factorize(poly)) {
        const Key key{power, static_cast<int>(factor.Size() - 1)};
        auto it = expected.find(key);
        if (it == expected.end()) {
          expected.emplace(key, factor);
        } else {
          it->second = std::move(it->second).Mul(factor).MakeMonic();
        }
      }

      std::map<Key, Poly> actual;
      for (const auto& [square_free_factor, power] :
           sff::SquareFreeFactorize(poly)) {
        for (const auto& [factor, degree] :
             ddf::naive::DistinctDegreeFactorize(square_free_factor)) {
          const Key key{power, degree};
          auto it = actual.find(key);
          if (it == actual.end()) {
            actual.emplace(key, factor);
          } else {
            it->second = std::move(it->second).Mul(factor).MakeMonic();
          }
        }
      }

      REQUIRE(actual == expected);
    }
  }
}
