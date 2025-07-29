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

  SECTION("Stress") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
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