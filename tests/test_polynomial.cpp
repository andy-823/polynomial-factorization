#include <random>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/simple_polynomial.hpp>

#include "generator.hpp"

using namespace factorization;  // NOLINT

TEST_CASE("SimplePolynomial") {
  std::mt19937 random_gen;

  SECTION("Add sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    {
      Poly first({1, 0, 1, 0, 1, 1});
      Poly second = first;

      REQUIRE((first + second).IsZero());
      REQUIRE((first - second).IsZero());
    }

    {
      Poly first({1, 0, 1, 0, 1, 1});
      Poly expected({0, 0, 1, 0, 1, 1});

      REQUIRE(first - Element(1) == expected);
      REQUIRE(first + Element(1) == expected);
    }

    {
      Poly first({1, 0, 1, 0, 1, 1});
      Poly second({1, 0, 1, 0, 0, 1});
      Poly expected({0, 0, 0, 0, 1});

      REQUIRE(first + second == expected);
    }
  }

  SECTION("Multiply sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    {
      Poly poly({1, 0, 1, 0, 1, 1});

      REQUIRE(poly / poly == Poly{1});
      REQUIRE((poly / poly).IsOne());

      REQUIRE((poly % Poly{1}).IsZero());
      REQUIRE((poly * Poly{0}).IsZero());
    }

    {
      Poly poly({3, 3, 3, 3});

      REQUIRE(poly * Element{2} == Poly({6, 6, 6, 6}));
      REQUIRE(poly / Element{3} == Poly({1, 1, 1, 1}));
    }
  }

  SECTION("Other methods sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    {
      Poly poly({1, 0, 1, 0, 1, 1});

      REQUIRE(poly.Derivative() == Poly({0, 0, 0, 0, 1}));
    }

    {
      Poly poly({3, 3, 3, 3});
      poly.MakeMonic();

      REQUIRE(poly == Poly({1, 1, 1, 1}));
    }

    {
      Poly first({1, 2, 3, 4, 5, 6, 7});
      Poly second({0, 1, 2, 3, 4, 5, 6});

      REQUIRE(second < first );

      first = Poly({1, 2, 3, 4, 5, 6, 7});
      second = Poly({1, 1, 3, 3, 4, 6, 7});

      REQUIRE(second < first );

      first = Poly({1, 2, 3});
      second = Poly({1, 2});

      REQUIRE(second < first);
    }
  }

  SECTION("Stress") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    constexpr int kTestsCount = 10000;

    for (int test = 0; test < kTestsCount; ++test) {
      auto first = GenPoly<Poly>(random_gen);
      auto second = GenPoly<Poly>(random_gen);

      {
        auto rem = first % second;
        auto div = first / second;

        REQUIRE(div * second + rem == first);
      }

      {
        auto sub = first - second;

        REQUIRE(sub + second == first);
      }
    }
  }
}