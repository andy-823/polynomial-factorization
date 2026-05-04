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

      REQUIRE(first.Add(second).IsZero());
      REQUIRE(first.Sub(second).IsZero());
    }

    {
      Poly first({1, 0, 1, 0, 1, 1});
      Poly expected({0, 0, 1, 0, 1, 1});

      REQUIRE(first.Add(Element(1)) == expected);
      REQUIRE(first.Add(Element(1)) == expected);
    }

    {
      Poly first({1, 0, 1, 0, 1, 1});
      Poly second({1, 0, 1, 0, 0, 1});
      Poly expected({0, 0, 0, 0, 1});

      REQUIRE(first.Add(second) == expected);
    }
  }

  SECTION("Multiply sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::SimplePolynomial<Element>;

    {
      Poly poly({1, 0, 1, 0, 1, 1});

      REQUIRE(poly.Div(poly) == Poly{1});
      REQUIRE(poly.Div(poly).IsOne());

      REQUIRE(poly.Rem(Poly{1}).IsZero());
      REQUIRE(poly.Mul(Poly{0}).IsZero());
    }

    {
      using Vec = std::vector<Element>;
      Poly poly(Vec(4, Element({1, 1, 0})));

      REQUIRE(poly.Mul(Element({0, 1, 0})) == Poly(Vec(4, Element({0, 1, 1}))));
      REQUIRE(poly.Div(Element({1, 1, 0})) == Poly(Vec(4, Element::One())));
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
      using Vec = std::vector<Element>;
      Poly poly(Vec(4, Element({1, 1, 0})));
      REQUIRE(poly.MakeMonic() == Poly(Vec(4, Element::One())));
    }

    {
      using Vec = std::vector<Element>;
      Poly first({Element({1, 0, 0}), Element({0, 1, 0}),
                  Element({1, 1, 0}), Element({0, 0, 1}),
                  Element({1, 0, 1}), Element({0, 1, 1}),
                  Element({1, 1, 1})});
      Poly second({Element({0, 0, 0}), Element({1, 0, 0}),
                   Element({0, 1, 0}), Element({1, 1, 0}),
                   Element({0, 0, 1}), Element({1, 0, 1}),
                   Element({0, 1, 1})});

      REQUIRE(second < first );

      first = Poly({Element({1, 0, 0}), Element({0, 1, 0}),
                       Element({1, 1, 0}), Element({0, 0, 1}),
                       Element({1, 0, 1}), Element({0, 1, 1}),
                       Element({1, 1, 1})});
      second = Poly({Element({1, 0, 0}), Element({1, 0, 0}),
                     Element({1, 1, 0}), Element({1, 1, 0}),
                     Element({0, 0, 1}), Element({0, 1, 1}),
                     Element({1, 1, 1})});

      REQUIRE(second < first );

      first = Poly({Element({1, 0, 0}), Element({0, 1, 0}),
                    Element({0, 1, 1})});
      second = Poly({Element({1, 0, 0}), Element({0, 1, 0})});

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
        auto rem = first.Rem(second);
        auto div = first.Div(second);

        REQUIRE(div.Mul(second).Add(rem) == first);
      }

      {
        auto sub = first.Sub(second);

        REQUIRE(sub.Add(second) == first);
      }
    }
  }
}