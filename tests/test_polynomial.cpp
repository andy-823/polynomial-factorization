#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/polynomial/naive_polynomial.hpp>

#include "generator.hpp"

using namespace factorization;         // NOLINT
using namespace std::chrono_literals;  // NOLINT

template <concepts::Polynom Poly, concepts::GaloisFieldElement Elem>
Poly MakePoly(std::vector<int> data) {
  std::vector<Elem> result;
  result.reserve(data.size());
  for (const auto& v : data) {
    result.emplace_back(Elem(v));
  }
  return Poly(std::move(result));
}

template <concepts::Polynom Reference, concepts::Polynom Verify,
          size_t kMaxSize, typename RandomGen>
void RunCompareTest(RandomGen& random_gen) {
  auto first = GenPoly<Reference, kMaxSize>(random_gen);
  auto second = GenPoly<Reference, kMaxSize>(random_gen);

  Verify generic_first(first.Get());
  Verify generic_second(second.Get());

  REQUIRE(generic_first.Mul(generic_second).Get() == first.Mul(second).Get());
  REQUIRE(generic_first.Div(generic_second).Get() == first.Div(second).Get());
  REQUIRE(generic_first.Rem(generic_second).Get() == first.Rem(second).Get());

  const auto second_modulus = generic_second.BuildModulus(generic_first.Size());
  REQUIRE(generic_first.Div(second_modulus).Get() == first.Div(second).Get());
  REQUIRE(generic_first.Rem(second_modulus).Get() == first.Rem(second).Get());

  const auto [generic_div, generic_rem] = generic_first.DivRem(generic_second);
  const auto [naive_div, naive_rem] = first.DivRem(second);
  REQUIRE(generic_div.Get() == naive_div.Get());
  REQUIRE(generic_rem.Get() == naive_rem.Get());

  const auto [modulus_div, modulus_rem] = generic_first.DivRem(second_modulus);
  REQUIRE(modulus_div.Get() == naive_div.Get());
  REQUIRE(modulus_rem.Get() == naive_rem.Get());
  REQUIRE(generic_first.Gcd(generic_second).MakeMonic().Get() ==
          first.Gcd(second).Get());
}

TEST_CASE("NaivePolynomial") {
  std::mt19937 random_gen;

  SECTION("Add sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 1, {1, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::NaivePolynomial<Element>;

    {
      Poly first = MakePoly<Poly, Element>({1, 0, 1, 0, 0, 1, 1});
      Poly second = first;

      REQUIRE(first.Add(second).IsZero());
      REQUIRE(first.Sub(second).IsZero());
    }

    {
      Poly first(MakePoly<Poly, Element>({1, 0, 1, 0, 1, 1}));
      Poly expected(MakePoly<Poly, Element>({0, 0, 1, 0, 1, 1}));

      REQUIRE(first.Add(Element(1)) == expected);
      REQUIRE(first.Add(Element(1)) == expected);
    }

    {
      Poly first(MakePoly<Poly, Element>({1, 0, 1, 0, 1, 1}));
      Poly second(MakePoly<Poly, Element>({1, 0, 1, 0, 0, 1}));
      Poly expected(MakePoly<Poly, Element>({0, 0, 0, 0, 1}));

      REQUIRE(first.Add(second) == expected);
    }
  }

  SECTION("Multiply sanity check") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::NaivePolynomial<Element>;

    {
      Poly poly(MakePoly<Poly, Element>({1, 0, 1, 0, 1, 1}));

      REQUIRE(poly.Div(poly) == MakePoly<Poly, Element>({1}));
      REQUIRE(poly.Div(poly).IsOne());

      REQUIRE(poly.Rem(MakePoly<Poly, Element>({1})).IsZero());
      REQUIRE(poly.Mul(MakePoly<Poly, Element>({0})).IsZero());
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
    using Poly = polynomial::NaivePolynomial<Element>;

    {
      Poly poly(MakePoly<Poly, Element>({1, 0, 1, 0, 1, 1}));

      REQUIRE(poly.Derivative() == MakePoly<Poly, Element>({0, 0, 0, 0, 1}));
    }

    {
      using Vec = std::vector<Element>;
      Poly poly(Vec(4, Element({1, 1, 0})));
      REQUIRE(poly.MakeMonic() == Poly(Vec(4, Element::One())));
    }

    {
      using Vec = std::vector<Element>;
      Poly first(Vec{Element({1, 0, 0}), Element({0, 1, 0}), Element({1, 1, 0}),
                     Element({0, 0, 1}), Element({1, 0, 1}), Element({0, 1, 1}),
                     Element({1, 1, 1})});
      Poly second(Vec{Element({0, 0, 0}), Element({1, 0, 0}),
                      Element({0, 1, 0}), Element({1, 1, 0}),
                      Element({0, 0, 1}), Element({1, 0, 1}),
                      Element({0, 1, 1})});

      REQUIRE(second < first);

      first =
          Poly(Vec{Element({1, 0, 0}), Element({0, 1, 0}), Element({1, 1, 0}),
                   Element({0, 0, 1}), Element({1, 0, 1}), Element({0, 1, 1}),
                   Element({1, 1, 1})});
      second =
          Poly(Vec{Element({1, 0, 0}), Element({1, 0, 0}), Element({1, 1, 0}),
                   Element({1, 1, 0}), Element({0, 0, 1}), Element({0, 1, 1}),
                   Element({1, 1, 1})});

      REQUIRE(second < first);

      first =
          Poly(Vec{Element({1, 0, 0}), Element({0, 1, 0}), Element({0, 1, 1})});
      second = Poly(Vec{Element({1, 0, 0}), Element({0, 1, 0})});

      REQUIRE(second < first);
    }
  }

  SECTION("Stress") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Poly = polynomial::NaivePolynomial<Element>;

    constexpr int kTestsCount = 10000;

    for (int test = 0; test < kTestsCount; ++test) {
      auto first = GenPoly<Poly, 128>(random_gen);
      auto second = GenPoly<Poly, 128>(random_gen);

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

TEST_CASE("GenericPolynomial") {
  std::mt19937 random_gen;

  SECTION("Small") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    constexpr int kTestsCount = 10000;

    for (int test = 0; test < kTestsCount; ++test) {
      RunCompareTest<NaivePoly, GenericPoly, 32>(random_gen);
    }
  }

  SECTION("Average") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    constexpr int kTestsCount = 10000;

    for (int test = 0; test < kTestsCount; ++test) {
      RunCompareTest<NaivePoly, GenericPoly, 256>(random_gen);
    }
  }

  SECTION("Big") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    constexpr int kTestsCount = 500;

    for (int test = 0; test < kTestsCount; ++test) {
      RunCompareTest<NaivePoly, GenericPoly, 1000>(random_gen);
    }
  }

  SECTION("Mul speed") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;
    using Timer = std::chrono::steady_clock;
    using Duration = std::chrono::milliseconds;

    auto get_duration_count = [](auto time) {
      return std::chrono::duration_cast<Duration>(time).count();
    };

    auto first_naive = GenPoly<NaivePoly, 100'000, kFixed>(random_gen);
    auto second_naive = GenPoly<NaivePoly, 100'000, kFixed>(random_gen);

    GenericPoly first_generic(first_naive.Get());
    GenericPoly second_generic(second_naive.Get());

    auto start_naive = Timer::now();
    auto res_naive = first_naive.Mul(second_naive);
    auto end_naive = Timer::now();

    auto res_generic = first_generic.Mul(second_generic);
    auto modulus = first_generic.BuildModulus(res_generic.Size());
    auto start_generic = Timer::now();
    REQUIRE(res_generic.Rem(modulus).IsZero());
    auto end_generic = Timer::now();

    std::cout << "Naive: " << get_duration_count(end_naive - start_naive)
              << " ms\n";
    std::cout << "Generic: " << get_duration_count(end_generic - start_generic)
              << " ms\n";
    REQUIRE(res_naive.Get() == res_generic.Get());
    REQUIRE(end_naive - start_naive > 5 * (end_generic - start_generic));
  }

  SECTION("GCD speed") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;
    using Timer = std::chrono::steady_clock;
    using Duration = std::chrono::milliseconds;

    auto get_duration_count = [](auto time) {
      return std::chrono::duration_cast<Duration>(time).count();
    };

    auto common = GenPoly<NaivePoly, 2'000, kFixed>(random_gen);
    auto first_tail = GenPoly<NaivePoly, 20'000, kFixed>(random_gen);
    auto second_tail = GenPoly<NaivePoly, 20'000, kFixed>(random_gen);
    auto first_naive = common.Mul(first_tail);
    auto second_naive = common.Mul(second_tail);

    GenericPoly first_generic(first_naive.Get());
    GenericPoly second_generic(second_naive.Get());

    auto start_naive = Timer::now();
    auto res_naive = first_naive.Gcd(second_naive);
    auto end_naive = Timer::now();

    auto start_generic = Timer::now();
    auto res_generic = first_generic.Gcd(second_generic).MakeMonic();
    auto end_generic = Timer::now();

    std::cout << "Naive GCD speed: " << get_duration_count(end_naive - start_naive)
              << " ms\n";
    std::cout << "Generic GCD speed: "
              << get_duration_count(end_generic - start_generic) << " ms\n";
    REQUIRE(res_generic.Get() == res_naive.Get());
  }
}

template <concepts::Polynom Poly, size_t kMaxPolySize, size_t kMaxModSize,
          typename RandomGen>
void RunCompModFrobeniusTest(RandomGen& random_gen) {
  using Element = typename Poly::Element;

  constexpr size_t kFieldSize =
      utils::BinPow(Element::FieldBase(), Element::FieldPower());
  constexpr size_t kTestsCount = 1000;

  for (size_t test = 0; test < kTestsCount; ++test) {
    Poly mod = GenPoly<Poly, kMaxModSize, kFixed>(random_gen).MakeMonic();

    const auto modulus = mod.BuildModulus(2 * mod.Size() + kMaxPolySize);
    const Poly x(std::vector<Element>{Element::Zero(), Element::One()});
    const Poly x_q = polynomial::BinPowMod(x, kFieldSize, modulus);

    const Poly poly = GenPoly<Poly, kMaxPolySize>(random_gen);
    const size_t t = static_cast<size_t>(std::ceil(std::sqrt(poly.Size())));

    const auto matrix = polynomial::BuildCompModMatrix(x_q, t, modulus);
    const Poly expected =
        polynomial::BinPowMod(poly.Rem(modulus), kFieldSize, modulus);
    const Poly actual = polynomial::CompMod(poly, matrix, modulus);

    REQUIRE(actual == expected);
  }
}

TEST_CASE("CompMod") {
  std::mt19937 random_gen;

  SECTION("GF_2^3") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }

  SECTION("GF_2^8") {
    using GaloisField = galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }

  SECTION("GF_2^16") {
    using GaloisField = galois_field::LogBasedField<2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::PolynomialEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }
}
