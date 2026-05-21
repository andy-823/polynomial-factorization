#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include <catch2/catch_test_macros.hpp>

#include <factorization/concepts.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/prime_ring.hpp>
#include <factorization/polynomial/common.hpp>
#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/polynomial/karatsuba_engine.hpp>
#include <factorization/polynomial/naive_polynomial.hpp>
#include <factorization/polynomial/ntt_engine.hpp>

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

template <typename Duration, typename Func>
typename Duration::rep Measure(Func&& func) {
  using Timer = std::chrono::steady_clock;

  const auto start = Timer::now();
  func();
  const auto finish = Timer::now();
  return std::chrono::duration_cast<Duration>(finish - start).count();
}

template <concepts::Polynom First, concepts::Polynom Second, int kGcdDeg,
          int kPolyDeg, typename RandomGen>
void CompareGcdSpeed(RandomGen& random_gen, const char* label1,
                     const char* label2) {
  using Duration = std::chrono::milliseconds;

  auto common = GenPoly<First, kGcdDeg, kFixed>(random_gen);
  auto tail1 = GenPoly<First, kPolyDeg, kFixed>(random_gen);
  auto tail2 = GenPoly<First, kPolyDeg, kFixed>(random_gen);
  auto first1 = common.Mul(tail1);
  auto second1 = common.Mul(tail2);

  Second first2(first1.Get());
  Second second2(second1.Get());

  First res1;
  const auto time1 = Measure<Duration>([&] {
    res1 = first1.Gcd(second1).MakeMonic();
  });
  Second res2;
  const auto time2 = Measure<Duration>([&] {
    res2 = first2.Gcd(second2).MakeMonic();
  });

  std::cout << label1 << ": " << time1 << " ms\n";
  std::cout << label2 << ": " << time2 << " ms\n";

  REQUIRE(res1.Get() == res2.Get());
}

template <concepts::Polynom First, concepts::Polynom Second, size_t kSize,
          typename RandomGen>
void CompareMulSpeed(RandomGen& random_gen, const char* label1,
                     const char* label2) {
  using Duration = std::chrono::microseconds;

  auto first1 = GenPoly<First, kSize, kFixed>(random_gen);
  auto second1 = GenPoly<First, kSize, kFixed>(random_gen);
  Second first2(first1.Get());
  Second second2(second1.Get());

  First result1;
  const auto time1 = Measure<Duration>([&] {
    result1 = first1.Mul(second1);
  });

  Second result2;
  const auto time2 = Measure<Duration>([&] {
    result2 = first2.Mul(second2);
  });

  std::cout << label1 << ": " << time1 << " us\n";
  std::cout << label2 << ": " << time2 << " us\n";
  REQUIRE(result2.Get() == result1.Get());
}

template <concepts::Polynom Poly, size_t kSize, typename RandomGen>
void RunMulSpeed(RandomGen& random_gen, const char* label) {
  auto first = GenPoly<Poly, kSize, kFixed>(random_gen);
  auto second = GenPoly<Poly, kSize, kFixed>(random_gen);

  const auto time = Measure<std::chrono::microseconds>([&] {
    (void)first.Mul(second);
  });

  std::cout << label << ": " << time << " us\n";
}

template <concepts::Polynom Poly, size_t kValueSize, size_t kDivSize,
          typename RandomGen>
void RunRemSpeed(RandomGen& random_gen, const char* label) {
  using Duration = std::chrono::microseconds;

  const auto value = GenPoly<Poly, kValueSize, kFixed>(random_gen);
  const auto divisor = GenPoly<Poly, kDivSize, kFixed>(random_gen).MakeMonic();
  const auto modulus = divisor.BuildModulus(value.Size());

  Poly plain_rem;
  const auto plain_time = Measure<Duration>([&] {
    plain_rem = value.Rem(divisor);
  });

  Poly precomputed_rem;
  const auto precomputed_time = Measure<Duration>([&] {
    precomputed_rem = value.Rem(modulus);
  });

  std::cout << label << "\n";
  std::cout << "Rem without precompute: " << plain_time << " us\n";
  std::cout << "Rem with precompute: " << precomputed_time << " us\n";
  REQUIRE(precomputed_rem == plain_rem);
}

TEST_CASE("GenericPolynomial") {
  std::mt19937 random_gen;

  SECTION("Small") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::KaratsubaEngine<Element>;
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
    using Engine = polynomial::KaratsubaEngine<Element>;
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
    using Engine = polynomial::KaratsubaEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    constexpr int kTestsCount = 500;

    for (int test = 0; test < kTestsCount; ++test) {
      RunCompareTest<NaivePoly, GenericPoly, 1000>(random_gen);
    }
  }

  SECTION("NTT") {
    using GaloisField = galois_field::PrimeRing<17>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::NttEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    constexpr int kTestsCount = 20;

    for (int test = 0; test < kTestsCount; ++test) {
      RunCompareTest<NaivePoly, GenericPoly, 700>(random_gen);
    }
  }

  SECTION("Karatsuba mul speed") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using KaratsubaPoly = polynomial::GenericPolynomial<Element, Engine>;

    CompareMulSpeed<NaivePoly, KaratsubaPoly, 100'000>(
        random_gen, "Naive mul", "NTT mul");
  }

  SECTION("NTT mul speed") {
    using GaloisField = galois_field::PrimeRing<17>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using KaraEngine = polynomial::KaratsubaEngine<Element>;
    using NttEngine = polynomial::NttEngine<Element>;
    using KaratsubaPoly = polynomial::GenericPolynomial<Element, KaraEngine>;
    using NttPoly = polynomial::GenericPolynomial<Element, NttEngine>;

    CompareMulSpeed<KaratsubaPoly, NttPoly, 1'000'000>(
        random_gen, "Karatsuba degree 1m", "NTT degree 1m");
  }

  SECTION("NTT raw speed") {
    using GaloisField = galois_field::PrimeRing<17>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NttEngine = polynomial::NttEngine<Element>;
    using NttPoly = polynomial::GenericPolynomial<Element, NttEngine>;

    RunMulSpeed<NttPoly, 16'384>(random_gen, "NTT ntt_size 32768");
  }

  SECTION("NTT rem speed") {
    using GaloisField = galois_field::PrimeRing<17>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NttPoly =
        polynomial::GenericPolynomial<Element, polynomial::NttEngine<Element>>;

    RunRemSpeed<NttPoly, 1'000'000, 500'000>(random_gen, "NTT Rem 1m");
  }

  SECTION("Karatuba GCD speed") {
    using GaloisField = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using GenericPoly = polynomial::GenericPolynomial<Element, Engine>;

    CompareGcdSpeed<NaivePoly, GenericPoly, 2'000, 20'000>(
        random_gen, "Naive GCD", "Karatsuba GCD");
  }

  SECTION("NTT GCD speed") {
    using GaloisField = galois_field::PrimeRing<17>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::NttEngine<Element>;
    using NaivePoly = polynomial::NaivePolynomial<Element>;
    using NttPoly = polynomial::GenericPolynomial<Element, Engine>;

    CompareGcdSpeed<NaivePoly, NttPoly, 2'000, 20'000>(random_gen, "Naive GCD",
                                                       "NTT GCD");
  }
}

template <concepts::Polynom Poly, size_t kMaxPolySize, size_t kMaxModSize,
          typename RandomGen>
void RunCompModFrobeniusTest(RandomGen& random_gen) {
  using Element = typename Poly::Element;
  constexpr auto kFieldBase = Element::FieldBase();
  constexpr auto kFieldPower = Element::FieldPower();
  constexpr auto kZero = Element::Zero();
  constexpr auto kOne = Element::One();

  constexpr size_t kFieldSize = utils::BinPow(kFieldBase, kFieldPower);
  constexpr size_t kTestsCount = 100;

  for (size_t test = 0; test < kTestsCount; ++test) {
    Poly mod = GenPoly<Poly, kMaxModSize, kFixed>(random_gen).MakeMonic();

    const auto modulus = mod.BuildModulus(2 * mod.Size() + kMaxPolySize);
    const Poly x(std::vector<Element>{kZero, kOne});
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
    using Engine = polynomial::KaratsubaEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }

  SECTION("GF_2^8") {
    using GaloisField =
        galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }

  SECTION("GF_2^16") {
    using GaloisField = galois_field::LogBasedField<
        2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GaloisField>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    RunCompModFrobeniusTest<Poly, 16, 16>(random_gen);
    RunCompModFrobeniusTest<Poly, 256, 128>(random_gen);
    RunCompModFrobeniusTest<Poly, 1024, 1024>(random_gen);
  }
}
