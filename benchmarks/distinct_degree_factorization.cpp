#include <chrono>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/solver/common.hpp>
#include <factorization/solver/distinct_degree_factorization.hpp>
#include <factorization/solver/square_free_factorization.hpp>
#include <tests/generator.hpp>

using namespace factorization;  // NOLINT

template <typename Poly>
std::map<int, Poly> Normalize(
    const std::vector<ddf::DistinctDegreeFactor<Poly>>& factors) {
  std::map<int, Poly> result;
  for (const auto& [factor, degree] : factors) {
    result.emplace(degree, factor);
  }
  return result;
}

template <typename GaloisField, int kTestsCount = 20, size_t kMaxSize = 1000>
void RunDistinctDegreeFactorizationBenchmark(const char* field_name) {
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Engine = polynomial::PolynomialEngine<Element>;
  using Poly = polynomial::GenericPolynomial<Element, Engine>;
  using Timer = std::chrono::steady_clock;
  using Duration = std::chrono::microseconds;

  std::mt19937 random_gen;

  std::vector<Poly> inputs;
  for (int test = 0; test < kTestsCount; ++test) {
    Poly poly = GenPoly<Poly, kMaxSize, kFixed>(random_gen).MakeMonic();
    if (poly.Size() <= 1) {
      continue;
    }
    for (const auto& [factor, _] : sff::SquareFreeFactorize(poly)) {
      if (factor.Size() > 1) {
        inputs.emplace_back(factor);
      }
    }
  }

  REQUIRE(!inputs.empty());

  auto measure = [&](auto factorize) {
    std::vector<std::map<int, Poly>> result;
    result.reserve(inputs.size());
    const auto start = Timer::now();
    for (const auto& input : inputs) {
      result.emplace_back(Normalize(factorize(input)));
    }
    const auto finish = Timer::now();
    return std::pair{
        std::move(result),
        std::chrono::duration_cast<Duration>(finish - start).count()};
  };

  auto [naive_result, naive_time] = measure([](const Poly& value) {
    return ddf::naive::DistinctDegreeFactorize(value);
  });

  std::map<std::string, long long> profile;
  std::map<std::string, Timer::time_point> profile_start;
  auto [ntl_like_result, ntl_like_time] = measure([&](const Poly& value) {
    ddf::ntl_like::DistinctDegreeFactorizer<Poly> factorizer(value);
    return factorizer.RunWithObserver([&](const char* stage, bool started) {
      if (started) {
        profile_start[stage] = Timer::now();
        return;
      }
      const auto finish = Timer::now();
      profile[stage] +=
          std::chrono::duration_cast<Duration>(finish - profile_start[stage])
              .count();
    });
  });

  REQUIRE(ntl_like_result == naive_result);

  const double ratio = ntl_like_time == 0
                           ? 0.0
                           : static_cast<double>(naive_time) / ntl_like_time;
  std::cout << "DDF " << field_name << ": inputs=" << inputs.size()
            << ", naive=" << naive_time << "us"
            << ", ntl_like=" << ntl_like_time << "us"
            << ", naive/ntl_like=" << ratio << '\n';
  for (const auto& [stage, time] : profile) {
    const double percent = ntl_like_time == 0
                               ? 0.0
                               : 100.0 * static_cast<double>(time) /
                                     static_cast<double>(ntl_like_time);
    std::cout << "  " << stage << ": " << time << "us (" << percent << "%)\n";
  }
}

TEST_CASE("DistinctDegreeFactorizationBenchmark") {
  // NOLINTBEGIN
  using GF2_3 = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
  using GF3_2 = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using GF2_16 = galois_field::LogBasedField<
      2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
  // NOLINTEND

  // RunDistinctDegreeFactorizationBenchmark<GF2_3, 5, 3000>("GF(2^3)");
  // RunDistinctDegreeFactorizationBenchmark<GF3_2, 5, 2000>("GF(3^2)");
  RunDistinctDegreeFactorizationBenchmark<GF2_16, 1, 5000>("GF(2^16)");
}
