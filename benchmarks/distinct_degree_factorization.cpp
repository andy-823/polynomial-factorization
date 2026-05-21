#include <chrono>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/prime_ring.hpp>
#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/polynomial/karatsuba_engine.hpp>
#include <factorization/polynomial/ntt_engine.hpp>
#include <factorization/solver/common.hpp>
#include <factorization/solver/distinct_degree_factorization.hpp>
#include <factorization/solver/square_free_factorization.hpp>
#include <tests/generator.hpp>

using namespace factorization;  // NOLINT

template <typename Poly, int kTestsCount = 20, size_t kMaxSize = 1000>
void RunDdfBenchmark(const char* field_name) {
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
    const auto start = Timer::now();
    for (const auto& input : inputs) {
      (void)factorize(input);
    }
    const auto finish = Timer::now();
    return std::chrono::duration_cast<Duration>(finish - start).count();
  };

  const auto own_lazy_time = measure([&](const Poly& value) {
    ddf::own_lazy::DistinctDegreeFactorizer<Poly> factorizer(value);
    return factorizer.Run();
  });
  std::map<std::string, long long> tree_profile;
  std::map<std::string, Timer::time_point> tree_profile_start;
  const auto own_tree_time = measure([&](const Poly& value) {
    ddf::own_tree::DistinctDegreeFactorizer<Poly> factorizer(value);
    return factorizer.RunWithObserver([&](const char* stage, bool started) {
      if (started) {
        tree_profile_start[stage] = Timer::now();
        return;
      }
      const auto finish = Timer::now();
      tree_profile[stage] += std::chrono::duration_cast<Duration>(
                                 finish - tree_profile_start[stage])
                                 .count();
    });
  });

  const double ratio = own_tree_time == 0
                           ? 0.0
                           : static_cast<double>(own_lazy_time) /
                                 static_cast<double>(own_tree_time);

  std::cout << "DDF " << field_name << ": inputs=" << inputs.size()
            << ", own_lazy=" << own_lazy_time << "us"
            << ", own_tree=" << own_tree_time << "us"
            << ", own_lazy/own_tree=" << ratio << '\n';
  for (const auto& [stage, time] : tree_profile) {
    const double percent = own_tree_time == 0
                               ? 0.0
                               : 100.0 * static_cast<double>(time) /
                                     static_cast<double>(own_tree_time);
    std::cout << "  own_tree " << stage << ": " << time << "us (" << percent
              << "%)\n";
  }
}

template <typename GaloisField, int kTestsCount = 20, size_t kMaxSize = 1000>
void RunDdfBenchmarkKaratsuba(const char* field_name) {
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Engine = polynomial::KaratsubaEngine<Element>;
  using Poly = polynomial::GenericPolynomial<Element, Engine>;

  RunDdfBenchmark<Poly, kTestsCount, kMaxSize>(field_name);
}

template <typename GaloisField, int kTestsCount = 20, size_t kMaxSize = 1000>
void RunDdfBenchmarkNtt(const char* field_name) {
  using Element = galois_field::FieldElementWrapper<GaloisField>;
  using Engine = polynomial::NttEngine<Element>;
  using Poly = polynomial::GenericPolynomial<Element, Engine>;

  RunDdfBenchmark<Poly, kTestsCount, kMaxSize>(field_name);
}

TEST_CASE("DdfBenchmarkKaratsuba") {
  // NOLINTBEGIN
  using GF2_3 = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
  using GF3_2 = galois_field::LogBasedField<3, 2, {2, 2, 1}>;
  using GF2_16 = galois_field::LogBasedField<
      2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
  // NOLINTEND

  // RunDistinctDegreeFactorizationBenchmark<GF2_3, 5, 3000>("GF(2^3)");
  // RunDistinctDegreeFactorizationBenchmark<GF3_2, 5, 2000>("GF(3^2)");
  RunDdfBenchmarkKaratsuba<GF2_16, 5, 3000>("GF(2^16)");
}

TEST_CASE("DdfBenchmarkNtt") {
  // NOLINTBEGIN
  using Z_17 = galois_field::PrimeRing<100'003>;
  // NOLINTEND
  RunDdfBenchmarkNtt<Z_17, 5, 10000>("GF(1e5+3)");
}