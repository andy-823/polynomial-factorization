#include <atomic>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <random>
#include <vector>

#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/prime_ring.hpp>

#include <factorization/polynomial/generic_polynomial.hpp>
#include <factorization/polynomial/karatsuba_engine.hpp>
#include <factorization/polynomial/ntt_engine.hpp>

#include <factorization/runtime/thread_pool.hpp>
#include <factorization/runtime/wait_group.hpp>

#include <factorization/solver/common.hpp>
#include <factorization/solver/distinct_degree_factorization.hpp>
#include <factorization/solver/square_free_factorization.hpp>

#include "factorization/concepts.hpp"
#include "generator.hpp"

using namespace factorization;  // NOLINT
using namespace std::chrono_literals;  // NOLINT
using Clock = std::chrono::steady_clock;
using Duration = std::chrono::microseconds;

constexpr static int kThreadCount = 1;

struct SimParams {
  std::vector<int> points;
  int run_count;
};

template <typename RandomGen>
class MultithreadRandomGen {
 public:
  using result_type = typename RandomGen::result_type;  // NOLINT

  void seed(result_type seed) {  // NOLINT
    lock_.lock();
    gen_.seed(seed);
    lock_.unlock();
  }

  result_type operator()() {
    lock_.lock();
    result_type result = gen_();
    lock_.unlock();
    return result;
  }

 private:
  RandomGen gen_;
  std::mutex lock_;
};

// assume RandomGen support multithreading
template <concepts::Polynom Poly, typename Solver, typename RandomGen>
int64_t RunPoint(int size, int run_count, RandomGen& random_gen) {
  std::atomic<int64_t> result = 0;

  runtime::ThreadPool rt(kThreadCount);

  rt.Start();
  for (int i = 0; i < run_count; ++i) {
    runtime::SubmitTask(&rt, [&] {
      Poly poly = GenPoly<Poly>(random_gen, size);
      Duration time = 0s;
      for (const auto& [poly, _] : sff::SquareFreeFactorize(poly)) {
        Solver solver(poly);
        auto start = Clock::now();
        (void)solver.Run();
        auto finish = Clock::now();
        time += std::chrono::duration_cast<Duration>(finish - start);
      }
      result.fetch_add(time.count());
    });
  }
  rt.Stop();

  return result;
}

template <concepts::Polynom Poly, typename Solver, typename RandomGen>
void RunSingleSolver(const char* solver_label, std::ostream& out,
                     const SimParams& params, RandomGen& random_gen) {
  out << solver_label << "\t";
  for (const auto& size : params.points) {
    auto total = RunPoint<Poly, Solver>(size, params.run_count, random_gen);
    double average = static_cast<double>(total) / params.run_count;
    out << std::setprecision(3) << std::fixed << average << "\t";
  }
  out << "\n";
}

template <concepts::Polynom Poly, typename RandomGen = std::mt19937_64>
void RunMultipleSolvers(std::ostream& out, const SimParams& params, const uint64_t seed = 0) {
  using ExactNtl = ddf::ntl_like::DistinctDegreeFactorizer<Poly, ddf::kExactNtl>;
  using ModifiedNtl = ddf::ntl_like::DistinctDegreeFactorizer<Poly, ddf::kSmallField>;
  using Lazy = ddf::own_lazy::DistinctDegreeFactorizer<Poly, ddf::kSmallField>;
  using Tree = ddf::own_tree::DistinctDegreeFactorizer<Poly, ddf::kSmallField>;

  MultithreadRandomGen<RandomGen> random_gen;

  random_gen.seed(seed);
  RunSingleSolver<Poly, ExactNtl>("exact_ntl", out, params,random_gen);

  random_gen.seed(seed);
  RunSingleSolver<Poly, ModifiedNtl>("modified_ntl", out, params, random_gen);

  random_gen.seed(seed);
  RunSingleSolver<Poly, Lazy>("lazy", out, params, random_gen);

  random_gen.seed(seed);
  RunSingleSolver<Poly, Tree>("tree", out, params, random_gen);
}

template <concepts::Polynom Poly, typename RandomGen = std::mt19937_64>
void Simulate(const char* label, std::ostream& out,
              const SimParams& params, const uint64_t seed = 0) {
  // print first line
  out << label << "\t";
  for (const auto& size : params.points) {
    out << size << "\t\t";
  }
  out << "\n";
  RunMultipleSolvers<Poly, RandomGen>(out, params, seed);
  out << "\n";
}

int main() {
  SimParams params;
  params.run_count = 1;
  params.points = {1000};

  std::ostream& out = std::cout;

  {
    // NOLINTNEXTLINE
    using GF2_8 = galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GF2_8>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    Simulate<Poly>("GF2^8", out, params);
  }

  {
    // NOLINTNEXTLINE
    using GF2_16 = galois_field::LogBasedField<2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
    using Element = galois_field::FieldElementWrapper<GF2_16>;
    using Engine = polynomial::KaratsubaEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    Simulate<Poly>("GF2^16", out, params);
  }

  {
    // NOLINTNEXTLINE
    using Z_p = galois_field::PrimeRing<100'003>;;
    using Element = galois_field::FieldElementWrapper<Z_p>;
    using Engine = polynomial::NttEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    Simulate<Poly>("Z_100'003", out, params);
  }

  {
    // NOLINTNEXTLINE
    using Z_p = galois_field::PrimeRing<2524775926340780033, uint64_t, __int128_t>;
    using Element = galois_field::FieldElementWrapper<Z_p>;
    using Engine = polynomial::NttEngine<Element>;
    using Poly = polynomial::GenericPolynomial<Element, Engine>;

    Simulate<Poly>("Z_2524775926340780033", out, params);
  }

  return 0;
}