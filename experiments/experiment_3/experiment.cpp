#include <atomic>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <mutex>
#include <random>
#include <string>

#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/naive_polynomial.hpp>
#include <factorization/runtime/thread_pool.hpp>
#include <factorization/runtime/wait_group.hpp>
#include <factorization/solver/berlekamp.hpp>
#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

#include "berlekamp.hpp"
#include "generator.hpp"

using namespace factorization;  // NOLINT

struct ExperimentParams {
  size_t test_value;
  size_t thread_count;
  size_t test_runs;
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

template <concepts::GaloisField Field, typename RandomGen>
void RunExperiment(std::ostream& out, const ExperimentParams& params,
                   RandomGen& gen) {
  using Element = galois_field::FieldElementWrapper<Field>;
  using Poly = polynomial::NaivePolynomial<Element>;

  constexpr int kFieldSize =
      utils::BinPow(Field::FieldBase(), Field::FieldPower());

  runtime::ThreadPool thread_pool(params.thread_count);
  runtime::WaitGroup wg;
  thread_pool.Start();

  std::vector<std::atomic<int64_t>> counts(params.test_value + 1);

  out << kFieldSize << "\t\t";
  wg.Add(params.test_runs);
  for (size_t test = 0; test < params.test_runs; ++test) {
    runtime::SubmitTask(&thread_pool, [&] {
      solver::BerlekampExperiment<Poly> solver;

      RandomGen local_gen;
      local_gen.seed(gen());

      Poly poly = GenPoly<Poly>(local_gen, params.test_value);
      for (const auto& [factor, power] : solver.Factorize(poly)) {
        counts[factor.Size() - 1].fetch_add(1);
      }
      wg.Done();
    });
  }
  wg.Wait();
  for (const auto& count : counts) {
    double avg = static_cast<double>(count.load());
    out << std::setprecision(3) << std::fixed << avg / params.test_runs << "\t";
  }
  out << "\n";

  thread_pool.Stop();
}

int main() {
  // https://www.partow.net/programming/polynomials/index.html
  // NOLINTBEGIN
  using GF2_1 = galois_field::LogBasedField<2, 1, {1, 1}>;
  using GF2_2 = galois_field::LogBasedField<2, 2, {1, 1, 1}>;
  using GF2_3 = galois_field::LogBasedField<2, 3, {1, 1, 0, 1}>;
  using GF2_4 = galois_field::LogBasedField<2, 4, {1, 1, 0, 0, 1}>;
  using GF2_5 = galois_field::LogBasedField<2, 5, {1, 0, 1, 0, 0, 1}>;
  using GF2_6 = galois_field::LogBasedField<2, 6, {1, 1, 0, 0, 0, 0, 1}>;
  using GF2_7 = galois_field::LogBasedField<2, 7, {1, 1, 0, 0, 0, 0, 0, 1}>;
  using GF2_8 = galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;
  using GF2_9 = galois_field::LogBasedField<2, 9, {1, 0, 0, 0, 1, 0, 0, 0, 0, 1}>;
  using GF2_10 = galois_field::LogBasedField<2, 10, {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1}>;
  using GF2_16 = galois_field::LogBasedField<2, 16, {1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1}>;
  // NOLINTEND

  constexpr int kRuns = 10'000;
  constexpr int kTestValue = 1000;
  constexpr int kThreads = 12;  // 20;
  const std::string path{"../../experiments/experiment_3/exp_3_out.txt"};

  ExperimentParams params;
  params.test_value = kTestValue;
  params.thread_count = kThreads;
  params.test_runs = kRuns;

  std::ofstream out_file;
  out_file.open(path);
  if (!out_file.is_open()) {
    std::cerr << "failed to open output file: " << path << "\n";
    return 0;
  }
  std::ostream& out = out_file;

  out << "\t\t\t";
  for (size_t i = 0; i <= kTestValue; ++i) {
    out << i << "\t\t";
  }
  out << "\n";

  MultithreadRandomGen<std::mt19937> gen;
  gen.seed(0);

  // out << "\n";
  gen.seed(0);

  RunExperiment<GF2_1>(out, params, gen);
  RunExperiment<GF2_2>(out, params, gen);
  // RunExperiment<GF2_3>(out, params, gen);
  RunExperiment<GF2_4>(out, params, gen);
  // RunExperiment<GF2_5>(out, params, gen);
  // RunExperiment<GF2_6>(out, params, gen);
  // RunExperiment<GF2_7>(out, params, gen);
  RunExperiment<GF2_8>(out, params, gen);
  // RunExperiment<GF2_9>(out, params, gen);
  RunExperiment<GF2_16>(out, params, gen);

  out_file.close();

  return 0;
}
