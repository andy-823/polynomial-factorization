#include <iomanip>
#include <iostream>
#include <fstream>
#include <mutex>
#include <random>

#include <factorization/galois_field/log_based_field.hpp>
#include <factorization/galois_field/field_element_wrapper.hpp>
#include <factorization/polynomial/simple_polynomial.hpp>
#include <factorization/parallel/thread_pool.hpp>
#include <factorization/parallel/wait_group.hpp>
#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

#include "berlekamp.hpp"
#include "generator.hpp"

using namespace factorization;  // NOLINT

struct ExperimentParams {
  size_t min_value;
  size_t max_value;
  size_t step;
  size_t thread_count;
  size_t test_runs;
};

template <typename RandomGen>
class MultithreadRandomGen {
 public:
  using result_type = typename RandomGen::result_type;

  void seed(result_type seed) {
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
void RunExperiment(std::ostream& out, const ExperimentParams& params, RandomGen& gen) {
  using Element = galois_field::FieldElementWrapper<Field>;
  using Poly = polynomial::SimplePolynomial<Element>;

  constexpr int kFieldSize = utils::BinPow(Field::FieldBase(), Field::FieldPower());

  parallel::ThreadPool runtime(params.thread_count);
  parallel::WaitGroup wg;
  runtime.Start();

  out << kFieldSize << "\t";
  for (size_t size = params.min_value; size < params.max_value; size += params.step) {
    wg.Add(params.test_runs);
    solver::BerlekampExperiment<Poly> solver;
    for (size_t test = 0; test < params.test_runs; ++test) {
      parallel::SubmitTask(&runtime, [&]{
        RandomGen local_gen;
        local_gen.seed(gen());

        Poly poly = GenPoly<Poly>(local_gen, size);
        solver.Factorize(poly);
        wg.Done();
      });
    }
    wg.Wait();
    double result = static_cast<double>(solver.GetMetricValue());
    result /= params.test_runs * size;
    out << std::setprecision(2) << std::fixed << result << "\t";
  }
  out << "\n";

  runtime.Stop();
}


int main() {
  // https://www.partow.net/programming/polynomials/index.html
  using GF2_1 = galois_field::LogBasedField<2, 1, {1, 1}>;
  using GF2_2 = galois_field::LogBasedField<2, 2, {1, 1, 1}>;
  using GF2_3 = galois_field::LogBasedField<2, 3, {1, 1, 0 , 1}>;
  using GF2_4 = galois_field::LogBasedField<2, 4, {1, 1, 0, 0, 1}>;
  using GF2_5 = galois_field::LogBasedField<2, 5, {1, 0, 1, 0, 0, 1}>;
  using GF2_6 = galois_field::LogBasedField<2, 6, {1, 1, 0, 0, 0, 0, 1}>;
  using GF2_7 = galois_field::LogBasedField<2, 7, {1, 1, 0, 0, 0, 0, 0, 1}>;
  using GF2_8 = galois_field::LogBasedField<2, 8, {1, 0, 1, 1, 1, 0, 0, 0, 1}>;

  constexpr int kRuns = 10000;
  constexpr int kMin = 10;
  constexpr int kMax = 100;
  constexpr int kStep = 2;
  constexpr int kThreads = 20;

  ExperimentParams params;
  params.min_value = kMin;
  params.max_value = kMax;
  params.step = kStep;
  params.thread_count = kThreads;
  params.test_runs = kRuns;

  std::ofstream out_file;
  out_file.open("/home/udoo/polynomial-factorization/exp_1_out.txt");
  if (!out_file.is_open()) {
    std::cout << "bad\n";
    return 0;
  }
  std::ostream& out = out_file; //std::cout;

  out << "\t";
  for (size_t i = kMin; i < kMax; i += kStep) {
    out << i << "\t";
  }
  out << "\n";

  MultithreadRandomGen<std::mt19937> gen;
  gen.seed(0);

  RunExperiment<GF2_1>(out, params, gen);
  RunExperiment<GF2_2>(out, params, gen);
  RunExperiment<GF2_3>(out, params, gen);
  RunExperiment<GF2_4>(out, params, gen);
  RunExperiment<GF2_5>(out, params, gen);
  RunExperiment<GF2_6>(out, params, gen);
  RunExperiment<GF2_7>(out, params, gen);
  RunExperiment<GF2_8>(out, params, gen);

  out_file.close();

  return 0;
}