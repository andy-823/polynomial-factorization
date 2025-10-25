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
#include <factorization/polynomial/simple_polynomial.hpp>
#include <factorization/parallel/thread_pool.hpp>
#include <factorization/parallel/wait_group.hpp>
#include <factorization/solver/berlekamp.hpp>
#include <factorization/concepts.hpp>
#include <factorization/utils.hpp>

#include "berlekamp.hpp"
#include "counting_field_element.hpp"
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
void RunDivisionsExperiment(std::ostream& out, const ExperimentParams& params, RandomGen& gen) {
  using Element = galois_field::CountingFieldElement<Field>;
  using Poly = polynomial::SimplePolynomial<Element>;

  constexpr int kFieldSize = utils::BinPow(Field::FieldBase(), Field::FieldPower());

  parallel::ThreadPool runtime(params.thread_count);
  parallel::WaitGroup wg;
  runtime.Start();

  out << kFieldSize << "\t";
  for (size_t size = params.min_value; size <= params.max_value; size += params.step) {
    std::atomic<uint64_t> actions_simple{0};
    std::atomic<uint64_t> actions_better{0};
    std::atomic<uint64_t> gauss{0};


    wg.Add(params.test_runs);
    for (size_t test = 0; test < params.test_runs; ++test) {
      parallel::SubmitTask(&runtime, [&]{    
        solver::BerlekampExperiment<Poly, false> solver_simple;
        solver::BerlekampExperiment<Poly, true> solver_better;

        RandomGen local_gen;
        local_gen.seed(gen());

        Poly poly = GenPoly<Poly>(local_gen, size);
        auto factors_simple = solver_simple.Factorize(poly);
        auto factors_better = solver_better.Factorize(poly);
        assert(factors_simple == factors_better);

        actions_simple.fetch_add(solver_simple.GetDivisionsActions());
        actions_better.fetch_add(solver_better.GetDivisionsActions());

        wg.Done();
      });
    }
    wg.Wait();
    double avg_simple = static_cast<double>(actions_simple.load());
    double avg_better = static_cast<double>(actions_better.load());
    avg_simple /= params.test_runs;
    avg_better /= params.test_runs;

    out << std::setprecision(1) << std::fixed << avg_simple << " ";
    out << std::setprecision(1) << std::fixed << avg_better << "\t";
  }
  out << "\n";

  runtime.Stop();
}

template <concepts::GaloisField Field, typename RandomGen>
void RunGaussExperiment(std::ostream& out, const ExperimentParams& params, RandomGen& gen) {
  using Element = galois_field::CountingFieldElement<Field>;
  using Poly = polynomial::SimplePolynomial<Element>;

  constexpr int kFieldSize = utils::BinPow(Field::FieldBase(), Field::FieldPower());

  parallel::ThreadPool runtime(params.thread_count);
  parallel::WaitGroup wg;
  runtime.Start();

  out << kFieldSize << "\t";
  for (size_t size = params.min_value; size <= params.max_value; size += params.step) {
    std::atomic<uint64_t> divisions{0};
    std::atomic<uint64_t> gauss{0};
    std::atomic<uint64_t> total{0};

    wg.Add(params.test_runs);
    for (size_t test = 0; test < params.test_runs; ++test) {
      parallel::SubmitTask(&runtime, [&]{    
        solver::BerlekampExperiment<Poly> solver;

        RandomGen local_gen;
        local_gen.seed(gen());

        Poly poly = GenPoly<Poly>(local_gen, size);
        auto factors = solver.Factorize(poly);

        divisions.fetch_add(solver.GetDivisionsActions());
        gauss.fetch_add(solver.GetGaussActions());
        total.fetch_add(solver.GetTotalActions());

        // auto check = solver::Berlekamp<Poly>().Factorize(poly);
        // auto cmp = [](auto first, auto second) {
        //   if (first.power != second.power) {
        //     return first.power < second.power;
        //   }
        //   return first.factor < second.factor;
        // };
        // std::sort(factors.begin(), factors.end(), cmp);
        // std::sort(check.begin(), check.end(), cmp);
        // assert(check == factors);

        wg.Done();
      });
    }
    wg.Wait();
    double avg_divisions = static_cast<double>(divisions.load());
    double avg_gauss = static_cast<double>(gauss.load());
    double avg_total = static_cast<double>(total.load());
    avg_divisions /= params.test_runs;
    avg_gauss /= params.test_runs;
    avg_total /= params.test_runs;

    out << std::setprecision(1) << std::fixed << avg_divisions << " ";
    out << std::setprecision(1) << std::fixed << avg_gauss << "\t";
    out << std::setprecision(1) << std::fixed << avg_total << "\t";
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
  using GF2_9 = galois_field::LogBasedField<2, 9, {1, 0, 0, 0, 1, 0, 0, 0, 0, 1}>;


  constexpr int kRuns = 100; //10'000;
  constexpr int kMin = 25;
  constexpr int kMax = 300; //400; //500;
  constexpr int kStep = 25;
  constexpr int kThreads = 12; //20;
  const std::string path{"../../experiments/experiment_2/exp_2_out.txt"};

  ExperimentParams params;
  params.min_value = kMin;
  params.max_value = kMax;
  params.step = kStep;
  params.thread_count = kThreads;
  params.test_runs = kRuns;

  std::ofstream out_file;
  out_file.open(path);
  if (!out_file.is_open()) {
    std::cout << "bad\n";
    return 0;
  }
  std::ostream& out = out_file; //std::cout;

  out << "\t";
  for (size_t i = kMin; i <= kMax; i += kStep) {
    out << i << "\t\t\t";
  }
  out << "\n";

  MultithreadRandomGen<std::mt19937> gen;
  gen.seed(0);

  // RunDivisionsExperiment<GF2_1>(out, params, gen);
  // RunDivisionsExperiment<GF2_2>(out, params, gen);
  // RunDivisionsExperiment<GF2_3>(out, params, gen);
  // RunDivisionsExperiment<GF2_4>(out, params, gen);
  // RunDivisionsExperiment<GF2_5>(out, params, gen);
  // RunDivisionsExperiment<GF2_6>(out, params, gen);
  // RunDivisionsExperiment<GF2_7>(out, params, gen);
  // RunDivisionsExperiment<GF2_8>(out, params, gen);
  // RunDivisionsExperiment<GF2_9>(out, params, gen);

  out << "\n";
  gen.seed(0);

  RunGaussExperiment<GF2_1>(out, params, gen);
  RunGaussExperiment<GF2_2>(out, params, gen);
  RunGaussExperiment<GF2_3>(out, params, gen);
  RunGaussExperiment<GF2_4>(out, params, gen);
  RunGaussExperiment<GF2_5>(out, params, gen);
  RunGaussExperiment<GF2_6>(out, params, gen);
  RunGaussExperiment<GF2_7>(out, params, gen);
  RunGaussExperiment<GF2_8>(out, params, gen);
  RunGaussExperiment<GF2_9>(out, params, gen);

  out_file.close();

  return 0;
}