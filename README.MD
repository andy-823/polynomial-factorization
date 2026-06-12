# Polynomial Factorization

This repository contains implementations of polynomial factorization algorithms
over finite fields, together with experiments and benchmarks used to inspect
their behavior.

## Build

The project uses CMake.

```bash
cmake -S . -B build
cmake --build build
```

This builds the library headers, tests, experiments, and benchmarks.

## Experiments

Experiment targets are defined in `experiments/`:

```bash
cmake --build build --target run_experiment_1
cmake --build build --target run_experiment_2
cmake --build build --target run_experiment_3
```

Run them from the `build` directory, because output paths in the experiment
sources are relative to it:

```bash
cd build
./experiments/run_experiment_1
./experiments/run_experiment_2
./experiments/run_experiment_3
```

The generated data is written to:

- `experiments/experiment_1/exp_1_out.txt`
- `experiments/experiment_2/exp_2_out.txt`
- `experiments/experiment_3/exp_3_out.txt`

The experiments use fixed random seeds, so results should be reproducible with
the same build, compiler, and experiment parameters.

## Experiment Data

`experiment_1` compares two versions of Berlekamp factorization and records
operation counts related to polynomial division and Gaussian elimination.

`experiment_2` collects more stable statistics for parts of Berlekamp
factorization using a larger number of runs.

`experiment_3` measures the average number of factors of each degree in
factorizations of random polynomials of a fixed degree.

Experiment parameters are currently set directly in
`experiments/experiment_*/experiment.cpp`: number of runs, degree ranges,
thread count, and selected fields.

## Plotting

Plots can be built from the generated `exp_*_out.txt` files. Plotting scripts
are expected to read the corresponding experiment output and save figures either
next to the data files or into a separate output directory.

If experiment parameters are changed in `experiment.cpp`, rebuild and rerun the
corresponding target before regenerating plots.

## Benchmarks

The distinct-degree factorization benchmark is built with:

```bash
cmake --build build --target benchmark_ddf
```

Run it from `build`:

```bash
cd build
./benchmarks/benchmark_ddf
```

## Tests

Tests are built with the project. They can be run directly from `build`:

```bash
cd build
./tests/test_polynomial
./tests/test_solver
./tests/test_galois_field
./tests/test_galois_field_element
```
