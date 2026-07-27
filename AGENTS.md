# Repository Guidelines

## Project Structure & Module Organization

MSTM is a Fortran 2023 multiple-sphere T-matrix solver. `src/` is the only implementation. `input_state.f90`, `input_parser.f90`, and `input_reporting.f90` form the input pipeline; `simulation_execution.f90` runs one case, while the `*_averaging.f90` modules coordinate averaged cases. `input_execution.f90`, `simulation_averaging.f90`, `input.f90`, and `special_functions.f90` are compatibility façades. Focused numerical modules contain Bessel, angular, translation, scattering, near-field, and GPFA code. `tests/` contains focused numerical tests, while `examples/` contains smoke-test inputs, position data, and plotting assets. The maintained manual is `docs/manual.md`; archived material is in `docs/reference/`.

## Build, Test, and Development Commands

GNU Fortran 15.2+ is required. Build and test serial mode with:

```sh
cmake -S . -B build/serial -DCMAKE_BUILD_TYPE=Release
cmake --build build/serial -j
ctest --test-dir build/serial --output-on-failure
```

Build MPI mode with `-DMSTM_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER=mpifort`. To enable serial/2-/4-rank equivalence tests, also pass `-DMSTM_SERIAL_REFERENCE_EXECUTABLE="$PWD/build/serial/mstm"`. Run a case as `build/serial/mstm examples/mstm-2022b-fig1.inp`; paths written by the solver are relative to the working directory.

For a diagnostic build, add `-DMSTM_ENABLE_RUNTIME_CHECKS=ON`. Warnings are enabled by default; CI also sets `-DMSTM_WARNINGS_AS_ERRORS=ON`.

Generate a line, function, and branch coverage report with `cmake -S . -B build/coverage -DMSTM_ENABLE_COVERAGE=ON`, followed by `cmake --build build/coverage --target coverage`. The report is written to `build/coverage/coverage/coverage.txt` and uploaded by CI.

## Coding Style & Naming Conventions

Use free-form, standard-conforming Fortran 2023 in `.f90` files. Use descriptive lowercase `snake_case` names for project code; reserve single letters for mathematical indices and short loops. Keep canonical GPFA and MINPACK routine names for provenance. Use three-space indentation and `implicit none` in every program unit. Run `fprettify --config-file .fprettify.rc --silent src/*.f90` before review. Clang-format does not support Fortran, and `.clang-format-ignore` protects these sources from its C++ parser. Prefer generic standard intrinsics (`abs`, `conjg`, `cmplx`) and Fortran 2023 degree intrinsics when inputs are degrees. Do not add compiler extensions or restore `-fallow-argument-mismatch`.

Define shared mathematical constants only in `src/constants.f90`; runtime coefficient tables belong in `src/numerical_tables.f90`. Use standard Fortran intrinsics where their argument domains match the algorithm, but retain the complex-argument Bessel implementations because the standard intrinsics accept real arguments only.

## Testing Guidelines

CTest covers analytical Mie identities, Bessel and GPFA kernels, concurrent kernels, translation and solver invariants, repeated state transitions, CLI behavior, averaging modes, nested spheres, plane surfaces, periodic lattices, near fields, and the December 2023 effective-medium case. MPI validation includes direct comparisons of an actual serial executable with two- and four-rank results. Every change must pass all configured serial tests; changes to MPI or shared numerical code must also pass an MPI-enabled build. For numerical changes, compare output against a known-good run and document tolerances, compiler flags, and MPI rank count. Add focused tests when existing cases do not cover the behavior.

## Commit & Pull Request Guidelines

Use a short imperative subject naming the affected behavior or module. Keep commits single-purpose. Pull requests should explain numerical impact, list exact build/test commands, identify output changes, and link related issues. Include plots only when calculated or notebook-generated results change. Never commit `build/`, module/object files, or generated `.dat` output.
