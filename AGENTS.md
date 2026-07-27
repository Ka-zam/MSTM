# Repository Guidelines

## Project Structure & Module Organization

MSTM is a Fortran 2023 multiple-sphere T-matrix solver. `src/` is the only implementation. Files are organized by module: `input_state.f90`, `input_parser.f90`, `input_execution.f90`, and `input_reporting.f90` form the input pipeline; focused numerical modules contain Bessel, angular, translation, scattering, near-field, and GPFA code. `input.f90` and `special_functions.f90` are compatibility façades. `tests/` contains focused numerical tests, while `examples/` contains smoke-test inputs, position data, and plotting assets. The maintained manual is `docs/manual.md`; archived material is in `docs/reference/`.

## Build, Test, and Development Commands

GNU Fortran 15.2 or newer is required. Build and test serial mode with:

```sh
cmake -S . -B build/serial -DCMAKE_BUILD_TYPE=Release
cmake --build build/serial -j
ctest --test-dir build/serial --output-on-failure
```

Build MPI mode with `-DMSTM_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER=mpifort`. Run a case as `build/serial/mstm examples/mstm-2022b-fig1.inp`; paths written by the solver are relative to the working directory.

## Coding Style & Naming Conventions

Use free-form, standard-conforming Fortran 2023 in `.f90` files. Use descriptive lowercase `snake_case` names for modules, procedures, types, and long-lived variables; reserve single-letter names for conventional mathematical indices and short loops. Use three-space indentation and `implicit none` in every program unit. Run `fprettify --config-file .fprettify.rc --silent src/*.f90` before review. Clang-format does not support Fortran, and `.clang-format-ignore` protects these sources from its C++ parser. Prefer generic standard intrinsics (`abs`, `conjg`, `cmplx`) and Fortran 2023 degree intrinsics when inputs are degrees. Do not add compiler extensions or restore `-fallow-argument-mismatch`.

Define shared mathematical constants only in `src/constants.f90`; runtime coefficient tables belong in `src/numerical_tables.f90`. Use standard Fortran intrinsics where their argument domains match the algorithm, but retain the complex-argument Bessel implementations because the standard intrinsics accept real arguments only.

## Testing Guidelines

CTest runs focused Bessel, GPFA, concurrent-kernel, and nested-sphere translation unit tests, the Figure 1 near-field case, and the December 2023 effective-medium case. Every change must pass all six in serial mode; changes to MPI or shared numerical code must also pass an MPI-enabled build. For numerical changes, compare output against a known-good run and document tolerances, compiler flags, and MPI rank count. Add focused tests when existing cases do not cover the behavior.

## Commit & Pull Request Guidelines

Use a short imperative subject naming the affected behavior or module. Keep commits single-purpose. Pull requests should explain numerical impact, list exact build/test commands, identify output changes, and link related issues. Include plots only when calculated or notebook-generated results change. Never commit `build/`, module/object files, or generated `.dat` output.
