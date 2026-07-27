# Repository Guidelines

## Project Structure & Module Organization

MSTM is a Fortran 2023 multiple-sphere T-matrix solver. `src/` is the only implementation: `main.f90` drives input and execution, `input.f90` owns configuration/reporting, and the remaining files contain numerical, scattering, FFT, random-configuration, and serial/MPI modules. `examples/` contains paper inputs, position data, and plotting assets. The maintained manual is `docs/manual.md`; archived papers and the former PDF manual are in `docs/reference/`. CMake helpers live in `cmake/`.

## Build, Test, and Development Commands

GNU Fortran 15.2 or newer is required. Build and test serial mode with:

```sh
cmake -S . -B build/serial -DCMAKE_BUILD_TYPE=Release
cmake --build build/serial -j
ctest --test-dir build/serial --output-on-failure
```

Build MPI mode with `-DMSTM_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER=mpifort`. Run a case as `build/serial/mstm examples/mstm-2022b-fig1.inp`; paths written by the solver are relative to the working directory.

## Coding Style & Naming Conventions

Use free-form, standard-conforming Fortran 2023 in `.f90` files. Match the existing lowercase `snake_case` identifiers, three-space indentation, and `implicit none` in every program unit. Run `fprettify --config-file .fprettify.rc --silent src/*.f90` before review. Clang-format does not support Fortran, and `.clang-format-ignore` protects these sources from its C++ parser. Prefer generic standard intrinsics (`abs`, `conjg`, `cmplx`) and standard I/O specifiers. Do not add compiler extensions or restore `-fallow-argument-mismatch`.

## Testing Guidelines

CTest runs the Figure 1 near-field case and the focused December 2023 effective-medium case. Every change must pass both in serial mode; changes to MPI or shared numerical code must also pass an MPI-enabled build. For numerical changes, compare relevant output against a known-good run and document tolerances, compiler flags, and MPI rank count. Add focused example inputs when existing cases do not cover the behavior.

## Commit & Pull Request Guidelines

Use a short imperative subject naming the affected behavior or module. Keep commits single-purpose. Pull requests should explain numerical impact, list exact build/test commands, identify output changes, and link related issues. Include plots only when calculated or notebook-generated results change. Never commit `build/`, module/object files, or generated `.dat` output.
