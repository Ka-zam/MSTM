# MSTM

MSTM is a Fortran multiple-sphere T-matrix solver for monochromatic electromagnetic scattering. It supports fixed and random orientations, nested spheres, plane boundaries, periodic lattices, near-field maps, effective-medium calculations, and serial or MPI execution.

The canonical implementation is in `src/`. It incorporates the features from the December 2023 split release and targets the Fortran 2023 dialect supported by GNU Fortran 15.2.

## Build

Configure and build the serial executable:

```sh
cmake -S . -B build/serial -DCMAKE_BUILD_TYPE=Release
cmake --build build/serial -j
```

For MPI, enable the MPI backend and select an MPI GNU Fortran wrapper if CMake does not find it automatically:

```sh
cmake -S . -B build/mpi -DCMAKE_BUILD_TYPE=Release \
  -DMSTM_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER=mpifort
cmake --build build/mpi -j
```

## Run and test

Pass an input file on the command line. With no argument, MSTM reads `mstm.inp` from the working directory.

```sh
build/serial/mstm examples/mstm-2022b-fig1.inp
ctest --test-dir build/serial --output-on-failure
mpiexec -n 4 build/mpi/mstm examples/mstm-2022b-fig1.inp
```

See [the user manual](docs/manual.md) for the input format, calculation modes, and option reference. Reproduction inputs and position data are in `examples/`; published papers and the archived PDF manual are in `docs/reference/`.

## Formatting

Install the development formatter and format all Fortran sources with:

```sh
python3 -m pip install -r requirements-dev.txt
fprettify --config-file .fprettify.rc --silent src/*.f90
```

The conventional `.f90` suffix identifies free-form source; CMake selects the Fortran 2023 language level independently. Clang-format does not parse Fortran, so `.clang-format` applies only to any supported-language files added later and `.clang-format-ignore` prevents accidental changes to `src/*.f90`.
