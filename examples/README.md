# Examples

The five `mstm-2022b-fig*.inp` files reproduce calculations from the 2022 paper in `docs/reference/`. Run them from a disposable working directory because output paths are relative to the process working directory.

```sh
mkdir -p build/example
cd build/example
../serial/mstm ../../examples/mstm-2022b-fig1.inp
```

`frac_agg.pos` and `random_in_cylinder_1000.pos` provide sphere positions for aggregate and FFT-acceleration examples. `plots mstm-2022b.nb` is the Mathematica notebook used to plot the paper results.

`effective-medium-smoke.inp` is a deliberately small regression case for the effective-medium fitting workflow added in December 2023. CTest runs it alongside the Figure 1 smoke test.
