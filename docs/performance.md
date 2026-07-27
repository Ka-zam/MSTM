# Performance Notes

These measurements identify algorithmic bottlenecks; they are not universal benchmarks. They were collected with GNU Fortran 15.2, an optimized serial build, OpenBLAS 0.3.32, and `OPENBLAS_NUM_THREADS=1`.

| Spheres | Unknowns | Pairwise BiCG | FFT BiCG | Dense direct | Direct memory |
|---:|---:|---:|---:|---:|---:|
| 100 | 1,600 | 0.116 s | 0.051 s | 0.334 s | 91 MB |
| 200 | 3,200 | 0.318 s | 0.067 s | 1.83 s | 334 MB |

The 200-sphere direct solve spent approximately 0.16 s assembling the interaction matrix, 1.50 s in LAPACK LU factorization, 0.10 s estimating its condition, and 0.02 s in backsolves. Dense LU is therefore already the limiting direct-solver operation. Its cubic time and quadratic storage make it a validation and small-system method; matrix-free iteration is the scalable path.

For this 3,200-unknown matrix, increasing OpenBLAS from one to two threads reduced factorization only from 1.50 s to 1.43 s while increasing total solve time; four and eight threads were slower still. Do not assume the BLAS default thread count is optimal—benchmark `OPENBLAS_NUM_THREADS=1` for small and medium direct problems.

The 200-sphere FFT run performed 12 interaction applications and 216 batched 3-D transforms. About 0.034 s of a 0.084 s sampled run was inside the transforms. The in-tree GPFA remains a significant target, but not the whole accelerated algorithm: node interpolation, exact neighbor interactions, and translation-matrix initialization also matter.

FFT translation is approximate. At automatic node order 2, the sampled total unpolarized efficiency differed from pairwise/direct translation by 0.04%, with a maximum polarized difference of about 0.24%. Node order 4 reduced the maximum difference to about 0.06%, but cost roughly the same as pairwise translation for this case. Always perform a node-order or pairwise convergence check for production calculations.

An FFTW3 microbenchmark executing the same 832 independent `30 x 30 x 4`
complex transforms took 0.017--0.019 s, excluding plan construction. This is about
twice as fast as the measured GPFA transform time, but would improve only the
transform portion of the complete FFT translation.

An optional FFTW3 backend is therefore reasonable future work. The preferred
design is a batched, planned backend retaining GPFA as the dependency-free
fallback. Replacing GPFA outright is not justified until FFTW3 is benchmarked
across larger cell grids and MPI configurations.
