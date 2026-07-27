# MSTM User Manual

Multiple-sphere T-matrix electromagnetic scattering solver

Version 4.0 · Fortran 2023 edition · July 2026

This Markdown manual replaces the Version 4.0 PDF dated 10 January 2023. The physical model and input conventions are retained, while build instructions and source organization describe the consolidated December 2023 implementation. Questions about the scientific model may be directed to the original author at `mstm@auburn.edu`.

## 1. Purpose and capabilities

MSTM calculates monochromatic electromagnetic scattering by multiple spherical domains using the multiple-sphere T-matrix method. Results are exact to the truncation error of the vector spherical wave-function (VSWF) expansions. MSTM can:

- calculate cross sections, asymmetry parameters, and far-field scattering matrices for fixed or random orientations;
- model non-overlapping exterior, nested, and optically active spheres;
- use plane-wave or Gaussian-beam illumination;
- include parallel plane boundaries separating media with different refractive indices;
- repeat a unit cell as a two-dimensional periodic lattice;
- calculate electric and magnetic near fields on regular grids;
- accelerate large, non-periodic systems with an FFT-based translation algorithm; and
- generate and average random configurations, including effective-medium estimates.

The serial backend and MPI distributed-memory backend implement the same numerical path.

## 2. Physical conventions

All values are dimensionless. Length is scaled by the free-space wave number $k_0=2\pi/\lambda_0$. Electric-field amplitudes are relative to the incident wave, and the time dependence is $\exp(-i\omega t)$.

Sphere surface $i$ has radius $a_i$, position $\mathbf r_i=(x_i,y_i,z_i)$, and a pair of complex refractive indices $(m_{iL},m_{iR})$. Equal left and right indices describe an isotropic material. Sphere surfaces may be nested but may not intersect.

$N_B$ parallel boundaries normal to $z$ divide the embedding medium into $N_B+1$ layers. The first boundary is fixed at $z=0$; subsequent positions are formed from `layer_thickness`. Exterior spheres must remain wholly inside one layer. A unit cell can repeat in $x$ and $y$ with periods $W_x,W_y$; repeated spheres must not overlap.

The incident propagation direction uses polar angle $\beta$ and azimuth $\alpha$. Random-orientation calculations require no plane boundaries and no periodicity.

## 3. Build and execution

### 3.1 Requirements

- CMake 3.24 or newer
- GNU Fortran 15.2 or newer
- an MPI implementation with a GNU Fortran wrapper for parallel builds

The build enforces the published Fortran 2023 dialect (`-std=f2023`).

### 3.2 Serial build

```sh
cmake -S . -B build/serial -DCMAKE_BUILD_TYPE=Release
cmake --build build/serial -j
ctest --test-dir build/serial --output-on-failure
```

### 3.3 MPI build

```sh
cmake -S . -B build/mpi -DCMAKE_BUILD_TYPE=Release \
  -DMSTM_ENABLE_MPI=ON -DCMAKE_Fortran_COMPILER=mpifort
cmake --build build/mpi -j
ctest --test-dir build/mpi --output-on-failure
```

### 3.4 Running

```sh
build/serial/mstm path/to/case.inp
mpiexec -n 8 build/mpi/mstm path/to/case.inp
build/serial/mstm --help
build/serial/mstm --version
```

`help`, `--help`, and `-h` print command-line usage. `version`, `--version`, and `-V` print the version generated from the CMake project metadata. With no argument, MSTM reads `mstm.inp` from the working directory. Input paths may be absolute or relative to that directory. Output paths in the input file are also resolved from the working directory. Multiples of four MPI ranks are often most efficient for configuration averaging.

## 4. Input format

Input normally consists of two-line option/value pairs. Option names are case-sensitive and must appear exactly as documented. The last occurrence wins when an option is repeated. Omitted options retain defaults.

Values may be integers, real values (`1.0` or `1.d0`), complex pairs (`(1.6d0,0.01d0)`), logicals (`t`, `f`, `.true.`, `.false.`), or unquoted strings.

```text
number_spheres
100
sphere_data_input_file
examples/frac_agg.pos
output_file
aggregate.dat
length_scale_factor
0.5d0
ref_index_scale_factor
(1.6d0,0.1d0)
random_orientation
t
end_of_options
```

`end_of_options` must terminate the final run.

### 4.1 Sphere data

`sphere_data_input_file` supplies at least `number_spheres` records. A record has 3–6 fields separated by whitespace, tabs, or commas:

```text
x  y  z  [radius]  [(isotropic_index)]  [(right_index)]
```

- Three fields use `length_scale_factor` as the common radius.
- Four fields use `ref_index_scale_factor` as the common refractive index.
- Five fields add an isotropic complex index.
- Six fields give left and right indices for optically active material.

Positions and explicit radii are multiplied by `length_scale_factor`. Indices are multiplied by `ref_index_scale_factor`.

Data can instead be embedded in the input:

```text
number_spheres
2
sphere_data
0.d0,0.d0,0.d0,1.d0,(1.5d0,0.d0)
0.d0,0.d0,2.d0,1.d0,(1.5d0,0.d0)
end_of_sphere_data
```

## 5. Common calculation recipes

### 5.1 Random orientation

The analytical T-matrix average is enabled with `random_orientation=t`:

```text
output_file
random.dat
number_spheres
100
sphere_data_input_file
examples/frac_agg.pos
length_scale_factor
0.5d0
ref_index_scale_factor
(1.6d0,0.1d0)
random_orientation
t
calculate_scattering_matrix
t
end_of_options
```

For a Monte Carlo average, set `random_orientation=f`, `incidence_average=t`, and choose `number_incident_directions`. Increase the sample count until reported properties converge.

### 5.2 Fixed orientation and nested spheres

```text
output_file
fixed.dat
number_spheres
2
sphere_data
0.d0,0.d0,0.d0,20.d0,(1.54d0,0.d0)
0.d0,0.d0,19.d0,1.d0,(0.15d0,3.d0)
end_of_sphere_data
incident_beta_deg
30.d0
incident_alpha_deg
0.d0
incident_frame
t
scattering_map_model
0
end_of_options
```

With `incident_frame=t`, zero scattering angle is the incident direction. With it false, angles use the target frame and $\theta=0$ is the positive $z$ axis.

### 5.3 Plane boundaries

The first boundary is always $z=0$. This example places a unit-radius air bubble in a glass layer of thickness 2:

```text
number_spheres
1
sphere_data
0.d0,0.d0,1.d0,1.d0,(1.d0,0.d0)
end_of_sphere_data
number_plane_boundaries
2
layer_thickness
2.d0
layer_ref_index
(1.d0,0.d0),(1.54d0,0.d0),(1.d0,0.d0)
end_of_options
```

No sphere may cross a boundary.

### 5.4 Periodic lattices

Add the following to repeat the bubble example with a $3\times3$ unit cell:

```text
periodic_lattice
t
cell_width
3.d0,3.d0
```

Cell widths must prevent overlap across cell edges. Periodic output includes lattice reflectance, absorptance, and transmittance; scattering matrices are reported at reciprocal-lattice directions.

### 5.5 Multiple and looped runs

`new_run` retains current values, applies subsequent changes, and appends results unless a new output name is supplied:

```text
new_run
incident_beta_deg
40.d0
length_scale_factor
2.d0
end_of_options
```

`loop_variable` accepts an option name followed by start, stop, and increment:

```text
loop_variable
incident_beta_deg
0.d0,180.d0,2.d0
end_of_options
```

Up to three loop levels may be nested. The first declaration is outermost.

### 5.6 Accelerated translation

For large systems with no plane boundaries or periodicity, set `fft_translation_option=t`. The convolution method can reduce the external-sphere interaction cost from approximately $N_S^2$ to $N_S\log N_S$. It works best for clusters with a reasonably uniform concentration. Compare an accelerated run with a conventional run when establishing tolerances.

## 6. Input option reference

Defaults are shown in parentheses where established by the current source.

### 6.1 System definition

| Option | Meaning |
|---|---|
| `number_spheres` | Number of sphere surfaces, $N_S$. |
| `sphere_data_input_file` | Position/radius/index file. Only the first $N_S$ records are read. |
| `length_scale_factor` | Multiplier for sphere geometry, layer thickness, and cell width (`1`). For geometry in micrometres, use $2\pi/\lambda_0$. |
| `ref_index_scale_factor` | Multiplier for indices read from sphere data; also supplies the common sphere index when none is present (`(1,0)`). |
| `x_shift`, `y_shift`, `z_shift` | Position shifts, multiplied by `length_scale_factor`. |
| `shifted_sphere` | Sphere to shift; zero shifts every sphere. Useful as a loop variable. |
| `check_positions` | Check sphere overlap and geometry (`t`). |
| `number_plane_boundaries` | Number of boundaries, $N_B$ (`0`, maximum `10`). |
| `layer_thickness` | $N_B-1$ layer thicknesses. |
| `layer_ref_index` | $N_B+1$ complex layer indices from bottom to top (`(1,0)`). |
| `medium_ref_index` | Embedding-medium index for generated configurations. |
| `periodic_lattice` | Enable infinite repetition in $x,y$ (`f`). |
| `cell_width` | One value for a square cell or two values for $W_x,W_y$. |
| `finite_lattice` | Use finite-lattice behavior instead of infinite periodic repetition. |
| `gaussian_beam_constant` | Inverse focal width; zero selects a plane wave (`0`). The localized approximation is intended for values up to about `0.2`. |
| `gaussian_beam_focal_point` | Unscaled focal-point coordinates (`0,0,0`). |

### 6.2 Numerical solution

| Option | Meaning |
|---|---|
| `mie_epsilon` | VSWF truncation tolerance (`1e-6`). A negative integer `-L` forces order `L`. |
| `solution_epsilon` | Iterative-solver error tolerance (`1e-6`). |
| `max_iterations` | Maximum solver iterations (`10000` in the current implementation). |
| `solution_method` | Solution algorithm (`iteration`). |
| `translation_epsilon` | Common-origin and T-matrix truncation tolerance (`1e-5`). |
| `max_t_matrix_order` | Safety limit on T-matrix/common-origin order (`120`). |
| `t_matrix_convergence_epsilon` | Relative extinction-change threshold for T-matrix convergence (`1e-6`). |
| `normalize_solution_error` | Normalize the reported iterative error. |
| `store_translation_matrix` | Cache translation matrices when memory permits. |
| `store_surface_matrix` | Cache plane-boundary surface matrices. |

Always repeat sensitive calculations with tighter tolerances. High expansion orders may encounter floating-point range limits.

### 6.3 FFT acceleration

| Option | Meaning |
|---|---|
| `fft_translation_option` | Enable FFT translation (`f`). |
| `min_fft_nsphere` | Minimum sphere count at which FFT translation activates (`200`). |
| `cell_volume_fraction` | Controls cubic-node spacing; zero selects an automatic estimate (`0`). |
| `d_cell` | Explicit node-cell size, overriding automatic selection. |
| `node_order` | Node expansion order; `-L` selects approximately $\lceil d\rceil+L$ (`-1`). |
| `neighbor_node_model` | Select treatment of near-node interactions. |
| `interaction_radius` | Limit for explicitly stored/interacted neighboring spheres. |

### 6.4 Files and reporting

| Option | Meaning |
|---|---|
| `output_file` | Final results file (`mstmtest.dat`). |
| `append_output_file` | Append instead of replacing an existing result file (`f`). |
| `run_file` | Intermediate progress file; blank uses standard output. |
| `copy_input_file` | Include input content in reporting (`f`). |
| `print_sphere_data` | Print individual sphere input and calculated properties (`t`). |
| `print_timings` | Print timing information (`t`). |
| `normalize_s11` | Normalize scattering-matrix output (`t`). |

### 6.5 Orientation and scattering maps

| Option | Meaning |
|---|---|
| `random_orientation` | Exact analytical orientation average (`f`); incompatible with boundaries and periodicity. |
| `incidence_average` | Monte Carlo average over incident directions (`f`). |
| `number_incident_directions` | Monte Carlo direction count (`16`). |
| `incident_alpha_deg` | Incident azimuth $\alpha$ in degrees (`0`). |
| `incident_beta_deg` | Incident polar angle $\beta$ in degrees (`0`). |
| `incident_sin_beta` | Generalized incidence, including evanescent excitation when greater than one. |
| `incident_direction` | Select propagation direction relative to layered geometry. |
| `incident_frame` | Base scattering angles on the incident rather than target frame (`f`). |
| `calculate_scattering_matrix` | Calculate the Stokes matrix (`t`). |
| `scattering_map_model` | `0`: angular cut; `1`: two-dimensional $(k_x,k_y)$ maps (`0`). |
| `scattering_map_increment` | Angular increment in degrees for model 0 (`1`). |
| `scattering_map_dimension` | Half-grid dimension for model 1 (`15`). |
| `azimuthal_average` | Analytically average over scattering azimuth (`f`). |
| `numerical_azimuthal_average` | Use numerical rather than analytical azimuthal averaging (`f`). |
| `single_origin_expansion` | Post-process scattering from one origin (`t`); required by analytical azimuthal averaging. |
| `calculate_up_down_scattering` | Report hemispherical up/down scattering (`t`). |

### 6.6 Near fields

| Option | Meaning |
|---|---|
| `calculate_near_field` | Enable electric and magnetic field maps (`f`); requires fixed orientation. |
| `near_field_minimum_border` | Minimum $(x,y,z)$ grid corner. |
| `near_field_maximum_border` | Maximum corner. Equal min/max values collapse an axis. |
| `near_field_step_size` | Grid spacing (`0.2`). |
| `near_field_calculation_model` | `1`: external total field and internal field; other values remove the external incident contribution (`1`). |
| `store_surface_vector` | Cache boundary/periodic contributions for faster evaluation (`t`). |
| `near_field_expansion_spacing` | Re-expansion grid spacing (`5`). |
| `near_field_expansion_order` | Re-expansion truncation order (`10`). |
| `near_field_output_file` | Field-map output (`nftest.dat`). |

The near-field file starts with run and geometry metadata, then grid dimensions. Each point has coordinates followed by complex electric and magnetic vector components for parallel and perpendicular incident polarizations (27 real columns total).

### 6.7 Random configuration generation

Generated radii have volume mean one before `length_scale_factor` is applied. A log-normal distribution is used:

$$P(a)=\frac{1}{\sigma a\sqrt{2\pi}}\exp\left[-\frac{(3\sigma^2+2\ln a)^2}{8\sigma^2}\right],\qquad
\int_0^\infty a^3P(a)\,da=1.$$

The generator starts from random sequential placement or a close-packed lattice, then simulates elastic hard-sphere motion.

| Option | Meaning |
|---|---|
| `random_configuration` | Generate a configuration (`f`). |
| `target_shape` | `0`: rectangular; `1`: cylindrical; `2`: spherical. |
| `target_dimensions` | Three half-dimensions. For a cylinder, the first is radius and third is half-thickness; for a sphere, the first is radius. |
| `target_width`, `target_thickness` | Shorthand width/radius and half-thickness controls. |
| `auto_target_radius` | Derive a spherical target radius from count and volume fraction (`f`). |
| `target_radius_padding` | Padding used when sampling an automatically sized target (`5`). |
| `sphere_volume_fraction` | Target particle volume fraction. High values (roughly above `0.4–0.5`) may fail to pack. |
| `periodic_bc` | Three logical flags for periodic target boundaries. Only flat boundaries can be periodic. |
| `wall_boundary_model` | Select hard-wall handling. |
| `psd_sigma` | Log-normal size-distribution width; zero is monodisperse. |
| `max_diffusion_simulation_time` | Maximum simulated displacement time (`5`). |
| `max_collisions_per_sphere` | Stop after this mean collision count (`3`). |
| `max_diffusion_cpu_time` | Wall-clock limit in seconds (`100`). |
| `print_random_configuration` | Write the sampled spheres (`f`). |
| `random_configuration_output_file` | Configuration output (`random_configuration.pos`). |
| `frozen_configuration` | Reuse rather than resample a configuration (`f`). |
| `number_components` | Number of polydisperse/material components (up to four). |
| `component_radii` | Component mean radii. |
| `component_number_fraction` | Component number fractions. |
| `component_ref_index` | Complex component indices. |

### 6.8 Configuration averaging and effective media

`configuration_average=t` generates multiple configurations and accumulates total and coherent fields. The diffuse scattering matrix is the averaged total matrix minus the matrix of the averaged (coherent) field. The same workflow can estimate particle-layer reflectance/absorptance and effective extinction or absorption.

| Option | Meaning |
|---|---|
| `configuration_average` | Enable averaging (`f`). |
| `number_configurations` | Number of independent samples. |
| `reflection_model` | Configure particle-layer reflection/absorption reporting (`f`). |
| `effective_medium_simulation` | For averaged spherical targets without a host sphere, fit a homogeneous effective sphere (`f`). |
| `fit_for_radius` | Fit effective radius as well as complex refractive index (`t`). |
| `random_configuration_host` | Add an enclosing host sphere to generated geometry (`f`). |
| `host_sphere_ref_index` | Complex index of the generated host sphere (`(1,0)`). |
| `random_configuration_host_model` | Select how the effective/host radius is defined (`1`). |
| `absorption_sample_radius` | Explicit radius used for absorption sampling. |
| `auto_absorption_sample_radius` | Derive the sampling radius from the target (`t`). |
| `absorption_sample_radius_fraction` | Automatic radius as a fraction of target radius (`0.8`). |

Effective-medium output includes the fitted complex index, fitted radius, an independent-scattering/radiative-transfer ratio, and fit status. When an averaged spherical target also requests near fields, a field-based effective index is reported. Treat these as model-dependent estimates and verify convergence with sphere count, configuration count, and target size.

### 6.9 Run control

| Option | Meaning |
|---|---|
| `new_run` | Execute the accumulated options, then begin another run retaining current settings. |
| `loop_variable` | Sweep the next named numerical parameter over start, stop, increment. |
| `end_of_options` | Execute the final run and stop reading input. Required at end of file. |

## 7. Output interpretation

The volume-equivalent cluster radius and area-mean sphere radius are

$$R_V=\left(\sum_i a_i^3\right)^{1/3},\qquad
R_A=\left(\frac{1}{N_{S,\mathrm{ext}}}\sum_{i\in\mathrm{ext}}a_i^2\right)^{1/2}.$$

The cross-section radius $a_{cs}$ defines efficiencies through $C=\pi a_{cs}^2Q$. Dimensional cross sections are obtained by dividing dimensionless values by $(2\pi/\lambda_0)^2$. For plane-wave or random-orientation calculations, $a_{cs}=R_V$. For Gaussian incidence, $a_{cs}=1/(\sqrt2 C_B)$. For periodic lattices, $\pi a_{cs}^2=W_xW_y$.

When sphere data is printed, `absorption` measures all absorption inside a surface. `volume absorption` subtracts absorption in nested spheres and therefore represents only the material immediately inside that surface.

Simple systems report extinction, absorption, and scattering efficiencies for unpolarized, parallel, and perpendicular incidence. Layered systems additionally separate forward and backward terms. Periodic systems report reflection, absorption, and transmission, which should sum to unity within numerical error.

With normalization enabled,

$$\int_0^{2\pi}\int_0^\pi S_{11}(\theta,\phi)\sin\theta\,d\theta\,d\phi=1.$$

For random orientation or analytical azimuthal averages, $S_{11}$ is independent of $\phi$ and $\int_0^\pi S_{11}(\theta)\sin\theta\,d\theta=1/(2\pi)$. Other matrix elements are scaled relative to $S_{11}$.

## 8. References

1. D. W. Mackowski and M. I. Mishchenko, “Calculation of the T matrix and the scattering matrix for ensembles of spheres,” *JOSA A* 13 (1996), 2266–2278.
2. D. W. Mackowski, “Exact solution for the scattering and absorption properties of sphere clusters on a plane surface,” *JQSRT* 109 (2007), 770–788.
3. D. W. Mackowski and M. I. Mishchenko, “Direct simulation of multiple scattering by discrete random media illuminated by Gaussian beams,” *Physical Review A* 83 (2011), 013804.
4. D. W. Mackowski, “The extension of Mie theory to multiple spheres,” in *The Mie Theory*, Springer (2012), 223–256.
5. D. W. Mackowski, “A general superposition solution for electromagnetic scattering by multiple spherical domains of optically active media,” *JQSRT* 133 (2014), 264–270.
6. D. W. Mackowski and L. Kolokolova, “Application of the multiple sphere superposition solution to large-scale systems of spheres via an accelerated algorithm,” 2022. See [the repository paper](reference/mackowski-mstm-2022b.pdf).

The archived 2023 PDF manual remains in [`docs/reference/mstm-manual-2021.pdf`](reference/mstm-manual-2021.pdf) for provenance.
