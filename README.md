# Data for "Low-Scaling G0W0 with Real-Space Resolution of the Identity for Thousands of Atoms"

R. Tyagi, M. Graml, J. Wilhelm
Institute of Theoretical Physics and Regensburg Center for Ultrafast Nanoscopy (RUN),
University of Regensburg, Germany

This archive contains the atomic geometries, the input and output files of all
calculations in the paper and its Supporting Information (SI), the numbers
extracted from them (CSV files), and the scripts that produce every figure and
SI table.

## Contents

```
README.md
data/
  gw100/                 GW100 benchmark, G0W0@PBE0/def2-TZVP
  acc24/                 ACC24 benchmark, G0W0@PBE0/def2-TZVP
  nanographene/          nanographenes C14L H4L+8
  csv/                   numbers extracted from the outputs by script/make_tables.py
script/                  make_tables.py, si_tables.py, plot_param.py
fig/                     all figures of the paper and the SI
```

## Data layout

### GW100 (`data/gw100/`)

93 of the 100 GW100 molecules (all molecules with elements up to Kr; the seven
molecules with fifth-row elements are excluded because the real-space grids of
Duchemin and Blase are available only up to Kr).

```
gw100/
  cp2k_rirs/<id>/        CP2K RI-RS
  cp2k_tensor/<id>/      CP2K tensor-based low-scaling GW (reference for RI-RS)
  fhi_aims/<id>/         FHI-aims, contour deformation (CD) (reference)
  cp2k_rirs/BASIS_def2_tzvp_all, cp2k_tensor/BASIS_def2_tzvp_all   basis set file used by all runs
  cp2k_rirs/submit.slurm, cp2k_tensor/submit.slurm                job scripts
```

### ACC24 (`data/acc24/`)

24 organic acceptor molecules, same structure as GW100 (`cp2k_rirs/<id>`,
`cp2k_tensor/<id>`, `fhi_aims/<id>`), plus the scan of the Tikhonov parameter
lambda (Figure S2):

```
acc24/
  cp2k_rirs_tikhonov/<id>/<k>/   RI-RS with lambda = 10^-k, k = 2, 4, 6, 10;
                                 k = 0 means no regularization (lambda = 0).
                                 lambda = 10^-8 is the default run in cp2k_rirs/<id>.
```

### Nanographenes (`data/nanographene/`)

Nanographenes C14L H4L+8 with L strips (N = 18 L + 8 atoms). The central carbon
atom of each zigzag edge carries two hydrogen atoms.

```
nanographene/
  <L>/                     L = 1, 2, 4, ..., 512: size scaling, G0W0+H@PBE (Figures 3, 4, 5)
    rirs_opt_mini/           RI-RS-Opt (r_AO = 6 A, r_W = 30 A, optimized grid), aug-SZV-MOLOPT-ae-mini, all L
    rirs_mini/               RI-RS, aug-SZV-MOLOPT-ae-mini, all L
    tensor_mini/             tensor-based GW, aug-SZV-MOLOPT-ae-mini, L = 1 ... 32
    rirs_avdz/               RI-RS, aug-cc-pVDZ, L = 1 ... 16
    tensor_avdz/             tensor-based GW, aug-cc-pVDZ, L = 1 ... 4
  fhi_aims_cp2k_rirs/      L = 1, 2, 4, 8: comparison with FHI-aims (Figure 2c, Tables S1, S2)
    cp2k_rirs/<L>/<basis>/   CP2K RI-RS
    fhi_aims/<L>/<basis>/    FHI-aims, analytic continuation (AC)
    fhi_aims/1/<basis>/ac/, fhi_aims/1/<basis>/cd/
                             L = 1 only: FHI-aims with AC and with CD (Table S1)
      <basis> = tzvp               G0W0@PBE, def2-TZVP
                avdz               G0W0@PBE, aug-cc-pVDZ
                avdz_hedin_shift   G0W0+H@PBE, aug-cc-pVDZ (Table S2)
  parameter_convergence_opt_mini/<scan>/<value>/
                           L = 16 (296 atoms), RI-RS-Opt, aug-SZV-MOLOPT-ae-mini (Figure S1).
                           One parameter is varied at a time; all others are kept at
                           r_c = 14 A, r_AO = 6 A, r_W = 30 A, EPS_FILTER = 1e-10,
                           RS_AO_RATIO = 6, CUTOFF_ATOMIC_CLUSTER = 3 A.
      <scan> (CP2K keyword)                         <value>
      cutoff_radius_ri      (CUTOFF_RADIUS_RI, r_c)       8, 10, 12, 14, 16 (A)
      cutoff_radius_rl_ao   (CUTOFF_RADIUS_RL_AO, r_AO)   4, 5, 6, 7 (A)
      cutoff_radius_rl_w    (CUTOFF_RADIUS_RL_W, r_W)     20, 25, 30, 35 (A)
      eps_filter            (EPS_FILTER)                  k for 10^-k: 6, 8, 10, 12
      rs_ao_ratio           (RS_AO_RATIO)                 4 ... 8 grid points per AO
      cutoff_atomic_cluster (CUTOFF_ATOMIC_CLUSTER)       1, 2, 3, 4 (A), cluster radius R_RS
```

### Files in a run folder

| File | Content |
|---|---|
| `input.inp` | CP2K input |
| `output.out` | CP2K output (contains the git revision, the RI-RS statistics and the timers) |
| `bandstructure_SCF_and_G0W0` | CP2K: per-level DFT and G0W0 energies (columns n, occ/vir, k, e_DFT, Sigma_c, Sigma_x, v_xc, e_G0W0, in eV) |
| `struc.xyz` | geometry (Angstrom) |
| `BASIS_*_all` | Gaussian basis sets (orbital and RI auxiliary) used by the run |
| `submit.slurm` | job script (nodes, MPI processes, OpenMP threads) |
| `control.in`, `geometry.in` | FHI-aims input |
| `submit.aims` | FHI-aims job script (GW100, ACC24) |

## Calculation details

### CP2K

- CP2K development version 2026.2, git revision **0e87843579** for all CP2K
  calculations.
- All-electron, non-periodic (`PERIODIC NONE`, multipole Poisson solver),
  GAPW scheme, plane-wave cutoff 600 Ry, relative cutoff 60 Ry.
- G0W0: 30 minimax time and frequency points, filter threshold
  EPS_FILTER = 1e-10 (atomic units), truncated Coulomb metric with
  r_c = CUTOFF_RADIUS_RI = 14 A, Pade analytic continuation with at most 16
  parameters fitted to all imaginary frequencies below 10 eV, quasiparticle
  equation solved with Newton's method.
- He and H2 of GW100: REGULARIZATION_MINIMAX 0.0, because their minimax grids
  span a very small energy range.
- RI-RS: atom-centered grids of Duchemin and Blase (optimized for def2-TZVP),
  Tikhonov parameter lambda = 1e-8. For GW100, ACC24 and the comparison with
  FHI-aims, no fixed r_AO, no r_W and no grid optimization are used.
- RI-RS-Opt (size scaling): r_AO = CUTOFF_RADIUS_RL_AO = 6 A,
  r_W = CUTOFF_RADIUS_RL_W = 30 A, RS_AO_RATIO = 6 grid points per AO for every
  atom (54 for C, 36 for H with aug-SZV-MOLOPT-ae-mini), grid positions
  optimized with L-BFGS-B on atomic clusters of radius CUTOFF_ATOMIC_CLUSTER = 3 A.
- Tensor-based GW: same DFT, basis sets, RI metric, time/frequency grids and
  filter threshold as RI-RS.
- Benchmarks: G0W0@PBE0 (25 % exact exchange), def2-TZVP with def2-TZVP-RIFIT.
- Nanographenes: G0W0@PBE (comparison with FHI-aims) and G0W0+H@PBE with
  Hedin's shift (`HEDIN_SHIFT`, size scaling); aug-cc-pVDZ with
  aug-cc-pVDZ-RIFIT, aug-SZV-MOLOPT-ae-mini with RI_aug-SZV-MOLOPT-ae-mini_tier1.
- PBE ground state of the nanographenes: dense diagonalization (ELPA) with
  Fermi-Dirac smearing, Broyden mixing, 100 added MOs, EPS_SCF = 1e-7.

### FHI-aims

- FHI-aims version 260530, the same Gaussian basis sets as CP2K (given in
  `control.in`), RI-V auxiliary basis generated by FHI-aims
  (`prodbas_threshold 1e-5`), non-relativistic all-electron Hamiltonian.
- GW100 and ACC24: contour deformation with 600 frequency points for the
  frontier levels (window `contour_def_gw`).
- Nanographenes: analytic continuation with a Pade model with 16 parameters
  from 100 imaginary frequencies. For L = 1, additionally contour deformation
  with 600 frequency points (`fhi_aims/1/<basis>/cd/`); CD and AC agree within
  0.2 meV (Table S1).

### Hardware

Noctua 2 (Paderborn Center for Parallel Computing), two AMD EPYC 7763
processors (128 cores) per node. The size-scaling runs use 8 MPI processes with
16 OpenMP threads per node on 1 to 18 nodes, with more nodes for larger L
(columns `mpi`, `threads`, `cores`, `nodes` of `ribbon_scaling.csv`). Costs are
given in core hours = wall time x number of cores.
