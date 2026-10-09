# Data for "Low-Scaling G0W0 with Real-Space Resolution of the Identity: Benchmark Accuracy for Systems with Thousands of Atoms"

R. Tyagi, M. Graml, R. Pasquier, J. Wilhelm

Atomic geometries, input and output files of all calculations in the paper,
the extracted tables, and the scripts that produce all figures.

## Layout

```
data/
  gw100/          GW100, G0W0@PBE0/def2-TZVP, one folder per molecule (GW100 index)
    cp2k_rirs/      CP2K RI-RS
    cp2k_tensor/    CP2K tensor-based low-scaling GW
    fhi_aims/       FHI-aims reference (contour deformation)
  acc24/          ACC24, same structure, plus
    cp2k_rirs_tikhonov/<k>/   RI-RS with Tikhonov parameter 10^-k (0 = no regularization)
  nanographene/   nanographenes C14L H4L+8
    <L>/            L = 1 ... 512 strips: size scaling, G0W0+H@PBE
      rirs_opt_mini   RI-RS-Opt (cutoffs + optimized grid), aug-SZV-MOLOPT-ae-mini
      rirs_mini       RI-RS, aug-SZV-MOLOPT-ae-mini
      tensor_mini     tensor-based GW, aug-SZV-MOLOPT-ae-mini
      rirs_avdz       RI-RS, aug-cc-pVDZ
      tensor_avdz     tensor-based GW, aug-cc-pVDZ
    fhi_aims_cp2k_rirs/   L = 1, 2, 4, 8: CP2K RI-RS vs FHI-aims (def2-TZVP, aug-cc-pVDZ, with/without Hedin's shift)
    parameter_convergence_opt_mini/   L = 16: scans of the numerical parameters
    parallel_scaling/                 L = 64: MPI / OpenMP scans (not shown in the paper)
  csv/            numbers extracted by script/make_tables.py
script/
  make_tables.py  reads all outputs and writes data/csv/*.csv
  plot_param.py   shared plot settings (imported by every figure script)
  fig2_validation.py ... fig5_scf_vs_gw.py, figS1_convergence.py,
  figS2_tikhonov.py   one script per figure, read data/csv, write fig/
  fig1_workflow.tex   Figure 1 (TikZ): pdflatex -output-directory=../fig fig1_workflow.tex
  si_tables.py    prints the rows of SI Tables S1 and S2
fig/              all figures of the paper and the SI (PDF, PNG)
```

Every CP2K run folder contains `input.inp`, `output.out` and, where available,
the per-level quasiparticle energies `bandstructure_SCF_and_G0W0`.

## Reproducing the tables and figures

```
cd script
python make_tables.py          # data/  -> data/csv/*.csv   (numpy, pandas)
python fig2_validation.py      # data/csv -> fig/            (matplotlib)
...                            # same for the other figure scripts
```

The definitions used for HOMO and LUMO are documented at the top of
`make_tables.py`.

## Notes

- GW100 molecule 37 (CBr4): the FHI-aims output file is not included; its
  HOMO and LUMO (-10.0341 and -0.5236 eV) are entered in `make_tables.py`
  (`AIMS_MANUAL`).
- Calculations: CP2K development version 2026.2 (git revision 0005d2399e for
  the nanographenes), FHI-aims version 260530, Noctua 2 (PC2, Paderborn).
