"""
SI tables S1 and S2 -- print the LaTeX table rows from data/csv/ribbon_validation.csv.

Table S1: L = 1 nanographene, G0W0@PBE, FHI-aims CD and AC vs CP2K RI-RS,
          def2-TZVP and aug-cc-pVDZ.  LUMO = QP energy of the KS LUMO orbital.
Table S2: L = 1, 2, 4, 8, G0W0+H@PBE, aug-cc-pVDZ, FHI-aims AC vs CP2K RI-RS.
          LUMO = "window" LUMO in both codes (lowest QP level of the virtual
          states within 2 eV of the KS LUMO, see make_tables.py), because for
          L = 1, 2 a diffuse state drops below the KS LUMO in both codes.

Copy the printed rows into si.tex.
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/ribbon_validation.csv")


def f(x):
    """Energy in eV with three decimals, minus sign in math mode, '--' if missing."""
    if np.isnan(x):
        return "--"
    s = f"{x:.3f}"
    return f"${s}$" if x < 0 else s


print("% Table S1")
for basis, name in (("tzvp", "def2-TZVP"), ("avdz", "aug-cc-pVDZ")):
    r = d[(d.L == 1) & (d.basis == basis)].iloc[0]
    for tag, lab in (("aims_cd", "FHI-aims CD"), ("aims_ac", "FHI-aims AC"), ("rirs", "CP2K RI-RS")):
        h, l = r[f"homo_{tag}"], r[f"lumo_{tag}"]
        print(f"{name} & {lab} & {f(h)} & {f(l)} & {f(l - h)} \\\\")

print("% Table S2")
for L in (1, 2, 4, 8):
    r = d[(d.L == L) & (d.basis == "avdz_hedin_shift")].iloc[0]
    ha, la, hr, lr = r.homo_aims_win, r.lumo_aims_win, r.homo_rirs_win, r.lumo_rirs_win
    print(f"{L} & {f(ha)} & {f(la)} & {f(la - ha)} & {f(hr)} & {f(lr)} & {f(lr - hr)} \\\\")
