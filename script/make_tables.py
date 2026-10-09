"""
Read every CP2K and FHI-aims output under data/ and write the numbers the
figures and tables need into data/csv/.  Run this once after new data arrive;
the figure scripts only read the CSV files.

Written CSV files:
    gw100.csv              HOMO/LUMO per molecule: FHI-aims, CP2K tensor, CP2K RI-RS
    acc24.csv              same for ACC24, plus RI-RS for every Tikhonov parameter
    ribbon_validation.csv  L = 1, 2, 4, 8 ribbons: FHI-aims (AC, CD) and CP2K RI-RS
    ribbon_scaling.csv     L = 1 ... 512: gap, basis sizes, grid size and timings
    parallel_scaling.csv   L = 64 ribbon on 8 nodes, varying MPI ranks / OMP threads
    convergence.csv        L = 16 ribbon, one numerical parameter scanned at a time

Two definitions of the LUMO are used, both applied identically to CP2K and
FHI-aims:

(1) "orbital" LUMO: the QP energy of the KS LUMO orbital.  Used for the
    molecular benchmarks (GW100, ACC24) and the def2-TZVP / aug-cc-pVDZ
    nanographene comparison, i.e. the same orbital is compared in both codes.

(2) "window" LUMO: the lowest QP energy among the virtual KS states whose KS
    energy lies at most WINDOW_EV = 2 eV above the KS LUMO.  Used for the
    HOMO-LUMO gaps of the size scaling and for the G0W0+H comparison.
    Why a window and not simply the lowest virtual QP level:
      - for L = 1, 2 a diffuse state 1.2-1.8 eV above the KS LUMO becomes the
        lowest QP level (G0W0 reorders the empty levels).  This is physical:
        FHI-aims finds the same level (L = 1, +H, aug-cc-pVDZ: 0.783 vs
        0.776 eV), and its Sigma_c is ordinary (about -0.8 eV);
      - for aug-cc-pVDZ and L = 8, 16, about 200 near-continuum states
        2.3 eV and more above the KS LUMO get Sigma_c of -4 to -9 eV and QP
        energies below the HOMO (down to -27 eV): failures of the analytic
        continuation.  The KS LUMO itself is fine (L = 8, +H: -1.195 eV vs
        -1.192 eV in FHI-aims).
    The 2 eV window keeps the first and drops the second in every run here.
    Note that the "G0W0 conduction band minimum" printed in output.out is
    the lowest virtual QP level WITHOUT a window and is therefore not used.

CP2K levels come from bandstructure_SCF_and_G0W0 (columns: n, (occ|vir), k,
e_DFT, Sigma_c, Sigma_x, v_xc, e_G0W0).  FHI-aims levels come from the LAST
table "GW quasi-particle energy levels" in output.out (the contour-deformation
runs print one table per self-consistency step of the QP equation).

The G0W0+H nanographene runs in fhi_aims_cp2k_rirs/cp2k_rirs/<L>/avdz_hedin_shift
are identical (same input, geometry and bandstructure file) to
nanographene/<L>/rirs_avdz; the latter is used as a fallback if a file is missing.
"""

import os
import re
import sys
import glob
import numpy as np
import pandas as pd

work_dir = sys.path[0]
data_dir = os.path.join(work_dir, "..", "data")
csv_dir = os.path.join(data_dir, "csv")
os.makedirs(csv_dir, exist_ok=True)

########## Readers ##########


WINDOW_EV = 2.0     # KS-energy window above the KS LUMO for the "window" LUMO


def pick_homo_lumo(levels, mode):
    """HOMO and LUMO from a list of (occupied?, e_KS, e_QP) tuples.

    mode = "orbital": LUMO = QP energy of the KS LUMO;
    mode = "window":  LUMO = lowest QP energy of the virtual states with
                      e_KS <= e_KS(LUMO) + WINDOW_EV (see the module docstring).
    """
    occ = [l for l in levels if l[0]]
    vir = [l for l in levels if not l[0]]
    if not occ or not vir:
        return np.nan, np.nan
    homo = occ[-1][2]
    if mode == "orbital":
        return homo, vir[0][2]
    e_lumo_ks = vir[0][1]
    return homo, min(l[2] for l in vir if l[1] <= e_lumo_ks + WINDOW_EV)


def read_cp2k_qp(run_dir, mode="orbital"):
    """Return (HOMO, LUMO) G0W0 energies in eV of a CP2K run, or (nan, nan).

    Read from bandstructure_SCF_and_G0W0; "(occ)" / "(vir)" is the KS
    occupation.  The last column is the QP energy (with Hedin's shift if the
    run used HEDIN_SHIFT).
    """
    path = os.path.join(run_dir, "bandstructure_SCF_and_G0W0")
    if not os.path.isfile(path):
        return read_cp2k_printed_edges(os.path.join(run_dir, "output.out"))
    levels = []
    for line in open(path):
        m = re.match(r"\s*\d+\s+\((occ|vir)\)\s+\d+\s+(.*)", line)
        if m:
            v = m.group(2).split()
            levels.append((m.group(1) == "occ", float(v[0]), float(v[-1])))
    return pick_homo_lumo(levels, mode)


def read_cp2k_printed_edges(out_file):
    """Fallback without bandstructure_SCF_and_G0W0: band edges from output.out.

    Only trusted if CP2K did not warn that "G0W0 reorders the empty levels";
    then the printed CBM is the QP energy of the KS LUMO and all three LUMO
    definitions coincide.  Otherwise nan.  (Used for the parallel-scaling and
    convergence runs, L = 64 and 16, none of which reorders.)
    """
    if not os.path.isfile(out_file):
        return np.nan, np.nan
    text = open(out_file, errors="replace").read()
    if "reorders the" in text:
        return np.nan, np.nan
    homo = re.findall(r"G0W0 valence band maximum \(eV\):\s+(\S+)", text)
    lumo = re.findall(r"G0W0 conduction band minimum \(eV\):\s+(\S+)", text)
    if not homo or not lumo:
        return np.nan, np.nan
    return float(homo[-1]), float(lumo[-1])


def read_aims_qp(out_file, mode="orbital"):
    """Return (HOMO, LUMO) quasiparticle energies in eV of an FHI-aims run."""
    if not os.path.isfile(out_file):
        return np.nan, np.nan
    lines = open(out_file, errors="replace").readlines()
    start = [i for i, l in enumerate(lines)
             if "GW quasi-particle energy levels" in l]
    if not start:
        return np.nan, np.nan
    levels = []
    for l in lines[start[-1]:]:
        p = l.split()
        # table rows: state, occ_num, e_gs, e_x^ex, e_xc^gs, e_c^nloc, e_qp
        if len(p) == 7 and p[0].isdigit():
            levels.append((float(p[1]) > 1e-3, float(p[2]), float(p[6])))
        if "Quasiparticle HOMO level" in l:
            break
    return pick_homo_lumo(levels, mode)


def read_timers(out_file):
    """Return {timer name: TOTAL TIME (maximum over ranks) in s} of a CP2K run."""
    t = {}
    inside = False
    for l in open(out_file, errors="replace"):
        if "T I M I N G" in l:
            inside = True
            continue
        if inside:
            if "PROGRAM ENDED" in l:
                break
            p = l.split()
            # SUBROUTINE  CALLS  ASD  SELF(avg)  SELF(max)  TOTAL(avg)  TOTAL(max)
            if len(p) == 7:
                try:
                    t[p[0]] = float(p[6])
                except ValueError:
                    pass
    return t


def grep_last_number(out_file, key, cast=float):
    """Last number on the last line of out_file that contains `key`.

    A trailing "%" (occupation lines) is ignored.
    """
    val = np.nan
    for l in open(out_file, errors="replace"):
        if key in l:
            p = l.replace("%", " ").split()
            try:
                val = cast(p[-1])
            except ValueError:
                pass
    return val


def read_run_info(out_file):
    """Size of the run (atoms, basis, grid) and of the machine (ranks, threads)."""
    return {
        "n_atoms": grep_last_number(out_file, "- Atoms:", int),
        "n_ao": grep_last_number(out_file, "Number of Gaussian basis functions for MOs", int),
        "n_ri": grep_last_number(out_file, "Number of auxiliary Gaussian basis functions", int),
        "n_grid": grep_last_number(out_file, "Total grid points used for RI-RS", int),
        "mpi": grep_last_number(out_file, "Total number of message passing processes", int),
        "threads": grep_last_number(out_file, "Number of threads for this process", int),
        "occ_phi": grep_last_number(out_file, "Percentage of non-zero matrix elements in Φ_μ(r_l)"),
        "occ_Z": grep_last_number(out_file, "Percentage of non-zero matrix elements in Z_lP"),
        # one sparsity pattern shared by the grid-by-grid matrices G_ll', chi_ll', W_ll'
        "occ_Gll": grep_last_number(out_file, "Percentage of non-zero matrix elements in G_ll'"),
        # chi_PQ is printed for the first time point only
        "occ_chiPQ": grep_last_number(out_file, "Percentage of non-zero matrix elements in χ_PQ"),
        "revision": open(out_file, errors="replace").read().split(
            "source code revision number:")[1].split()[0]
        if "source code revision number:" in open(out_file, errors="replace").read() else "",
    }


# CP2K timers of the GW steps.  'scf' is everything before the band-structure
# post-processing, i.e. the PBE ground state including the setup of the
# Hamiltonian; 'total' is the full GW part.
STEP_TIMERS = {
    'total':   'gw',
    'grid':    'ri_rs_grid_assembler',      # includes the grid optimisation
    'Z_lP':    'compute_Z_lP',
    'chi':     'get_mat_chi_Gamma_tau',
    'W':       'compute_W',
    'Sigma_x': 'compute_Sigma_x',
    'Sigma_c': 'compute_Sigma_c',
    'QP':      'compute_QP_energies',
}


def step_times(out_file):
    """Wall time in s of each GW step, plus the SCF and the whole run."""
    t = read_timers(out_file)
    row = {f"t_{k}": t.get(v, np.nan) for k, v in STEP_TIMERS.items()}
    row["t_scf"] = t.get("qs_energies", np.nan) - t.get("qs_energies_properties", 0.0)
    row["t_cp2k"] = t.get("CP2K", np.nan)
    return row


########## GW100 and ACC24 ##########

def formula(xyz):
    """Sum formula from an xyz file, elements in order of first appearance."""
    sym = [l.split()[0] for l in open(xyz).readlines()[2:] if l.strip()]
    out = ""
    for s in dict.fromkeys(sym):
        n = sym.count(s)
        out += s + (str(n) if n > 1 else "")
    return out


# FHI-aims G0W0@PBE0 values (eV) given by hand where the output file is not
# (yet) in data/.  Remove an entry once its output.out is added; values from
# the file always take precedence.
#   GW100 37 = CBr4: FHI-aims numbers provided by Ritaj, 2026-10-08.
AIMS_MANUAL = {("gw100", 37): (-10.0341, -0.5236)}

for name in ("gw100", "acc24"):
    rows = []
    ids = sorted(int(d) for d in os.listdir(f"{data_dir}/{name}/fhi_aims"))
    for i in ids:
        row = {"id": i, "formula": formula(f"{data_dir}/{name}/fhi_aims/{i}/struc.xyz")}
        row["homo_aims"], row["lumo_aims"] = read_aims_qp(
            f"{data_dir}/{name}/fhi_aims/{i}/output.out")
        if np.isnan(row["homo_aims"]) and (name, i) in AIMS_MANUAL:
            row["homo_aims"], row["lumo_aims"] = AIMS_MANUAL[(name, i)]
        for code in ("cp2k_tensor", "cp2k_rirs"):
            tag = code.replace("cp2k_", "")
            row[f"homo_{tag}"], row[f"lumo_{tag}"] = read_cp2k_qp(
                f"{data_dir}/{name}/{code}/{i}")
        if name == "acc24":
            # Tikhonov scan: folder k means lambda = 10^-k, folder 0 means lambda = 0
            for k in (0, 2, 4, 6, 8, 10, 12):
                row[f"homo_tikh{k}"], row[f"lumo_tikh{k}"] = read_cp2k_qp(
                    f"{data_dir}/{name}/cp2k_rirs_tikhonov/{i}/{k}")
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(f"{csv_dir}/{name}.csv", index=False, float_format="%.4f")
    print(f"wrote data/csv/{name}.csv  ({len(df)} molecules, "
          f"{df.isna().any(axis=1).sum()} with a missing value)")

########## Ribbons: FHI-aims vs CP2K (L = 1, 2, 4, 8) ##########

base = f"{data_dir}/nanographene/fhi_aims_cp2k_rirs"
rows = []
for L in (1, 2, 4, 8):
    for basis in ("tzvp", "avdz", "avdz_hedin_shift"):
        row = {"L": L, "basis": basis}
        # FHI-aims: L = 1 was run with both analytic continuation (ac/) and
        # contour deformation (cd/); the longer ribbons with AC only.
        aims = f"{base}/fhi_aims/{L}/{basis}"
        if os.path.isdir(f"{aims}/ac"):
            row["homo_aims_ac"], row["lumo_aims_ac"] = read_aims_qp(f"{aims}/ac/output.out")
            row["homo_aims_cd"], row["lumo_aims_cd"] = read_aims_qp(f"{aims}/cd/output.out")
        else:
            row["homo_aims_ac"], row["lumo_aims_ac"] = read_aims_qp(f"{aims}/output.out")
            row["homo_aims_cd"], row["lumo_aims_cd"] = np.nan, np.nan
        cp2k = f"{base}/cp2k_rirs/{L}/{basis}"
        if not os.path.isfile(f"{cp2k}/bandstructure_SCF_and_G0W0"):
            # +H runs: identical to the size-scaling run nanographene/<L>/rirs_avdz
            cp2k = f"{data_dir}/nanographene/{L}/rirs_avdz"
        row["homo_rirs"], row["lumo_rirs"] = read_cp2k_qp(cp2k)
        # "window" LUMO in both codes (see module docstring)
        aims_out = f"{aims}/ac/output.out" if os.path.isdir(f"{aims}/ac") else f"{aims}/output.out"
        row["homo_aims_win"], row["lumo_aims_win"] = read_aims_qp(aims_out, mode="window")
        row["homo_rirs_win"], row["lumo_rirs_win"] = read_cp2k_qp(cp2k, mode="window")
        rows.append(row)
df = pd.DataFrame(rows)
for c in ("aims_ac", "aims_cd", "rirs", "aims_win", "rirs_win"):
    df[f"gap_{c}"] = df[f"lumo_{c}"] - df[f"homo_{c}"]
df.to_csv(f"{csv_dir}/ribbon_validation.csv", index=False, float_format="%.4f")
print("wrote data/csv/ribbon_validation.csv")

########## Ribbons: size scaling (L = 1 ... 512) ##########

rows = []
for L in (1, 2, 4, 8, 16, 32, 64, 128, 256, 512):
    for flavour in ("rirs_opt_mini", "rirs_mini", "tensor_mini", "rirs_avdz", "tensor_avdz"):
        run = f"{data_dir}/nanographene/{L}/{flavour}"
        out = f"{run}/output.out"
        if not os.path.isfile(out):
            continue
        row = {"L": L, "flavour": flavour}
        row.update(read_run_info(out))
        row["homo"], row["lumo"] = read_cp2k_qp(run, mode="window")
        row["homo_orb"], row["lumo_orb"] = read_cp2k_qp(run, mode="orbital")
        row["gap"] = row["lumo"] - row["homo"]
        row.update(step_times(out))
        rows.append(row)
df = pd.DataFrame(rows)
df["cores"] = df["mpi"] * df["threads"]
df["nodes"] = df["cores"] / 128           # Noctua 2: 128 cores per node

# Reference cost of the PBE ground state (SCF) for every L.  The RI-RS runs
# for large L restarted from an already converged wave function, so their
# own SCF time is not the cost of a ground-state calculation.  Sources (set
# by Ritaj, 2026-10-06):
#   L = 1 ... 64    tensor_mini/output.out   (SCF from scratch)
#   L = 128         rirs_mini/scf.out        (separate run, SCF part only)
#   L = 256         rirs_mini/output.out
#   L = 512         rirs_mini/scf.out        (separate SCF-only run of the
#                   RI-RS + cutoffs setup)
# SCF time = qs_energies - qs_energies_properties, i.e. everything before the
# band-structure / GW post-processing; core hours with the cores of that run.
SCF_SOURCE = {L: f"{L}/tensor_mini/output.out" for L in (1, 2, 4, 8, 16, 32, 64)}
SCF_SOURCE.update({128: "128/rirs_mini/scf.out", 256: "256/rirs_mini/output.out",
                   512: "512/rirs_mini/scf.out"})
for L, src in SCF_SOURCE.items():
    out = f"{data_dir}/nanographene/{src}"
    t = step_times(out)["t_scf"]
    info = read_run_info(out)
    m = df.L == L
    df.loc[m, "t_scf_ref"] = t
    df.loc[m, "cores_scf_ref"] = info["mpi"] * info["threads"]
    df.loc[m, "scf_ref_source"] = src
df.to_csv(f"{csv_dir}/ribbon_scaling.csv", index=False, float_format="%.4f")
print("wrote data/csv/ribbon_scaling.csv")

########## Parallel scaling (L = 64, 1160 atoms, 8 nodes) ##########

rows = []
for kind in ("mpi", "omp"):
    for d in sorted(glob.glob(f"{data_dir}/nanographene/parallel_scaling/{kind}/*")):
        out = f"{d}/output.out"
        row = {"scan": kind, "value": int(os.path.basename(d))}
        row.update(read_run_info(out))
        row["homo"], row["lumo"] = read_cp2k_qp(d, mode="window")
        row.update(step_times(out))
        rows.append(row)
df = pd.DataFrame(rows)
df["cores"] = df["mpi"] * df["threads"]
df = df.sort_values(["scan", "cores"])
df.to_csv(f"{csv_dir}/parallel_scaling.csv", index=False, float_format="%.4f")
print("wrote data/csv/parallel_scaling.csv")

########## Parameter convergence (L = 16, 296 atoms) ##########

# folder name -> (CP2K keyword, how the folder name maps to the value)
SCANS = {
    "cutoff_radius_ri":    ("CUTOFF_RADIUS_RI", lambda s: float(s)),
    "cutoff_radius_rl_ao": ("CUTOFF_RADIUS_RL_AO", lambda s: float(s)),
    "cutoff_radius_rl_w":  ("CUTOFF_RADIUS_RL_W", lambda s: float(s)),
    "eps_filter":          ("EPS_FILTER", lambda s: 10.0 ** (-int(s))),
    "grid_percentage":     ("N_POINT_PERCENTAGE", lambda s: float(s)),
}
rows = []
for scan, (keyword, to_value) in SCANS.items():
    for d in glob.glob(f"{data_dir}/nanographene/parameter_convergence_opt_mini/{scan}/*"):
        out = f"{d}/output.out"
        if not os.path.isfile(out):          # cutoff_radius_rl_ao/7 was never run
            continue
        row = {"scan": scan, "keyword": keyword,
               "value": to_value(os.path.basename(d))}
        row.update(read_run_info(out))
        row["homo"], row["lumo"] = read_cp2k_qp(d, mode="window")
        row["gap"] = row["lumo"] - row["homo"]
        row.update(step_times(out))
        rows.append(row)
df = pd.DataFrame(rows).sort_values(["scan", "value"])
df.to_csv(f"{csv_dir}/convergence.csv", index=False, float_format="%.6g")
print("wrote data/csv/convergence.csv")
