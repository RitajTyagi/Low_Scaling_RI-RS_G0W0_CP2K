"""
SI figure S1 -- convergence of the RI-RS G0W0+H@PBE HOMO-LUMO gap and of the
G0W0 wall time with the numerical parameters, for the L = 16 ribbon
(C224H72, 296 atoms, aug-SZV-MOLOPT-ae-mini, 32 MPI x 16 threads).

One parameter is scanned at a time.  r_c is plotted up to 13 A.  The grey
band marks +-10 meV around the gap at the production value, which we take as
the convergence criterion.  The others are NOT the same for all
scans (the scans were run while the production setup was being fixed); the
fixed values are written into the caption and summarized here:

  scan                      r_c (A)  EPS_FILTER  r_AO (A)  r_W (A)  grid
  (a) CUTOFF_RADIUS_RI       scan    1e-11        6        off      full
  (b) CUTOFF_RADIUS_RL_AO    14      1e-11        scan     off      full
  (c) CUTOFF_RADIUS_RL_W     12      1e-11        6        scan     full
  (d) EPS_FILTER             12      scan         6        30       full
  (e) N_POINT_PERCENTAGE     12      1e-10        6        35       scan

Production setup (Figures 3, 4): r_c = 12 A, EPS_FILTER = 1e-10, r_AO = 6 A,
r_W = 35 A, 30 % of the grid points.

Each panel: gap (left axis, filled circles) and G0W0 wall time (right axis,
open squares).  The production value is marked by a vertical band.

Data: data/csv/convergence.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/convergence.csv")

PANELS = [  # scan, x label, production value, log-x?
    # r_c is shown up to 13 A only (runs up to 16 A exist in the data)
    ("cutoff_radius_ri",    r"$r_\mathrm{c}$ of the RI metric (Å)", 12, False),
    ("cutoff_radius_rl_ao", r"$r_\mathrm{AO}$ (Å)", 6, False),
    ("cutoff_radius_rl_w",  r"$r_W$ (Å)", 35, False),
    ("eps_filter",          r"filter threshold $\epsilon$", 1e-10, True),
    ("grid_percentage",     "grid points kept (%)", 30, False),
]

C_GAP, C_TIME = C_RIRS_OPT, C_NEUTRAL   # gap of RI-RS-Opt runs; time is not a method
BAND = tint(C_RIRS_OPT, 0.80)
CONV = tint(INK, 0.88)       # +-10 meV band around the production gap

RC_MAX = 13.0                # largest r_c shown in (a)
TOL = 0.010                  # convergence criterion (eV), drawn as a band
GAP_YLIM = (3.85, 4.25)      # same gap axis (eV) in every panel

# two rows: (a) (b) (c) on top, (d) (e) centred below.  A 2 x 6 grid with
# every panel spanning two columns lets the lower row sit in the middle.
fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(2, 6, wspace=1.3, hspace=0.38)
SLOTS = [gs[0, 0:2], gs[0, 2:4], gs[0, 4:6], gs[1, 1:3], gs[1, 3:5]]
axes = [fig.add_subplot(sl) for sl in SLOTS]
LEFT, RIGHT = (0, 3), (2, 4)       # panels that carry the left / right y labels

for k, (scan, xlabel, prod, logx) in enumerate(PANELS):
    ax = axes[k]
    s = d[d.scan == scan].sort_values("value")
    if scan == "cutoff_radius_ri":
        s = s[s.value <= RC_MAX]
    g_prod = s[np.isclose(s.value, prod)].gap.iloc[0]
    ax.axhspan(g_prod - TOL, g_prod + TOL, color=CONV, zorder=0, lw=0)
    ax.plot(s.value, s.gap, '-o', color=C_GAP, ms=9)
    tx = ax.twinx()
    tx.plot(s.value, s.t_total, '--s', color=C_TIME, ms=9, mfc='none', mew=2)
    ax.axvspan(prod * (0.8 if logx else 1) - (0 if logx else 0.15 * np.ptp(s.value) / 4),
               prod * (1.25 if logx else 1) + (0 if logx else 0.15 * np.ptp(s.value) / 4),
               color=BAND, zorder=0, lw=0)
    if logx:
        ax.set_xscale('log')
        ax.invert_xaxis()          # tighter threshold to the right
    ax.set_xlabel(xlabel)
    # a tick at every scanned value (r_AO also at the missing 7 A)
    if logx:
        ax.set_xticks(s.value.values)
        ax.set_xticklabels([f"$10^{{{int(round(np.log10(v)))}}}$" for v in s.value],
                           fontsize=20)
        ax.minorticks_off()
    elif scan == "cutoff_radius_rl_ao":
        ax.set_xticks([5, 6, 7, 8])
    else:
        ax.set_xticks(s.value.values)
        ax.set_xticklabels([f"{v:g}" for v in s.value])
    if k in LEFT:
        ax.set_ylabel(r"$G_0W_0$+H gap (eV)", color=C_GAP)
    if k in RIGHT:
        tx.set_ylabel(r"$G_0W_0$ wall time (s)", color=C_TIME)
    ax.tick_params(axis='y', colors=C_GAP, right=False)
    tx.tick_params(axis='y', colors=C_TIME, direction='in')
    tx.set_ylim(0, 1.25 * s.t_total.max())
    lo, hi = s.gap.min(), s.gap.max()
    # one gap axis for all panels, ticks every 0.05 eV; the lower end is set
    # by (e) (10 % of the grid points: 3.888 eV)
    ax.set_ylim(*GAP_YLIM)
    ax.set_yticks(np.arange(GAP_YLIM[0], GAP_YLIM[1] + 1e-6, 0.05))
    panel_label(ax, "abcde"[k])
    # gap change relative to the tightest setting of the scan (meV)
    ref = s.gap.iloc[0] if (logx) else s.gap.iloc[-1]
    print(f"({'abcde'[k]}) {scan}: gap at production "
          f"{s[np.isclose(s.value, prod)].gap.iloc[0]:.3f} eV, "
          f"tightest {ref:.3f} eV, range {1000 * (hi - lo):.0f} meV")

save(fig, "figS1_convergence")
