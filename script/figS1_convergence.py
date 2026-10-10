"""
SI figure S1 -- convergence of the RI-RS G0W0+H@PBE HOMO-LUMO gap and of the
G0W0 wall time with the numerical parameters, for the L = 16 ribbon
(C224H72, 296 atoms, aug-SZV-MOLOPT-ae-mini, 32 MPI x 16 threads).

all scans start from the same production setup and vary one
parameter at a time (CP2K revision 0e87843579):

  r_c = 14 A, EPS_FILTER = 1e-10, r_AO = 6 A, r_W = 30 A,
  RS_AO_RATIO = 6, CUTOFF_ATOMIC_CLUSTER = 3 A

The grey band marks +-10 meV around the gap of this production setup (the
same run appears in every scan), which we take as the convergence criterion.
Points below the gap axis (r_AO = 4 A, cluster radius 1 A: 3.20 eV) get no
marker; the line runs to them and is cut at the lower edge of the axis.

Each panel: gap (left axis, filled circles) and G0W0 wall time (right axis,
open squares).  The production value is marked by a vertical band.

Data: data/csv/convergence.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/convergence.csv")

PANELS = [  # scan, x label, production value, log-x?
    ("cutoff_radius_ri",      r"$r_\mathrm{c}$ of the RI metric (Å)", 14, False),
    ("cutoff_radius_rl_ao",   r"$r_\mathrm{AO}$ (Å)", 6, False),
    ("cutoff_radius_rl_w",    r"$r_W$ (Å)", 30, False),
    ("eps_filter",            r"filter threshold $\epsilon$", 1e-10, True),
    ("rs_ao_ratio",           "grid points per AO", 6, False),
    ("cutoff_atomic_cluster", "grid-optimization cluster radius (Å)", 3, False),
]

C_GAP, C_TIME = C_RIRS_OPT, C_NEUTRAL   # gap of RI-RS-Opt runs; time is not a method
BAND = tint(C_RIRS_OPT, 0.80)
CONV = tint(INK, 0.88)       # +-10 meV band around the production gap

TOL = 0.010                  # convergence criterion (eV), drawn as a band
GAP_YLIM = (3.95, 4.20)      # same gap axis (eV) in every panel
# gap of the production setup (identical run in every scan)
G_REF = d[(d.scan == "rs_ao_ratio") & np.isclose(d.value, 6)].gap.iloc[0]

# two rows of three panels
fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(2, 3, wspace=0.45, hspace=0.38)
axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
LEFT, RIGHT = (0, 3), (2, 5)       # panels that carry the left / right y labels

for k, (scan, xlabel, prod, logx) in enumerate(PANELS):
    ax = axes[k]
    s = d[d.scan == scan].sort_values("value")
    ax.axhspan(G_REF - TOL, G_REF + TOL, color=CONV, zorder=0, lw=0)
    # points below the gap axis (r_AO = 4 A, cluster radius 1 A: 3.20 eV): the
    # line runs to them and is cut at the axis, but no marker is drawn
    inside = s.gap >= GAP_YLIM[0]
    ax.plot(s.value, s.gap, '-', color=C_GAP)
    ax.plot(s.value[inside], s.gap[inside], 'o', color=C_GAP, ms=9)
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
    # one gap axis for all panels, ticks every 0.05 eV
    ax.set_ylim(*GAP_YLIM)
    ax.set_yticks(np.arange(GAP_YLIM[0], GAP_YLIM[1] + 1e-6, 0.05))
    panel_label(ax, "abcdef"[k])
    # gap change relative to the tightest setting of the scan (meV)
    ref = s.gap.iloc[0] if (logx) else s.gap.iloc[-1]
    print(f"({'abcdef'[k]}) {scan}: gap at production "
          f"{s[np.isclose(s.value, prod, rtol=1e-6, atol=0)].gap.iloc[0]:.3f} eV, "
          f"tightest {ref:.3f} eV, range {1000 * (hi - lo):.0f} meV")

save(fig, "figS1_convergence")
