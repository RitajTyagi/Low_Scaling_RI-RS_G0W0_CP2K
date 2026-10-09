"""
Main figure 4 -- G0W0+H@PBE on nanographene ribbons C14L H4L+8, L = 1 ... 512
(26 ... 9224 atoms).  Same layout as the poster (run-nanoribbon/poster/mkfigs.py):
x axis = number of strips L on a log2 scale with the number of atoms written
under every tick, and scaling exponents given in L (N = 18 L + 8, so for long
ribbons an exponent in L is the exponent in N).

Five setups are compared (folder name in data/nanographene/<L>/):
    tensor_avdz     CP2K tensor-based low-scaling GW,  aug-cc-pVDZ
    rirs_avdz       CP2K RI-RS,                        aug-cc-pVDZ
    tensor_mini     CP2K tensor-based low-scaling GW,  aug-SZV-MOLOPT-ae-mini
    rirs_mini       CP2K RI-RS,                        aug-SZV-MOLOPT-ae-mini
    rirs_opt_mini   CP2K RI-RS + r_AO = 6 A, r_W = 35 A, grid optimized to 30 %
                                                       aug-SZV-MOLOPT-ae-mini
The tensor-based curves are drawn as a thick pale "halo" underneath, so that
an RI-RS curve lying on top of them stays visible.

(a) HOMO-LUMO gap = "window" LUMO minus HOMO (lowest QP level of the virtual
    states within 2 eV of the KS LUMO, see make_tables.py).
(b) Cost of the full G0W0 step in core hours (wall time x cores of that run),
    SCF not included.  Dashed lines are power-law fits L^a to the first three
    and to the last three points of every setup; the legend gives
    "first three -> last three", and the last-three exponent is written at
    the end of each curve.  Setups with fewer than six points (aug-cc-pVDZ)
    get one fit, to their last three points.
(c) Cost of the individual steps of RI-RS + cutoffs + optimized grid, with
    "first three / last three" exponents, also for the total G0W0.  The SCF
    is NOT shown here; its comparison with G0W0 is Figure 5.  The SCF is NOT taken
    from the same run (those restarted from a converged wave function for
    large L) but from the reference runs listed in make_tables.py (column
    scf_ref_source): tensor_mini for L <= 64, rirs_mini/scf.out for L = 128
    and 512, rirs_mini/output.out for L = 256.

Machine: Noctua 2, 2 x AMD EPYC 7763 per node (128 cores), 8 MPI ranks x 16
threads per node, 1 ... 18 nodes.

Data: data/csv/ribbon_scaling.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/ribbon_scaling.csv")
for c in [c for c in d.columns if c.startswith("t_")]:
    d["ch_" + c[2:]] = d[c] * d.cores / 3600.0          # core hours
# SCF reference run (own cores), see make_tables.py
d["ch_scf"] = d.t_scf_ref * d.cores_scf_ref / 3600.0

L_ALL = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
ATOMS = {L: 18 * L + 8 for L in L_ALL}                    # C14L H4L+8

SETUPS = [  # (flavour, label, colour, marker, drawn as halo?)
    ("tensor_avdz",   "Tensor, aug-cc-pVDZ",          C_TENSOR,   's', True),
    ("rirs_avdz",     "RI-RS, aug-cc-pVDZ",           C_RIRS,     's', False),
    ("tensor_mini",   "Tensor, mini",                 C_TENSOR,   'o', True),
    ("rirs_mini",     "RI-RS, mini",                  C_RIRS,     'o', False),
    ("rirs_opt_mini", "RI-RS-Opt, mini",        C_RIRS_OPT, 'o', False),
]


def style(col, mk, halo):
    """Line/marker style.  The tensor reference is a thick pale line with
    open markers underneath; the colour is pre-blended with white instead of
    using alpha, because the figure is saved with a transparent background."""
    if halo:
        return dict(color=tint(col, 0.55), lw=8.0, marker=mk, mfc='none',
                    mec=col, mew=2.2, ms=18 if mk == 's' else 15, zorder=2)
    return dict(color=col, lw=2.6, marker=mk, mfc=col, mec='white', mew=1.2,
                ms=11 if mk == 's' else 9, zorder=4)


def slope(x, y):
    """Exponent a of a power-law fit y = c x^a (least squares in log-log)."""
    return np.polyfit(np.log10(x), np.log10(y), 1)[0]


def strips_axis(ax, xmax=640):
    """log2 axis in L at the bottom, number of atoms N = 18 L + 8 on top."""
    ax.set_xscale('log', base=2)
    ax.set_xticks(L_ALL)
    ax.set_xticklabels([f"$2^{{{int(np.log2(L))}}}$" for L in L_ALL], fontsize=22)
    ax.minorticks_off()
    ax.set_xlim(0.8, xmax)
    ax.set_xlabel("Number of strips $L$")
    # top axis: same tick positions, labelled with the number of atoms
    top = ax.secondary_xaxis('top')
    top.set_xticks(L_ALL)
    top.set_xticklabels([f"{ATOMS[L]:,}" for L in L_ALL], fontsize=18, rotation=45,
                        ha='left', rotation_mode='anchor')
    top.minorticks_off()
    top.tick_params(direction='in')
    top.set_xlabel("Number of atoms $N$", labelpad=10)
    ax.tick_params(axis='x', top=False)


def ribbon_L3():
    """Geometry of the L = 3 nanographene, built from the L = 4 xyz file.

    The ribbons are periodic along y with one C14H4 strip per 4.264 A; strip
    k of the L = 4 file is centred at y = k * 4.264 A.  We delete strip 2 and
    shift everything above it down by one period, which gives L = 3 with the
    same CH2-terminated ends.  Returns symbols, coordinates and the y range
    of the middle strip (the one that is repeated L - 2 times).
    """
    period = 4.264
    lines = open(f"{data_dir}/nanographene/4/struc.xyz").read().split("\n")[2:]
    rows = [l.split() for l in lines if l.strip()]
    sym = np.array([r[0] for r in rows])
    xyz = np.array([[float(v) for v in r[1:4]] for r in rows])
    lo, hi = 2 * period - period / 2, 2 * period + period / 2
    keep = ~((xyz[:, 1] > lo) & (xyz[:, 1] < hi))
    sym, xyz = sym[keep], xyz[keep].copy()
    xyz[xyz[:, 1] >= hi, 1] -= period
    return sym, xyz, (period - period / 2, period + period / 2)


def draw_ribbon(ax):
    """Ball-and-stick drawing of the L = 3 ribbon, long axis horizontal, with
    the repeated middle strip in brackets labelled L - 2."""
    sym, xyz, (y0, y1) = ribbon_L3()
    u, v = xyz[:, 1], xyz[:, 0]               # long axis horizontal
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    rcov = np.where(sym == "C", 0.76, 0.31)
    bonded = np.triu(d < 1.25 * (rcov[:, None] + rcov[None, :]), k=1)
    for i, j in zip(*np.where(bonded)):
        ax.plot([u[i], u[j]], [v[i], v[j]], color='#6e6e6e', lw=3.0,
                solid_capstyle='round', zorder=1)
    for el, col, size in (("C", '#3a3a3a', 150), ("H", 'white', 70)):
        m = sym == el
        ax.scatter(u[m], v[m], s=size, c=col, edgecolors='#3a3a3a',
                   linewidths=1.0, zorder=2)
    # brackets around the middle strip
    ylo, yhi = v.min() - 0.9, v.max() + 0.9
    for x, sgn in ((y0, 1), (y1, -1)):
        ax.plot([x + sgn * 0.45, x, x, x + sgn * 0.45], [ylo, ylo, yhi, yhi],
                color=INK, lw=2.4, zorder=3)
    ax.text(y1 + 0.2, ylo - 0.1, "$L-2$", ha='left', va='center',
            fontsize=24, color=INK)
    ax.set_aspect('equal')
    ax.axis('off')


fig, axes = plt.subplots(1, 3, figsize=(30, 10))
fig.subplots_adjust(wspace=0.27)

########## (a) gap ##########
ax = axes[0]
for fl, lab, col, mk, halo in SETUPS:
    s = d[d.flavour == fl].dropna(subset=["gap"]).sort_values("L")
    ax.plot(s.L, s.gap, label=lab, clip_on=False, **style(col, mk, halo))
g512 = d[(d.flavour == "rirs_opt_mini") & (d.L == 512)].gap.iloc[0]
ax.annotate(f"{g512:.3f} eV", xy=(512, g512), xytext=(-4, 22),
            textcoords='offset points', ha='right', fontsize=20, color=C_RIRS_OPT)
strips_axis(ax)
ax.set_ylabel("HOMO-LUMO gap (eV)")
ax.set_ylim(3.6, 9.4)
ax.legend(loc='upper right', fontsize=17, handlelength=2.4)
# structure of the nanographene in the empty middle-right part of (a)
ins = ax.inset_axes([0.40, 0.26, 0.58, 0.40])
draw_ribbon(ins)
panel_label(ax, "a")

########## (b) total G0W0 cost with first-3 / last-3 fits ##########
ax = axes[1]
handles = []
for fl, lab, col, mk, halo in SETUPS:
    s = d[d.flavour == fl].sort_values("L")
    x, y = s.L.values.astype(float), s.ch_total.values
    a_hi = slope(x[-3:], y[-3:])
    if len(x) >= 6:
        a_lo = slope(x[:3], y[:3])
        txt = f"{lab}\n$L^{{{a_lo:.2f}}} \\rightarrow L^{{{a_hi:.2f}}}$"
        # (points, extend left by, extend right by, write exponent at end?)
        fits = [(x[:3], y[:3], 1 / 2.2, 1.0, False), (x[-3:], y[-3:], 1.0, 1.35, True)]
    else:
        txt = f"{lab}\n$L^{{{a_hi:.2f}}}$"
        fits = [(x[-3:], y[-3:], 1.0, 1.9, True)]
    print(f"(b) {fl:14s} " + (f"first 3: L^{a_lo:.2f}  " if len(x) >= 6 else "") +
          f"last 3: L^{a_hi:.2f}")
    h, = ax.plot(x, y, label=txt, clip_on=False, **style(col, mk, halo))
    handles.append(h)
    for xr, yr, f0, f1, write in fits:
        a, b = np.polyfit(np.log10(xr), np.log10(yr), 1)
        xx = np.geomspace(xr[0] * f0, xr[-1] * f1, 30)
        ax.plot(xx, 10 ** b * xx ** a, ls=(0, (4, 3)), lw=1.8, color=col, zorder=3)
        if write:
            ax.annotate(f"$L^{{{a:.2f}}}$", xy=(xx[-1], 10 ** b * xx[-1] ** a),
                        xytext=(6, 2), textcoords='offset points', fontsize=21,
                        color=col, fontweight='bold', va='center')
strips_axis(ax, xmax=1700)       # room for the exponents at the end of the fits
ax.set_yscale('log')
ax.set_ylim(0.15, 1.5e6)
ax.set_ylabel(r"$G_0W_0$ computational cost (core hours)")
ax.legend(handles=handles, loc='lower right', fontsize=17, handlelength=2.4,
          labelspacing=0.6)
panel_label(ax, "b")

########## (c) steps of RI-RS + cutoffs, first-3 / last-3 ##########
ax = axes[2]
s = d[d.flavour == "rirs_opt_mini"].sort_values("L")
x = s.L.values.astype(float)
# the SCF is shown separately, in Figure 5 (fig5_scf_vs_gw.py)
for step, ls, mk, lw in (("total", '-', 'o', 3.4),):
    y = s["ch_" + step].values
    a_lo, a_hi = slope(x[:3], y[:3]), slope(x[-3:], y[-3:])
    print(f"(c) {step:8s} first 3: L^{a_lo:.2f}  last 3: L^{a_hi:.2f}")
    ax.plot(x, y, ls=ls, marker=mk, color=C_STEP[step], lw=lw, ms=9,
            mfc=C_STEP[step] if step == "total" else 'white', mec=C_STEP[step],
            mew=2.0, zorder=6, clip_on=False,
            label=f"{STEP_LABEL[step]}   $L^{{{a_lo:.2f}}}$ / $L^{{{a_hi:.2f}}}$")
for step in ("grid", "Z_lP", "chi", "W", "Sigma_x", "Sigma_c", "QP"):
    y = s["ch_" + step].values
    a_lo, a_hi = slope(x[:3], y[:3]), slope(x[-3:], y[-3:])
    print(f"(c) {step:8s} first 3: L^{a_lo:.2f}  last 3: L^{a_hi:.2f}")
    ax.plot(x, y, marker='o', color=C_STEP[step], mec='white', mew=1.0, ms=8,
            lw=2.2, zorder=4, clip_on=False,
            label=f"{STEP_LABEL[step]}   $L^{{{a_lo:.2f}}}$ / $L^{{{a_hi:.2f}}}$")
strips_axis(ax)
ax.set_yscale('log')
ax.set_ylim(1e-3, 3e6)
ax.set_ylabel("Computational cost (core hours)")
ax.legend(loc='upper left', ncol=2, fontsize=17, columnspacing=1.0,
          handlelength=1.8, labelspacing=0.45, bbox_to_anchor=(0.0, 0.93))
panel_label(ax, "c")

save(fig, "fig4_size_scaling")
