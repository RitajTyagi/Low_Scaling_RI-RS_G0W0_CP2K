"""
Main figure 4 -- G0W0+H@PBE on nanographenes C14L H4L+8, L = 1 ... 512
(26 ... 9224 atoms), one-column layout:
  (a) structure (drawn with 11 strips; the L - 2 middle strips bracketed)
  (b) HOMO-LUMO gap, number of atoms N = 18 L + 8 on the top axis
  (c) cost of the full G0W0 step in core hours (wall time x cores, SCF not
      included); grey dashed lines: power-law fits L^a through the three
      largest systems of every setup, exponent written above the last point
  (d) cost of the individual steps of RI-RS-Opt and of the total G0W0, with
      the exponents of the three largest systems in the legend

Setups (folder name in data/nanographene/<L>/):
    tensor_avdz     CP2K tensor-based low-scaling GW,  aug-cc-pVDZ
    tensor_mini     CP2K tensor-based low-scaling GW,  aug-SZV-MOLOPT-ae-mini
    rirs_avdz       CP2K RI-RS,                        aug-cc-pVDZ
    rirs_mini       CP2K RI-RS,                        aug-SZV-MOLOPT-ae-mini
    rirs_opt_mini   CP2K RI-RS-Opt (cutoffs + optimized grid), mini
Colour = method (tensor red, RI-RS navy, RI-RS-Opt teal), marker = basis
(square aug-cc-pVDZ, circle mini); tensor curves dashed with large open
markers so that the RI-RS curves on top of them stay visible.

Machine: Noctua 2, 2 x AMD EPYC 7763 per node (128 cores), 8 MPI ranks x 16
threads per node, 1 ... 18 nodes.

Data: data/csv/ribbon_scaling.csv
"""

from plot_param import *
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter

d = pd.read_csv(f"{csv_dir}/ribbon_scaling.csv")
for c in [c for c in d.columns if c.startswith("t_")]:
    d["ch_" + c[2:]] = d[c] * d.cores / 3600.0          # core hours
d["ch_Sigma"] = d["ch_Sigma_x"] + d["ch_Sigma_c"]

L_ALL = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
ATOMS = {L: 18 * L + 8 for L in L_ALL}                    # C14L H4L+8

# (flavour, method, basis)
SETUPS = [
    ("tensor_avdz",   "tensor", "avdz"),
    ("tensor_mini",   "tensor", "mini"),
    ("rirs_avdz",     "rirs",   "avdz"),
    ("rirs_mini",     "rirs",   "mini"),
    ("rirs_opt_mini", "opt",    "mini"),
]
METHOD = {  # colour, line style, marker filled?, marker size, line width
    "tensor": dict(col=C_TENSOR,   ls=(0, (4, 2)), filled=False, ms=11, lw=1.6),
    "rirs":   dict(col=C_RIRS,     ls='-',         filled=True,  ms=6,  lw=1.6),
    "opt":    dict(col=C_RIRS_OPT, ls='-',         filled=True,  ms=7,  lw=2.4),
}
METHOD_LABEL = {"tensor": "Tensor", "rirs": "RI-RS", "opt": "RI-RS-Opt"}
MARKER = {"avdz": 's', "mini": 'o'}
BASIS_LABEL = {"avdz": "aug-cc-pVDZ", "mini": "aug-SZV-MOLOPT-ae-mini"}

# colours of the RI-RS-Opt steps: none of them is used for a method
STEPS = [  # (column, label, colour)
    ("grid",  "Grid optimization",                          '#8E6BBF'),
    ("Z_lP",  r"$Z_{lP}$",                                   '#4DBBD5'),
    ("chi",   r"$\chi(i\tau)$",                              '#E89F0C'),
    ("W",     r"$W(i\tau)$",                                 '#D0458F'),
    ("Sigma", r"$\Sigma^\mathrm{x}+\Sigma^\mathrm{c}(i\tau)$", '#7A8B2E'),
    ("QP",    "QP",                                          '#5A5A5A'),
]
C_FIT = '#9A9A9A'


def style(method, basis, clip=False):
    m = METHOD[method]
    mk = 'D' if method == "opt" else MARKER[basis]
    return dict(color=m["col"], ls=m["ls"], lw=m["lw"], marker=mk,
                ms=m["ms"] if mk != 'D' else 7,
                mfc=m["col"] if m["filled"] else 'white',
                mec=m["col"] if not m["filled"] else 'white',
                mew=1.8 if not m["filled"] else 0.6, clip_on=clip)


def slope(x, y):
    return np.polyfit(np.log10(x), np.log10(y), 1)


def fit_line(ax, x, y, text):
    """Grey dashed line of the power-law fit through the last three points,
    extended a little on both sides, and the exponent above the last point."""
    a, b = slope(x[-3:], y[-3:])
    xx = np.array([x[-3] / 1.35, x[-1] * 1.35])
    ax.plot(xx, 10 ** (a * np.log10(xx) + b), color=C_FIT, lw=1.3,
            ls=(0, (3, 2)), zorder=1)
    ax.annotate(text.format(a=a), xy=(x[-1], y[-1]), xytext=(0, 3),
                textcoords='offset points', fontsize=FG - 1, ha='center',
                va='bottom', color=INK, zorder=6)
    return a


def x_axis(ax, bottom=True):
    ax.set_xscale('log', base=2)
    ax.set_xticks(L_ALL)
    ax.minorticks_off()
    ax.set_xlim(0.8, 700)
    ax.set_axisbelow(True)
    ax.grid(True, which="major", axis="both", color="#d6d6d6", lw=0.8)
    ax.tick_params(axis='both', labelsize=FT)
    if bottom:
        ax.set_xticklabels([str(L) for L in L_ALL])
        ax.set_xlabel("Number of strips $L$", fontsize=FL)
    else:
        ax.tick_params(axis='x', labelbottom=False)


def log_y(ax):
    ax.set_yscale('log')
    ax.set_ylim(1e-3, 1e6)
    ax.yaxis.set_major_locator(LogLocator(base=10, numticks=20))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_locator(LogLocator(base=10, subs=(), numticks=20))


# drawn at 8 in width, shown at one column (3.37 in): FT = 19 -> ~8 pt
FT, FL, FG = 19, 21, 17
fig = plt.figure(figsize=(8, 19.0))
gs = fig.add_gridspec(5, 1, height_ratios=[0.30, 0.13, 0.8, 1.0, 1.0], hspace=0.08)
ax_s = fig.add_subplot(gs[0])
ax_b = fig.add_subplot(gs[2])
ax_c = fig.add_subplot(gs[3], sharex=ax_b)
ax_d = fig.add_subplot(gs[4], sharex=ax_b)

########## (a) structure ##########
def draw_ribbon_long(ax, n_strip=11):
    """Ball-and-stick drawing of the ribbon with n_strip strips, built from
    the L = 4 xyz file by repeating its second strip; the n_strip - 2 middle
    strips are bracketed and labelled L - 2."""
    period = 4.264
    lines = open(f"{data_dir}/nanographene/4/rirs_opt_mini/struc.xyz").read().split("\n")[2:]
    rows = [l.split() for l in lines if l.strip()]
    sym0 = np.array([r[0] for r in rows])
    xyz0 = np.array([[float(v) for v in r[1:4]] for r in rows])
    y = xyz0[:, 1]
    first, middle, last = y < period / 2, (y >= period / 2) & (y < 1.5 * period), y >= 2.5 * period
    parts = [(sym0[first], xyz0[first])]
    for j in range(n_strip - 2):
        parts.append((sym0[middle], xyz0[middle] + [0, j * period, 0]))
    parts.append((sym0[last], xyz0[last] + [0, (n_strip - 4) * period, 0]))
    sym = np.concatenate([a for a, _ in parts])
    xyz = np.concatenate([b for _, b in parts])
    u, v = xyz[:, 1], xyz[:, 0]               # long axis horizontal
    dd = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=-1)
    rcov = np.where(sym == "C", 0.76, 0.31)
    bonded = np.triu(dd < 1.25 * (rcov[:, None] + rcov[None, :]), k=1)
    for i, j in zip(*np.where(bonded)):
        ax.plot([u[i], u[j]], [v[i], v[j]], color='#6e6e6e', lw=1.7,
                solid_capstyle='round', zorder=1)
    for el, col, size in (("C", '#3a3a3a', 42), ("H", 'white', 22)):
        m = sym == el
        ax.scatter(u[m], v[m], s=size, c=col, edgecolors='#3a3a3a',
                   linewidths=0.7, zorder=2)
    y0, y1 = period / 2, (n_strip - 1.5) * period
    ylo, yhi = v.min() - 0.9, v.max() + 0.9
    for x, sgn in ((y0, 1), (y1, -1)):
        ax.plot([x + sgn * 0.45, x, x, x + sgn * 0.45], [ylo, ylo, yhi, yhi],
                color=INK, lw=2.2, zorder=3)
    ax.text(y1 + 0.5, ylo + 0.3, "$L-2$", ha='left', va='top',
            fontsize=FG + 3, color=INK)
    ax.set_aspect('equal')
    ax.axis('off')


draw_ribbon_long(ax_s)
ax_s.text(0.5, 1.0, "C$_{14L}$H$_{4L+8}$", transform=ax_s.transAxes,
          fontsize=FL, ha='center', va='bottom', color=INK)

########## (b) gap ##########
ax = ax_b
for fl, meth, bas in SETUPS:
    s_ = d[d.flavour == fl].dropna(subset=["gap"]).sort_values("L")
    ax.plot(s_.L, s_.gap, zorder=3 if meth == "tensor" else 4, **style(meth, bas))
g512 = d[(d.flavour == "rirs_opt_mini") & (d.L == 512)].gap.iloc[0]
ax.annotate(f"{g512:.3f} eV", xy=(512, g512), xytext=(8, 12),
            textcoords='offset points', ha='right', va='bottom',
            fontsize=FG, color=C_RIRS_OPT)
x_axis(ax, bottom=False)
ax.set_ylabel("HOMO-LUMO gap (eV)", fontsize=FL)
ax.set_ylim(3.5, 9.0)
t = ax.secondary_xaxis('top')
t.set_xticks(L_ALL)
t.set_xticklabels([f"{ATOMS[L]:,}" for L in L_ALL], fontsize=FT - 3,
                  rotation=45, ha='left', rotation_mode='anchor')
t.minorticks_off()
t.tick_params(direction='in')
t.set_xlabel("Number of atoms $N$", fontsize=FL, labelpad=8)
ax.tick_params(axis='x', top=False)
# legends: colour/line = method, marker = basis
h_meth = [Line2D([], [], **{**style(m, "mini"), "marker": None}, label=METHOD_LABEL[m])
          for m in ("tensor", "rirs", "opt")]
h_bas = [Line2D([], [], ls='none', marker=MARKER[b], ms=8, mfc=INK_GREY,
                mec='white', label=BASIS_LABEL[b]) for b in ("avdz", "mini")]
h_bas.append(Line2D([], [], ls='none', marker='D', ms=7, mfc=C_RIRS_OPT,
                    mec='white', label="mini, opt. grid"))
h_bas.append(Line2D([], [], ls='none', marker='o', ms=10, mfc='white',
                    mec=INK_GREY, mew=1.8, label="open: tensor"))
ax.legend(handles=h_meth + h_bas[:2], loc='upper right', fontsize=FG,
          frameon=False, handlelength=2.4,
          labelspacing=0.3, borderpad=0.4)

########## (c) total G0W0 cost ##########
ax = ax_c
for fl, meth, bas in SETUPS:
    s_ = d[d.flavour == fl].sort_values("L")
    x, y = s_.L.values.astype(float), s_.ch_total.values
    ax.plot(x, y, zorder=3 if meth == "tensor" else 4, **style(meth, bas))
    a = fit_line(ax, x, y, "$L^{{{a:.1f}}}$")
    print(f"(c) {fl:14s} last 3: L^{a:.2f}")
x_axis(ax, bottom=False)
log_y(ax)
ax.set_ylabel(r"$G_0W_0$ cost (core hours)", fontsize=FL)

########## (d) steps of RI-RS-Opt ##########
ax = ax_d
s_ = d[d.flavour == "rirs_opt_mini"].sort_values("L")
x = s_.L.values.astype(float)
ax.plot(x, s_.ch_total.values, zorder=5, label="RI-RS-Opt total",
        **style("opt", "mini"))
a = fit_line(ax, x, s_.ch_total.values, "$L^{{{a:.1f}}}$")
print(f"(d) total    last 3: L^{a:.2f}")
for step, lab, col in STEPS:
    y = s_["ch_" + step].values
    ax.plot(x, y, color=col, lw=1.8, marker='o', ms=6, mec='white', mew=0.6,
            clip_on=False, zorder=4, label=lab)
    a, b = slope(x[-3:], y[-3:])
    print(f"(d) {step:8s} last 3: L^{a:.2f}")
x_axis(ax, bottom=True)
log_y(ax)
ax.set_ylabel("Cost of the steps (core hours)", fontsize=FL)
# exponents of the steps in the legend (six fit lines would hide the data)
hs, ls_ = ax.get_legend_handles_labels()
ls_ = [l if l == "RI-RS-Opt total" else
       f"{l}  ($L^{{{slope(x[-3:], s_['ch_' + st].values[-3:])[0]:.1f}}}$)"
       for l, st in zip(ls_, ["total"] + [s[0] for s in STEPS])]
ax.legend(hs, ls_, loc='upper left', fontsize=FG - 2, frameon=False,
          handlelength=1.5,
          labelspacing=0.25, borderpad=0.3, ncol=2, columnspacing=0.8,
          handletextpad=0.4)

# (a) at the same x as (b)-(d), level with the top of the structure panel
fig.canvas.draw()
_pb = ax_b.get_position()
_top = fig.transFigure.inverted().transform_bbox(ax_s.get_tightbbox(fig.canvas.get_renderer())).y1
fig.text(_pb.x0 - 0.17 * _pb.width, _top, "(a)", ha='left', va='top',
         fontsize=FL + 2, color=INK)
for ax, l in ((ax_b, "b"), (ax_c, "c"), (ax_d, "d")):
    ax.text(-0.17, 1.0, f"({l})", transform=ax.transAxes, ha='left',
            va='top', fontsize=FL + 2, color=INK)
save(fig, "fig4_size_scaling")
