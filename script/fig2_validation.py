"""
Main figure 2 -- numerical precision of RI-RS G0W0 against FHI-aims.

(a) GW100 and (b) ACC24: mean signed (MSE), mean absolute (MAE) and maximum
absolute (MAX) deviation of the G0W0@PBE0 HOMO and LUMO from FHI-aims, for
both CP2K implementations: the tensor-based low-scaling GW (RI with the
truncated Coulomb metric) and RI-RS (this work).  The two CP2K bars are
almost identical; that is the point of the panel: the real-space
factorization adds no error on top of the RI-tC scheme it fits to.
(c) Nanographene ribbons L = 1, 2, 4, 8: G0W0@PBE HOMO and LUMO from FHI-aims
and CP2K RI-RS.

All with def2-TZVP (RI: def2-TZVP-RIFIT in CP2K, RI-V with the FHI-aims
auxiliary basis in FHI-aims).  LUMO = QP energy of the SCF LUMO orbital in
both codes (see make_tables.py).

FHI-aims reference: contour deformation for GW100 and ACC24; for the ribbons
analytic continuation (contour deformation for L = 1 differs by < 0.2 meV,
SI Table S1).

Data: data/csv/{gw100,acc24,ribbon_validation}.csv
"""

from plot_param import *

########## (a), (b): deviation statistics ##########


def stats(dev):
    """MSE, MAE, MAX of a deviation array (meV)."""
    return dev.mean(), np.abs(dev).mean(), np.abs(dev).max()


# full page width: drawn at 18 in, shown at ~7 in
fig = plt.figure(figsize=(18, 7.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.28)

panel_axes = []
for k, name in enumerate(("gw100", "acc24")):
    ax = fig.add_subplot(gs[0, k])
    d = pd.read_csv(f"{csv_dir}/{name}.csv").dropna(
        subset=["homo_aims", "lumo_aims", "homo_rirs", "lumo_rirs",
                "homo_tensor", "lumo_tensor"])

    labels = ["MSE", "MAE", "MAX"]
    width = 0.36
    x0 = 0.0
    xt, xl = [], []
    for state in ("homo", "lumo"):
        for j, (code, col, lab) in enumerate((("tensor", C_TENSOR, "CP2K Tensor"),
                                              ("rirs", C_RIRS, "CP2K RI-RS"))):
            dev = 1000 * (d[f"{state}_{code}"] - d[f"{state}_aims"])
            vals = stats(dev)
            xs = x0 + 1.3 * np.arange(3) + (j - 0.5) * width
            bars = ax.bar(xs, vals, width, color=col,
                          label=lab if state == "homo" else None)
            # number on top of every RI-RS bar only, to keep it readable
            # value centred over each tensor/RI-RS pair; the two bars agree
            # to < 1 meV, so one number per pair is enough
            if code == "rirs":
                for xx, v in zip(xs, vals):
                    ax.text(xx - width / 2, max(v, 0) + 1.5, f"{abs(v):.0f}" if abs(v) < 0.5 else f"{v:.0f}",
                            ha='center', va='bottom', fontsize=20)
        xt += list(x0 + 1.3 * np.arange(3))
        xl += labels
        ax.text(x0 + 1.3, -0.17, state.upper(), transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=24)
        x0 += 4.6

    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=18, rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='x', which='both', top=False, bottom=False)
    ax.set_ylim(-5, 100)
    ax.set_ylabel("Error w.r.t. CD-$G_0W_0$ FHI-aims (meV)", fontsize=21)
    if name == "acc24":        # (a) has no free space; one legend serves both
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.86), fontsize=19)
    ax.text(0.03, 0.96, {"gw100": "GW100", "acc24": "ACC24"}[name], transform=ax.transAxes,
            ha='left', va='top', fontsize=24, color=INK)
    panel_axes.append(ax)

########## (c): ribbons ##########

ax = fig.add_subplot(gs[0, 2])
d = pd.read_csv(f"{csv_dir}/ribbon_validation.csv")
d = d[d.basis == "tzvp"].sort_values("L")

for state in ("lumo", "homo"):
    ax.plot(d.L, d[f"{state}_aims_ac"], 's', mfc='none', mec=C_AIMS, ms=17,
            mew=2.0, label="FHI-aims" if state == "homo" else None)
    ax.plot(d.L, d[f"{state}_rirs"], 'o', color=C_RIRS, ms=8,
            label="CP2K RI-RS" if state == "homo" else None)
    ax.plot(d.L, d[f"{state}_rirs"], '-', color=C_RIRS, lw=1.5, zorder=0)
# HOMO / LUMO labels next to the curves (LUMO above, HOMO below the data)
ax.text(1.05, -1.65, "LUMO", ha="left", va="center",
        fontsize=20, color=INK_GREY)
ax.text(1.05, -5.55, "HOMO", ha="left", va="center",
        fontsize=20, color=INK_GREY)

dev_h = 1000 * np.abs(d.homo_rirs - d.homo_aims_ac)
dev_l = 1000 * np.abs(d.lumo_rirs - d.lumo_aims_ac)
print("(c) |RI-RS - FHI-aims| HOMO (meV):", np.round(dev_h.values, 1),
      " LUMO:", np.round(dev_l.values, 1))

ax.set_xscale('log', base=2)
ax.set_xticks([1, 2, 4, 8])
ax.set_xticklabels(["1", "2", "4", "8"])
ax.minorticks_off()
ax.set_xlim(0.8, 9.6)             # L = 8 close to the right edge
ax.set_ylim(-8.4, 2.6)
ax.set_xlabel(r"Number of strips $L$", fontsize=21)
ax.set_ylabel(r"$G_0W_0$@PBE energy (eV)", fontsize=21)   # same size as (a), (b)
ax.legend(loc='upper right', bbox_to_anchor=(1.0, 0.91), fontsize=16, handletextpad=0.3, borderaxespad=0.3, labelspacing=0.25)
ax.text(0.06, 0.96, "Nanoribbon", transform=ax.transAxes, ha='left', va='top',
        fontsize=24, color=INK)
panel_axes.append(ax)

# inset: absolute deviation from FHI-aims, the quantity that cannot be read
# off the eV scale of the main panel
ins = ax.inset_axes([0.45, 0.33, 0.52, 0.22])   # between the HOMO and LUMO curves
ins.plot(d.L, dev_h, 'o-', color=C_RIRS, ms=6, lw=1.6, label="HOMO")
ins.plot(d.L, dev_l, 'o--', color=C_RIRS, mfc='white', ms=6, lw=1.6, label="LUMO")
ins.set_xscale('log', base=2)
ins.set_xticks([1, 2, 4, 8])
ins.set_xticklabels(["1", "2", "4", "8"], fontsize=14)
ins.minorticks_off()
ins.set_xlim(0.8, 9.6)
ins.set_ylim(0, 13)
ins.set_yticks([0, 5, 10])
ins.tick_params(labelsize=14)
ins.set_ylabel(r"$|\Delta|$ (meV)", fontsize=16, labelpad=4)
ins.legend(loc="upper left", fontsize=13, handlelength=1.6, borderaxespad=0.2)

outer_labels(fig, panel_axes, "abc", fontsize=24)
save(fig, "fig2_validation")
