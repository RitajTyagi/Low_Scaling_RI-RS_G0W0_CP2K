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

GW100 molecule 37 (CBr4) is dropped: the FHI-aims run did not finish.

Data: data/csv/{gw100,acc24,ribbon_validation}.csv
"""

from plot_param import *

########## (a), (b): deviation statistics ##########


def stats(dev):
    """MSE, MAE, MAX of a deviation array (meV)."""
    return dev.mean(), np.abs(dev).mean(), np.abs(dev).max()


fig = plt.figure(figsize=(24, 8.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.85], wspace=0.28)

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
            xs = x0 + np.arange(3) + (j - 0.5) * width
            bars = ax.bar(xs, vals, width, color=col,
                          label=lab if state == "homo" else None)
            # number on top of every RI-RS bar only, to keep it readable
            # value centred over each tensor/RI-RS pair; the two bars agree
            # to < 1 meV, so one number per pair is enough
            if code == "rirs":
                for xx, v in zip(xs, vals):
                    ax.text(xx - width / 2, max(v, 0) + 1.5, f"{abs(v):.0f}" if abs(v) < 0.5 else f"{v:.0f}",
                            ha='center', va='bottom', fontsize=20)
        xt += list(x0 + np.arange(3))
        xl += labels
        ax.text(x0 + 1, -0.12, state.upper(), transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=24)
        x0 += 3.8

    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xticks(xt)
    ax.set_xticklabels(xl, fontsize=19)
    ax.tick_params(axis='x', which='both', top=False, bottom=False)
    ax.set_ylim(-5, 100)
    ax.set_ylabel("Deviation from CD-$G_0W_0$ FHI-aims (meV)", fontsize=21)
    if name == "acc24":        # (a) has no free space; one legend serves both
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.86), fontsize=19)
    panel_label(ax, "abc"[k], {"gw100": "GW100", "acc24": "ACC24"}[name])

########## (c): ribbons ##########

ax = fig.add_subplot(gs[0, 2])
d = pd.read_csv(f"{csv_dir}/ribbon_validation.csv")
d = d[d.basis == "tzvp"].sort_values("L")

for state, lab in (("lumo", "LUMO"), ("homo", "HOMO")):
    ax.plot(d.L, d[f"{state}_aims_ac"], 's', mfc='none', mec=C_AIMS, ms=17,
            mew=2.0, label="FHI-aims" if state == "homo" else None)
    ax.plot(d.L, d[f"{state}_rirs"], 'o', color=C_RIRS, ms=8,
            label="CP2K RI-RS" if state == "homo" else None)
    ax.plot(d.L, d[f"{state}_rirs"], '-', color=C_RIRS, lw=1.5, zorder=0)
    ax.text(10.5, d[f"{state}_rirs"].iloc[-1], lab, va='center', ha='left',
            fontsize=20, color=INK_GREY)

maxdev = 1000 * np.nanmax(np.abs(np.r_[d.homo_rirs - d.homo_aims_ac,
                                       d.lumo_rirs - d.lumo_aims_ac]))
print(f"(c) max |RI-RS - FHI-aims| over HOMO and LUMO, L = 1-8: {maxdev:.1f} meV")

ax.set_xscale('log', base=2)
ax.set_xticks([1, 2, 4, 8])
ax.set_xticklabels(["1", "2", "4", "8"])
ax.minorticks_off()
ax.set_xlim(0.8, 22)
ax.set_ylim(-8.4, 2.6)      # headroom so '(c)' does not touch the L = 1 LUMO
ax.set_xlabel(r"Number of strips $L$", fontsize=21)
ax.set_ylabel(r"$G_0W_0$@PBE energy (eV)", fontsize=21)   # same size as (a), (b)
ax.legend(loc='center right', bbox_to_anchor=(1.0, 0.45), fontsize=19)
panel_label(ax, "c", "Nanoribbon")

save(fig, "fig2_validation")
