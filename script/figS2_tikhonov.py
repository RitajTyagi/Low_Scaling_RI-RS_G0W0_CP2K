"""
SI figure S2 -- effect of the Tikhonov parameter lambda on RI-RS G0W0@PBE0
HOMO and LUMO energies of the ACC24 set, def2-TZVP, against FHI-aims (contour
deformation).

lambda enters the fit of Z_lP as  D' = d D d + lambda * 1,  where d D d is
the Jacobi-scaled grid overlap matrix whose diagonal is 1.  lambda = 0 means
no regularization (Cholesky still succeeds for these molecules); 1e-8 is the
CP2K default and the value used everywhere else in the paper.

Bars: mean signed (MSE), mean absolute (MAE) and maximum absolute (MAX)
deviation in meV, for (a) HOMO and (b) LUMO.

Data: data/csv/acc24.csv (columns homo_tikh<k>, lumo_tikh<k>, lambda = 10^-k)
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/acc24.csv")
KS = [0, 2, 4, 6, 8, 10, 12]                       # folder names
LABELS = ["0"] + [f"$10^{{-{k}}}$" for k in KS[1:]]
# three shades of the RI-RS colour, light to full
METRICS = [("MSE", tint(C_RIRS, 0.65)), ("MAE", tint(C_RIRS, 0.30)), ("MAX", C_RIRS)]

# HOMO above, LUMO below, same y axis
fig, axes = plt.subplots(2, 1, figsize=(14, 13), sharey=True)
fig.subplots_adjust(hspace=0.35)

for ax, state, letter in ((axes[0], "homo", "a"), (axes[1], "lumo", "b")):
    rows = []
    for k in KS:
        dev = 1000 * (d[f"{state}_tikh{k}"] - d[f"{state}_aims"])
        rows.append((dev.mean(), dev.abs().mean(), dev.abs().max()))
    rows = np.array(rows)
    x = np.arange(len(KS))
    w = 0.26
    for j, (name, col) in enumerate(METRICS):
        ax.bar(x + (j - 1) * w, rows[:, j], w, color=col, label=name)
    # mark the default
    ax.axvspan(KS.index(8) - 0.5, KS.index(8) + 0.5, color=tint(C_RIRS_OPT, 0.8),
               zorder=0, lw=0)
    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(LABELS)
    ax.tick_params(axis='x', top=False, bottom=False)
    ax.set_xlabel(r"Tikhonov parameter $\lambda$")
    panel_label(ax, letter, state.upper())
    for k, r in zip(KS, rows):
        print(f"{state} lambda=1e-{k:<2d}: MSE {r[0]:6.1f} MAE {r[1]:6.1f} MAX {r[2]:6.1f} meV")

for ax in axes:
    ax.set_ylabel("Deviation from FHI-aims (meV)", fontsize=22)
axes[0].set_ylim(-3, 32)
axes[0].legend(loc='upper right', fontsize=20)

save(fig, "figS2_tikhonov")
