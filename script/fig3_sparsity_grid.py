"""
Main figure 3 -- sparsity of the grid quantities for the nanographenes C14L H4L+8, L = 1 ... 512, RI-RS-Opt setup
(aug-SZV-MOLOPT-ae-mini, r_AO = 6 A, r_W = 35 A, optimized grid with 30 % of
the tabulated points).

Occupation (percentage of non-zero matrix elements, after filtering with
    EPS_FILTER = 1e-10) of
      Phi_mu(r_l)            grid x AO
      Z_lP                   grid x RI
      chi_ll'                grid x grid  (CP2K prints one common sparsity
                                          pattern for G_ll', chi_ll', W_ll')
    as printed by CP2K.  The dotted line is 1/L: an occupation that falls as
    1/L means that the number of non-zero elements grows only linearly with L.
(The grid sizes, tabulated vs optimized, are only quoted in the text:
 35 and 10.5 grid points per AO for L = 512.)

Data: data/csv/ribbon_scaling.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/ribbon_scaling.csv")
opt = d[d.flavour == "rirs_opt_mini"].sort_values("L")
full = d[d.flavour == "rirs_mini"].sort_values("L")
L = opt.L.values.astype(float)

fig, ax = plt.subplots(figsize=(7.5, 5.6))


def strips_axis(ax):
    ax.set_xscale('log', base=2)
    ax.set_xticks(L)
    ax.set_xticklabels([f"$2^{{{int(np.log2(l))}}}$" for l in L], fontsize=17)
    ax.minorticks_off()
    ax.set_xlim(0.75, 700)
    ax.set_xlabel("Number of strips $L$", fontsize=21)


########## (a) occupation ##########
MATS = [  # column, label, colour, marker
    ("occ_phi",   r"$\phi_\mu(\mathbf{r}_l)$", C_STEP['grid'],    'o'),
    ("occ_Z",     r"$Z_{lP}$",                  C_STEP['Z_lP'],    's'),
    # CP2K prints one common sparsity pattern for G_ll', chi_ll' and W_ll'
    ("occ_Gll",   r"$\chi_{ll'}$",              C_STEP['Sigma_c'], 'D'),
]
for col, lab, c, mk in MATS:
    ax.plot(L, opt[col], marker=mk, color=c, ms=8, lw=2.2, mec='white', mew=1.0,
            label=lab, clip_on=False)
    print(f"(a) {col:9s} " + "  ".join(f"{v:6.2f}" for v in opt[col]))
# guide: occupation proportional to 1/L, anchored at the Z_lP value for L = 512
ll = np.array([16, 700.0])
ax.plot(ll, opt.occ_Z.values[-1] * 512 / ll, ':', color=INK, lw=1.6)
# (the dotted 1/L guide is explained in the caption)
strips_axis(ax)
ax.set_yscale('log')
ax.set_ylim(0.3, 200)
# plain numbers instead of powers of ten on the y axis
ax.set_yticks([1, 10, 100])
ax.set_yticklabels(["1", "10", "100"])
ax.set_ylabel("Non-zero elements (%)", fontsize=21)
ax.legend(loc='lower left', fontsize=15, handlelength=1.8)

save(fig, "fig3_sparsity_grid")
