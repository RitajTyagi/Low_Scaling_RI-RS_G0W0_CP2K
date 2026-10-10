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
      chi_PQ                 RI x RI
      Sigma^c_munu(i tau)    AO x AO  (first time point; Sigma^c(-i tau) has
                                       the same sparsity within 0.05 %)
    as printed by CP2K.  Dashed lines are power-law fits to the three largest
    systems; all five follow 1/L, i.e., the number of non-zero elements grows
    only linearly with L.
(The grid sizes, tabulated vs optimized, are only quoted in the text:
 35 and 10.5 grid points per AO for L = 512.)

Data: data/csv/ribbon_scaling.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/ribbon_scaling.csv")
opt = d[d.flavour == "rirs_opt_mini"].sort_values("L")
full = d[d.flavour == "rirs_mini"].sort_values("L")
L = opt.L.values.astype(float)

# drawn at 6 in width; in the paper it is shown at one column (3.4 in), so a
# font size of 15 here ends up at about 8.5 pt
fig, ax = plt.subplots(figsize=(6.0, 4.6))

MATS = [  # column, label, colour, marker
    ("occ_phi",   r"$\phi_\mu(\mathbf{r}_l)$", C_STEP['grid'],    'o'),
    ("occ_Z",     r"$Z_{lP}$",                  C_STEP['Z_lP'],    's'),
    # CP2K prints one common sparsity pattern for G_ll', chi_ll' and W_ll'
    ("occ_Gll",   r"$\chi_{ll'}$",              C_STEP['W'],       'D'),
    # the auxiliary-basis matrix, for comparison with the grid quantities
    ("occ_chiPQ", r"$\chi_{PQ}$",               C_STEP['chi'],     '^'),
    # the self-energy in the AO basis, Sigma^c(i tau) at the first time point
    # (Sigma^c(-i tau) is equally sparse)
    ("occ_Sigma", r"$\Sigma^c_{\mu\nu}(i\tau)$",  C_STEP['Sigma_c'], 'v'),
]
for col, lab, c, mk in MATS:
    y = opt[col].values
    ax.plot(L, y, marker=mk, color=c, ms=7, lw=2.0, mec='white', mew=0.8,
            label=lab, clip_on=False, zorder=3)
    # power-law fit to the three largest systems, drawn as a dashed line from
    # the third-last point to slightly beyond the last one, labelled L^a
    a_, b_ = np.polyfit(np.log(L[-3:]), np.log(y[-3:]), 1)
    xx = np.geomspace(L[-3], L[-1] * 1.3, 20)
    ax.plot(xx, np.exp(b_) * xx ** a_, '--', color=c, lw=1.6, zorder=2)
    print(f"{col:9s} " + "  ".join(f"{v:6.2f}" for v in y) + f"   last-3 exponent {a_:.2f}")

ax.set_xscale('log', base=2)
ax.set_xticks(L)
ax.set_xticklabels([f"$2^{{{int(np.log2(l))}}}$" for l in L], fontsize=16)
ax.minorticks_off()
ax.set_xlim(0.75, 900)   # a little room to see the tails of the dashed fits
ax.set_xlabel("Number of strips $L$", fontsize=19)
ax.set_yscale('log')
ax.set_ylim(0.2, 200)
# plain numbers instead of powers of ten on the y axis
ax.set_yticks([1, 10, 100])
ax.set_yticklabels(["1", "10", "100"], fontsize=16)
ax.set_ylabel("Non-zero elements (%)", fontsize=19)
# all five fits have exponent -1.0, so one label for all dashed lines
# just above the end of the chi_PQ curve
y_pq = opt["occ_chiPQ"].values
ax.text(np.sqrt(L[-2] * L[-1]), 2.6 * np.sqrt(y_pq[-2] * y_pq[-1]), r"$\propto 1/L$",
        ha="center", va="bottom", fontsize=22, color=INK)
ax.legend(loc='lower left', fontsize=15, handlelength=1.6, labelspacing=0.3)

save(fig, "fig3_sparsity_grid")
