"""
Main figure 5 -- how the total cost of a G0W0+H@PBE calculation on the
nanographenes is split between the PBE ground state (SCF) and the G0W0 step
(RI-RS-Opt), for L = 1 ... 512.

One stacked bar per L, normalised to 100 %: the lower part is the SCF, the
upper part G0W0.  The total cost (SCF + G0W0, core hours) is written above
every bar.  Normalised bars are used because the totals span four orders of
magnitude: on a linear axis the small systems vanish, and on a log axis the
upper segment of a stacked bar looks shorter than its share.

SCF cost: separate SCF calculations started from scratch, see make_tables.py
(column scf_ref_source): tensor_mini for L <= 64, rirs_mini/scf.out for
L = 128 and 512, rirs_mini/output.out for L = 256.
G0W0 cost: the `gw` timer of the RI-RS-Opt run (rirs_opt_mini/output.out).
Core hours = wall time x cores of the respective run.

Data: data/csv/ribbon_scaling.csv
"""

from plot_param import *

d = pd.read_csv(f"{csv_dir}/ribbon_scaling.csv")
d = d[d.flavour == "rirs_opt_mini"].sort_values("L")
scf = (d.t_scf_ref * d.cores_scf_ref / 3600.0).values      # core hours
gw = (d.t_total * d.cores / 3600.0).values
tot = scf + gw
L = d.L.values

C_SCF, C_GW = C_NEUTRAL, C_RIRS_OPT

fig, ax = plt.subplots(figsize=(7.5, 5.6))   # drawn at ~2x column width
x = np.arange(len(L))
ax.bar(x, 100 * scf / tot, 0.7, color=tint(C_SCF, 0.35), label="SCF (PBE)")
ax.bar(x, 100 * gw / tot, 0.7, bottom=100 * scf / tot, color=C_GW,
       label=r"$G_0W_0$ (RI-RS-Opt)")

# total core hours above every bar (what the numbers are is said in the caption)
for xi, t in zip(x, tot):
    txt = f"{t:,.0f}" if t >= 10 else f"{t:.1f}"
    ax.text(xi, 102, txt, ha='center', va='bottom', fontsize=17, rotation=90)

for Li, s_, g_ in zip(L, scf, gw):
    print(f"L = {Li:3d}: SCF {s_:9.1f}  G0W0 {g_:9.1f}  total {s_ + g_:9.1f} core h, "
          f"G0W0 share {100 * g_ / (s_ + g_):3.0f} %")

ax.set_xticks(x)
ax.set_xticklabels([f"$2^{{{int(np.log2(l))}}}$" for l in L], fontsize=19)
ax.tick_params(axis='x', top=False, bottom=False)
ax.set_xlabel("Number of strips $L$")
ax.set_ylabel("Share of total cost (%)")
ax.set_ylim(0, 128)
ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.set_xlim(-0.6, len(L) - 0.4)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=17)

save(fig, "fig5_scf_vs_gw")
