"""
Shared matplotlib settings, paths and small file readers for every figure.

Every figure script starts with  `from plot_param import *`  so that all
panels share one look (serif/STIX fonts, inward ticks, same colours).

Folder layout of the repo:
    data/          raw CP2K / FHI-aims runs, and data/csv/ written by make_tables.py
    script/        this file, make_tables.py and one script per figure
    fig/           every figure, as transparent PDF + PNG
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os

########## Paths ##########
# work_dir is the folder holding this file, i.e. <repo>/script
work_dir = sys.path[0]
data_dir = os.path.join(work_dir, "..", "data")
csv_dir = os.path.join(data_dir, "csv")
fig_dir = os.path.join(work_dir, "..", "fig")

########## Plotting parameters ##########
plt.rcParams['axes.labelsize'] = 26
plt.rcParams['axes.titlesize'] = 24
plt.rcParams.update({'font.size': 22})
plt.rcParams['axes.linewidth'] = 1.6
plt.rcParams['xtick.major.pad'] = 6
plt.rcParams['ytick.major.pad'] = 6
plt.rcParams['xtick.major.size'] = 7
plt.rcParams['ytick.major.size'] = 7
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['xtick.major.width'] = 1.6
plt.rcParams['ytick.major.width'] = 1.6
plt.rcParams['axes.labelpad'] = 12
plt.rcParams['lines.linewidth'] = 2.5
plt.rcParams['lines.markersize'] = 9
plt.rcParams['legend.frameon'] = False
plt.rcParams['figure.figsize'] = [15, 10]
mpl.rcParams["font.family"] = "serif"
# "STIX" alone is not found on this machine; the installed names are these
mpl.rcParams["font.serif"] = ["STIXGeneral", "Times New Roman", "DejaVu Serif"]
mpl.rcParams["mathtext.fontset"] = "stix"
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
# keep text as text (not outlines) in the PDF, so it stays editable/searchable
mpl.rcParams['pdf.fonttype'] = 42

########## Colours ##########
# One colour per method / code, used the same way in EVERY figure of the paper
# and the SI.  Palette: the "npg" (Nature Publishing Group) palette as used in
# Nature-family journals (ggsci), warm red for the reference implementation
# we compare against, navy for RI-RS, teal for RI-RS-Opt, warm brown for the
# FHI-aims reference.
C_TENSOR = '#E64B35'      # CP2K tensor-based low-scaling GW (warm red)
C_RIRS = '#3C5488'        # CP2K RI-RS (this work, navy)
C_RIRS_OPT = '#00A087'    # CP2K RI-RS-Opt: cutoffs + optimised grid (teal)
C_AIMS = '#7E6148'        # FHI-aims reference (warm brown)
C_NEUTRAL = '#8C8C8C'     # quantities that are not a method (e.g. wall time)
INK = '#111111'           # panel labels, primary text
INK_GREY = '#55544f'      # secondary annotation text

# One colour per GW step, shared by the scaling figures (same npg family)
C_STEP = {
    'total':   '#111111',
    'scf':     '#7f7f7f',
    'grid':    '#8491B4',
    'Z_lP':    '#4DBBD5',
    'chi':     '#F39B7F',
    'W':       '#00A087',
    'Sigma_x': '#B09C85',
    'Sigma_c': '#E64B35',
    'QP':      '#3C5488',
}

# Labels of the GW steps; mathtext keeps them plain but correct
STEP_LABEL = {
    'total': r'total $G_0W_0$',
    'scf':   'SCF (PBE)',
    'grid':  'grid',
    'Z_lP':  r'$Z_{lP}$',
    'chi':   r'$\chi(i\tau)$',
    'W':     r'$W(i\tau)$',
    'Sigma_x': r'$\Sigma^\mathrm{x}$',
    'Sigma_c': r'$\Sigma^\mathrm{c}(i\tau)$',
    'QP':    'QP',
}


########## Small helpers ##########

def tint(color, frac):
    """Blend `color` `frac` of the way towards white and return it opaque.

    The figures are saved with a transparent background, and a
    semi-transparent patch drawn on a transparent background keeps its full
    saturation in the file.  Pre-blending with white here instead of passing
    alpha= keeps every fill looking the same whatever it is placed on.
    """
    rgb = np.array(mpl.colors.to_rgb(color))
    return tuple(rgb + frac * (1.0 - rgb))


def panel_label(ax, letter, extra="", fontsize=24, x=0.025, y=0.955):
    """Put "(a)" in the top-left corner INSIDE the axes.

    Panels carry only their letter; what they show belongs in the caption.
    `extra` adds a short identifying label (e.g. the data set name).
    """
    ax.text(x, y, f"({letter})" + ("  " + extra if extra else ""),
            transform=ax.transAxes, ha='left', va='top',
            fontsize=fontsize, color=INK)


def outer_labels(fig, axes, letters, fontsize=24, pad=0.004, one_row=False):
    """Put "(a)", "(b)", ... at the top-left corner of each panel's full
    extent, i.e. outside the plot box and including the axis labels and
    tick labels (and a top axis, if present)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    bbs = [inv.transform_bbox(ax.get_tightbbox(r)) for ax in axes]
    for ax, bb, letter in zip(axes, bbs, letters):
        # same height for all panels in one row (rows = equal plot-box top)
        top = ax.get_position().y1
        if one_row:          # caller says: all these panels form one row
            y = max(b.y1 for b in bbs)
        else:
            y = max(b.y1 for a, b in zip(axes, bbs) if abs(a.get_position().y1 - top) < 0.02)
        fig.text(bb.x0 - pad, y + pad, f"({letter})", ha='left', va='bottom',
                 fontsize=fontsize, color=INK)


def fit_exponent(x, y):
    """Slope of log(y) against log(x), i.e. the exponent a in y ~ x^a."""
    a, b = np.polyfit(np.log(x), np.log(y), 1)
    return a, np.exp(b)


def save(fig, name):
    """Save a figure as transparent PDF and PNG into <repo>/fig."""
    os.makedirs(fig_dir, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", {"dpi": 300})):
        fig.savefig(os.path.join(fig_dir, f"{name}.{ext}"),
                    bbox_inches='tight', pad_inches=0.05,
                    transparent=True, **kw)
    print(f"wrote {os.path.basename(fig_dir)}/{name}.pdf and .png")
