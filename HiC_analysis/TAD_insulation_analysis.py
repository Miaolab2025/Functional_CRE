## The insulation score algorithm in cooltools (v0.4.1) was used to determine the location of the TAD boundary for each sample, and the location and number of TADs. 


# conda activate PYTHON3.9
# python
# example: 
    # rs1652965[14:34,618,269 (GRCh37)] PHD3[34,393,433 - 34,420,280] range[34Mb-35Mb]
    # CRC-N_hic.mcool  CRC-A_hic.mcool  CRC-C_hic.mcool

# Import python package for working with cooler files and tools for analysis
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import cooler
import bioframe
from matplotlib.colors import LogNorm
from matplotlib.ticker import EngFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cooltools import insulation
from cooltools.lib.numutils import adaptive_coarsegrain, interp_nan

# Set parameters
data_dir = '~/user/HiC_data/'
samples = ['CRC-N_hic', 'CRC-A_hic', 'CRC-C_hic']
chrom = 'chr14'
region_start = 34_000_000
resolution = 10_000
window = 20 * resolution  # insulation window size
region = (chrom, region_start, region_start + 5 * window)  # 50 bins region

# Set plotting style
plt.rcParams['font.size'] = 12
sns.set_style("white")
cmap_hic = 'fall'
cmap_ins = plt.cm.plasma
norm = LogNorm(vmin=0.0001, vmax=0.01)
bp_formatter = EngFormatter('b')  # formats axis ticks (like 34M)

# Tick formatter helper
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x', rotation=45)

# Custom 45° Hi-C triangle plot
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    n = matrix_c.shape[0]
    start_pos_vector = [start + resolution * i for i in range(n + 1)]
    import itertools
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(np.array([(i[1], i[0])
                                for i in itertools.product(start_pos_vector[::-1],
                                                           start_pos_vector)]), t)
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im

# Main plotting loop for each sample
for sample in samples:
    print(f"Processing sample: {sample}")
    
    # Load .cool file
    clr_path = os.path.join(data_dir, f"{sample}.mcool::resolutions/{resolution}")
    clr = cooler.Cooler(clr_path)

    # Compute insulation score
    insulation_table = insulation(clr, [window], verbose=True)
    insul_region = bioframe.select(insulation_table, region)

    # Fetch balanced matrix and apply adaptive coarse graining
    raw_matrix = clr.matrix(balance=True).fetch(region)
    raw_matrix_unbalanced = clr.matrix(balance=False).fetch(region)
    cg = adaptive_coarsegrain(raw_matrix, raw_matrix_unbalanced, cutoff=3, max_levels=8)
    cgi = interp_nan(cg)  # interpolate NaNs to clean up

    # Plot Hi-C matrix and insulation profile
    fig, ax = plt.subplots(figsize=(18, 6))
    im = pcolormesh_45deg(ax, cgi, start=region[1], resolution=resolution, norm=norm, cmap=cmap_hic)
    ax.set_aspect(0.5)
    ax.set_ylim(0, 6 * window)
    format_ticks(ax, rotate=False)

    # Add colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="8%", pad=0.12, aspect=6)
    plt.colorbar(im, cax=cax)

    # Plot insulation score below triangle
    ins_ax = divider.append_axes("bottom", size="50%", pad=0., sharex=ax)
    ins_ax.set_prop_cycle(plt.cycler("color", cmap_ins(np.linspace(0, 1, 5))))
    ins_ax.plot(insul_region[['start', 'end']].mean(axis=1),
                insul_region[f'log2_insulation_score_{window}'],
                label=f'{sample} - {window} bp')
    ins_ax.legend(bbox_to_anchor=(0., -1), loc='lower left', ncol=4)
    format_ticks(ins_ax, y=False, rotate=False)
    ax.set_xlim(region[1], region[2])

    # Save to PDF
    out_pdf = f"./{sample}_chr14_TAD.pdf"
    plt.savefig(out_pdf, dpi=600, bbox_inches='tight')
    plt.close()
    print(f"Saved plot to {out_pdf}")
