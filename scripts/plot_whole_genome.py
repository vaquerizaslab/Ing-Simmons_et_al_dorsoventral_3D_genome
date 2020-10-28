import fanc
import fanc.plotting as fancplot
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
from matplotlib.colors import LogNorm

hic_file = snakemake.input[0]
outfile = snakemake.output[0]

plt.rcParams['figure.figsize'] = [4, 4]
plt.rcParams['figure.dpi'] = 300
mpl.rcParams.update({'font.size': 12})

# add Nan colour to colourmap, here light grey
cmap = copy.copy(mpl.cm.get_cmap("germany"))
cmap.set_bad('0.9')

hic = fanc.load(hic_file)
m = hic.matrix()

# ticks at ends of chromosomes
tick_locations = [0] + [v[1] for v in hic.chromosome_bins.values()]
#chromosome name labels at middle of chromosomes
label_locations = [(v[0]+v[1])/2 for v in hic.chromosome_bins.values()]
labels = [k for k in hic.chromosome_bins.keys()]

fig = plt.figure()
ax = plt.axes()
# set appropriate vmax and vmin
plt.imshow(m, cmap='germany', norm = LogNorm(vmin=1e-5, vmax=1e-1))

# lots of axis formatting
ax.yaxis.set_major_locator(plt.FixedLocator(tick_locations))
ax.yaxis.set_major_formatter(plt.NullFormatter())
ax.yaxis.set_tick_params(which='minor', length=0)
ax.yaxis.set_ticks(label_locations, minor=True)
ax.yaxis.set_ticklabels(labels, minor=True)

ax.xaxis.set_major_locator(plt.FixedLocator(tick_locations))
ax.xaxis.set_major_formatter(plt.NullFormatter())
ax.xaxis.set_tick_params(which='minor', length=0)
ax.xaxis.set_ticks(label_locations, minor=True)
ax.xaxis.set_ticklabels(labels, minor=True)
cb = plt.colorbar()
cb.ax.minorticks_off()

fig.savefig(outfile, bbox_inches='tight')