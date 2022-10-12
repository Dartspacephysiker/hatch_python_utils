########################################
# AXISNESS
# AXI1  (20190611): Forandre fargene til sekundær axis
# AXI2  (20190711): Gridspec!
# AXI3  (20190823): Sekundær akse med forskjellige tick-merker
# AXI4  (20200122): Change tick labels
# AXI5  (20200407): Axis label formatting with StrMethodFormatter
# AXI6  (20200407): Convert log_10values on x axis to 10**(log_10 values)
# AXI7  (20200630): Nice scientific notation (base 10)
# AXI8  (20200701): Remove top and right spines
# AXI9  (20201021): Three-column gridspec WITH colorbar, with y-axis labels only visible for the leftmost column
# AXI10 (20210112): Add a second x axis on the SAME SIDE as the original x axis, with custom labels
# AXI11 (20210223): Nice rotation of datetime labels on x axis (use fig.autofmt_xdate()!)
# AXI12 (20210316): Polarsubplot filled_cells example, with equal-area grid and colorbar
# AXI13 (20220721): Three-column gridspec with THREE colorbars at bottom

########################################
# LABELS
# LAB1  (20220106): When panels have shared axes, hide one panel's axis labels

########################################
# BOXPLOTS
# BP1   (20190904): Multi boxplots, change colors, etc.
########################################
# COLOR MAPS
# COL1  (20190910): Get colormap values for an array
########################################
# LEGENDS
# LEG1  (20190604): Multiple axes, single legend
# LEG2  (20191105): Gjør legend-alpha 1
########################################
# SPACING
# SPA1  (20190723): fig.suptitle with plt.tight_layout()
########################################
# SUBPLOTS
# SUB1  (20200222): Adjust subplots
# SUB2  (20200222): Preserve x scale for subplots of different size
# SUB3  (20200222): Preserve y scale for subplots of different size
########################################
# TEXT
# TEX1  (20190904): Text annotation example
########################################
# COLORBARS
# CLB1  (20201201): Scientific notation with colorbars
########################################
# COLORS
# COL1 (20210426): Lighten og darken a color

from matplotlib import pyplot as plt
import numpy as np


def fabdata():
    x = np.arange(0, 100, 1)
    y0 = np.random.normal(size=100)
    y1 = np.random.gamma(size=100)
    y2 = np.random.negative_binomial(5, 0.5, size=100)
    return x, y0, y1, y2

########################################
# AXISNESS
########################################

## AXI1  (20190611): Forandre fargene til sekundær axis


x, y0, y1, y2 = fabdata()

fig, ax = plt.subplots(1, 1, figsize=(16, 10))

ax.grid()

junk = ax.set_xlabel("Time (s)")

l0 = ax.plot(x, y0, label="normal")

c2 = 'orange'
ax2 = ax.twinx()
l1 = ax2.plot(x, y1, label="gamma")

ax2.axes.spines['right'].set_color(c2)
ax2.yaxis.label.set_color(c2)
ax2.tick_params(axis='y', colors=c2)


## AXI2  (20190711): Gridspec!

fig = plt.figure(constrained_layout=True)
gs = fig.add_gridspec(2, 3)
ax00 = fig.add_subplot(gs[0, :-1])
yunk = ax00.set_title(uplegTitle)
ax01 = fig.add_subplot(gs[0, -1:], sharey=ax00)
ax10 = fig.add_subplot(gs[1, :-1], sharex=ax00, sharey=ax00)
yunk = ax10.set_title(downlegTitle)
ax11 = fig.add_subplot(gs[1, -1:], sharex=ax01, sharey=ax10)

ax00.grid()
ax01.grid()
ax10.grid()
ax11.grid()

## AXI3  (20190823): Sekundær akse med forskjellige tick-merker
ax2 = ax1.twiny()
ax1Xs = ax1.get_xticks()

ax2Xs = []
for X in ax1Xs:
    ax2Xs.append(X * 2)

ax2.set_xticks(ax1Xs)
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(ax2Xs)

## AXI4  (20200122): Change tick labels

# V1
axnewticks = [0,10,20]
axnewticklabels = ['Its','WildN','Happy']

ax.set_xticks(axnewticks)
ax.set_xbound((axnewticks[0],axnewticks[1]))
ax.set_xticklabels(axnewticklabels)

# V2
axcurticks = ax.get_xticks()

axnewticklabels = []
for X in axcurticks:
    axnewticklabels.append(X * 2)
ax.set_xticklabels(axnewticklabels)


## AXI5  (20200407): Axis label formatting with StrMethodFormatter
data = [10.**np.arange(-3,1),np.arange(-3,1)]

fig,ax = plt.subplots(1,1)
_ = ax.plot(data[0],data[1])
ax.set_xscale('log')
ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.2g}'))

## AXI6  (20200407): Convert log_10values on x axis to 10**(log_10 values)
@mpl.ticker.FuncFormatter
def little_formatter(x, pos):
    return "{:.2g}".format(10**x) if np.isclose(x %1,0) else ""

ax.xaxis.set_major_formatter(little_formatter)

## AXI7  (20200630): Nice scientific notation (base 10)
_ = ax.ticklabel_format(axis='y',style='sci',scilimits=(4,4))
_ = ax.yaxis.major.formatter._useMathText = True  # This chunk forces 10**x notation, as I recall

## AXI8  (20200701): Remove top and right spines
# Hide the right and top spines
_ = ax.spines['right'].set_visible(False)
_ = ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
_ = ax.yaxis.set_ticks_position('left')
_ = ax.xaxis.set_ticks_position('bottom')

## AXI9  (20201021): Three-column gridspec WITH colorbar, with y-axis labels only visible for the leftmost column
# From journal__20201020__ELF_despun__L_and_R_polarization.ipynb

fig = plt.figure(figsize=(18,10),constrained_layout=True)
ncols=3
colwidth = 10
gs = fig.add_gridspec(1, ncols*colwidth+1)
ax00 = fig.add_subplot(gs[0, :colwidth])
ax01 = fig.add_subplot(gs[0, colwidth:colwidth*2], sharex=ax00,sharey=ax00)
ax02 = fig.add_subplot(gs[0, colwidth*2:colwidth*3], sharex=ax00, sharey=ax00)
axes = [ax00,ax01,ax02]
cax = fig.add_subplot(gs[0,colwidth*3:])

xstuff = np.arange(0,10)
ystuff = np.arange(0,20)
xstuff,ystuff = np.meshgrid(xstuff,ystuff)
datsize = (xstuff.shape[0]-1,xstuff.shape[1]-1)
exampdat = np.random.normal(size=datsize)


from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
levels = MaxNLocator(nbins=20).tick_values(exampdat.min(), exampdat.max())
# levels = np.arange(-35.5,35.5,1)

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('seismic')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

for iax,ax in enumerate(axes):
    showdat = np.random.normal(size=datsize)

    title = f'ELF L/R Polarization'

    im = ax.pcolormesh(xstuff, ystuff, showdat, cmap=cmap, norm=norm)
    # _ = ax.set_title(title)
    
    if iax == 0:
        _ = ax.set_ylabel('Frequency [Hz]')
    else:
        # _ = ax.yaxis.set_ticklabels([])        
        _ = [ylab.set_visible(False) for ylab in ax.yaxis.get_ticklabels()]

    _ = ax.set_xlabel('Time [s]')


## AXI10 (20210112): Add a second x axis on the SAME SIDE as the original x axis, with custom labels
def add_2nd_xaxis(ax,newlabels,posbelow=-0.25):
    ax2 = ax.twiny()
    # Add some extra space for the second axis at the bottom
    # fig.subplots_adjust(bottom=0.2)
    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")
    
    
    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", posbelow))
    # Turn on the frame for the twin axis, but then hide all 
    # but the bottom spine
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)
    # as @ali14 pointed out, for python3, use this
    # and for python2, use this
    #for sp in ax2.spines.itervalues():
    for sp in ax2.spines.values():
        sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xticks(ax.xaxis.get_ticklocs())
    ax2.xaxis.set_major_formatter(ax.xaxis.get_major_formatter())
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticklabels(newlabels)
    # ax2.set_xlabel(r"MLat [deg]")
    return ax2


## AXI11 (20210223): Nice rotation of datetime labels on x axis
# Just use fig.autofmt_xdate()!

## AXI12 (20210316): Polarsubplot filled_cells example, with equal-area grid and colorbar
from pysymmetry.visualization import grids

#Make the binning grid
LOWLAT=50
M0 = 4  #4
dr = 2  #2
mlat_, mlt_, mltres_ = grids.equal_area_grid(dr = dr, M0 = M0, N = int((90-LOWLAT)//dr))
ncells = len(mlt_)

#Construct the binned average dict    
grid = {'mlt':mlt_, 'mlat':mlat_, 'mltres':mltres_, 'cmlt':mlt_ + mltres_/2., 'cmlat':mlat_ + dr/2.}

grid['bias'] = np.sin(np.deg2rad((90-grid['mlat'])/(90-LOWLAT)*90))

vmin,vmax = grid['bias'].min(),grid['bias'].max()
vmin,vmax = 0,1
nlevels = 21
levels = np.linspace(vmin,vmax,nlevels)
cmap = 'Blues'
cmap = plt.get_cmap(cmap)

fig,ax = plt.subplots(1,1)
pax = Polarsubplot(ax)
pax.filled_cells(grid['mlat'],grid['mlt'],2,grid['mltres'],grid['bias'],levels=levels,cmap=cmap)

cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True),cmap=cmap))

## AXI13 (20220721): Three-column gridspec with THREE colorbars at bottom
# From pct_diff_ns_currents_for_varying_conditions__plot.py

fig = plt.figure(figsize=(8,4),constrained_layout=True)
ncols=3
nrows = 1
colwidth = 10
rowwidth = 20

padrow = 0
gs = fig.add_gridspec(nrows=nrows*rowwidth+1+padrow, ncols=ncols*colwidth)
ax00 = fig.add_subplot(gs[padrow:rowwidth+padrow,           :colwidth])
ax01 = fig.add_subplot(gs[padrow:rowwidth+padrow, colwidth  :colwidth*2], sharex=ax00,sharey=ax00)
ax02 = fig.add_subplot(gs[padrow:rowwidth+padrow, colwidth*2:colwidth*3], sharex=ax00, sharey=ax00)
cax00 = fig.add_subplot(gs[nrows*rowwidth+padrow:,          :colwidth])
cax01 = fig.add_subplot(gs[nrows*rowwidth+padrow:,colwidth  :colwidth*2])
cax02 = fig.add_subplot(gs[nrows*rowwidth+padrow:,colwidth*2:colwidth*3])

plt.setp(ax01.get_yticklabels(), visible=False)
plt.setp(ax02.get_yticklabels(), visible=False)

fig.suptitle(r"$\phi_{ca}$"+f": {ca_orientation}")

axes = [ax00,ax01,ax02]
caxes = [cax00,cax01,cax02]


########################################
# LABELS
########################################

# LAB1  (20220106): When panels have shared axes, hide one panel's axis labels

plt.setp(ax00.get_xticklabels(), visible=False)

########################################
# BOXPLOTS
########################################

# BP1   (20190904): Multi boxplots, change colors, etc.

# SEE journal__20190814__standalone_JGR_Figur_7.ipynb in /SPENCEdata/Research/sandbox_and_journals/journals/batch_jobs/Strangeway_2005/.
# Funksjoner update_bp osv. er nyttige

########################################
# Color maps
# COL1  (20190910): Get colormap values for an array
plotstat = np.random.uniform(0, 1000, 100)
plotmax, plotmin = np.max(plotstat[np.isfinite(plotstat)]), np.min(
    plotstat[np.isfinite(plotstat)])
cmapvals = cmap(plt.Normalize(plotmin, plotmax)(plotstat))

########################################
# LEGENDS
########################################

# LEG1 (20190604) : Multiple axes, single legend


x = np.arange(0, 100, 1)
y0 = np.random.normal(size=100)
y1 = np.random.gamma(size=100)
y2 = np.random.negative_binomial(5, 0.5, size=100)

fig5, ax = plt.subplots(1, 1, figsize=(16, 10), sharex=True)

ax.grid()

junk = ax.set_xlabel("Time (s)")

l0 = ax.plot(x, y0, label="normal")
l1 = ax.plot(x, y1, label="gamma")

ax2 = ax.twinx()
l2 = ax2.plot(x, y2, label="neg_binom(5,0.5)")

lines = l0+l1+l2
lines = [l0]+[l1]+[l2]                # IF SCATTER INSTEAD OF PLOT
llabs = [l.get_label() for l in lines]
leg = ax.legend(lines, llabs, loc=0)

# LEG2  (20191105): Gjør legend-alpha 1
leg = ax.get_legend()
for lh in leg.legendHandles:
    try:
        lh._legmarker.set_alpha(1)
    except:
        lh.set_alpha(1) # IF SCATTER INSTEAD OF PLOT

########################################
# SPACING
########################################

# SPA1  (20190723): fig.suptitle with plt.tight_layout()
# Just do dis!
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

########################################
# SUBPLOTS
########################################

# SUB1  (20200222): Adjust subplots

# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
# left = 0.125  # the left side of the subplots of the figure
# right = 0.9   # the right side of the subplots of the figure
# bottom = 0.1  # the bottom of the subplots of the figure
# top = 0.9     # the top of the subplots of the figure
# wspace = 0.2  # the amount of width reserved for space between subplots,
#               # expressed as a fraction of the average axis width
# hspace = 0.2  # the amount of height reserved for space between subplots,
#               # expressed as a fraction of the average axis height

# SUB2  (20200222): Preserve x scale for subplots of different size
short_ax_lim = 4.
f, (top_ax, bottom_ax) = plt.subplots(2,1,figsize=[9,5])
bottom_ax.set_xlim(0,5)
top_ax.set_xlim(0,short_ax_lim)

# first find the pixel where the bottom x-axis equals the desired limit.
short_ax_stop_pix = bottom_ax.transData.transform((short_ax_lim,0))  # (x,y) in pixels at data point (4,0)


pix_to_fig = f.transFigure.inverted() # Transformation to figure coordinates from display coordinates
short_ax_stop_fig = pix_to_fig.transform(short_ax_stop_pix)  # (x,y) in figure space (0,1)

top_orig_position_px = top_ax.bbox.corners()  # (ll, ul, lr, ur)
top_orig_anchor_px = top_orig_position_px[0]  # this is lower left corner.
top_orig_anchor_fig = pix_to_fig.transform(top_orig_anchor_px)  #convert to figure space
top_x_anchor, top_y_anchor = top_orig_anchor_fig

top_width = short_ax_stop_fig[0] - top_x_anchor
new_pos = (top_x_anchor, top_y_anchor,  top_width, top_ax.get_position().height)

top_ax.set_position(new_pos)

# SUB3  (20200222): Preserve y scale for subplots of different size
# PRESERVE Y SCALE
ref_ax_lim = (0,5.)
short_ax_lim = (0,4.)
f, (short_ax, ref_ax) = plt.subplots(2,1,figsize=[9,5])
ref_ax.set_ylim(ref_ax_lim)
short_ax.set_ylim(short_ax_lim)

# first find the pixel where the "tall" (reference) y-axis equals the desired limit.
# short_ax_stop_pix = ref_ax.transData.transform(short_ax_lim)  # (x,y) in pixels at data point
ref_ax_stop_pix = ref_ax.transData.transform(short_ax_lim)  # (x,y) in pixels at data point

pix_to_fig = f.transFigure.inverted() # Transformation to figure coordinates from display coordinates
# short_ax_stop_fig = pix_to_fig.transform(short_ax_stop_pix)  # (x,y) in figure space (0,1)
ref_ax_stop_fig = pix_to_fig.transform(ref_ax_stop_pix)  # (x,y) in figure space (0,1)

short_orig_position_px = short_ax.bbox.corners()  # (ll, ul, lr, ur)
short_orig_ll_px = short_orig_position_px[0]  # this is lower left corner.
short_orig_ll_fig = pix_to_fig.transform(short_orig_ll_px)  #convert to figure space
short_x_ll, short_y_ll = short_orig_ll_fig

# Use to regne ut short_height
ref_orig_position_px = ref_ax.bbox.corners()  # (ll, ul, lr, ur)
ref_orig_ll_px = ref_orig_position_px[0]  # this is lower left corner.
ref_orig_ll_fig = pix_to_fig.transform(ref_orig_ll_px)  #convert to figure space
ref_x_ll, ref_y_ll = ref_orig_ll_fig

short_height = ref_ax_stop_fig[1] - ref_y_ll

new_pos = (short_x_ll, short_y_ll,  short_ax.get_position().width, short_height)

short_ax.set_position(new_pos)
########################################
# TEXT
########################################

# TEX1  (20190904): Text annotation example
fontsize = 35
normtextpos = (0.02, 0.86)       # Normalized axis coordinates
junkera = ax.text(*normtextpos, 'a',
                  horizontalalignment='center', #  ['center' | 'right' | 'left' ]
                  verticalalignment='top',  # 	[ 'center' | 'top' | 'bottom' | 'baseline' ]
                  weight='bold', #[ 'normal' | 'bold' | 'heavy' | 'light' | 'ultrabold' | 'ultralight']
                  transform=ax.transAxes,
                  fontsize=fontsize)

########################################
# COLORBARS
########################################

# CLB1  (20201201): Scientific notation with colorbars
        cb.formatter.set_scientific(True)
        cb.formatter.set_powerlimits((6,6))
        cb.update_ticks()

########################################
# COLORS
########################################

# COL1 (20210426): Lighten og darken a color

# Frå https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib/49601444
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
