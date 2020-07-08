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

# AXI1  (20190611): Forandre fargene til sekundær axis


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


# AXI2  (20190711): Gridspec!

fig3 = plt.figure(constrained_layout=True)
gs = fig3.add_gridspec(2, 3)
f3_ax00 = fig3.add_subplot(gs[0, :-1])
yunk = f3_ax00.set_title(uplegTitle)
f3_ax01 = fig3.add_subplot(gs[0, -1:], sharey=f3_ax00)
f3_ax10 = fig3.add_subplot(gs[1, :-1], sharex=f3_ax00, sharey=f3_ax00)
yunk = f3_ax10.set_title(downlegTitle)
f3_ax11 = fig3.add_subplot(gs[1, -1:], sharex=f3_ax01, sharey=f3_ax10)

f3_ax00.grid()
f3_ax01.grid()
f3_ax10.grid()
f3_ax11.grid()

# AXI3  (20190823): Sekundær akse med forskjellige tick-merker
ax2 = ax1.twiny()
ax1Xs = ax1.get_xticks()

ax2Xs = []
for X in ax1Xs:
    ax2Xs.append(X * 2)

ax2.set_xticks(ax1Xs)
ax2.set_xbound(ax1.get_xbound())
ax2.set_xticklabels(ax2Xs)

# AXI4  (20200122): Change tick labels

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


# AXI5  (20200407): Axis label formatting with StrMethodFormatter
data = [10.**np.arange(-3,1),np.arange(-3,1)]

fig,ax = plt.subplots(1,1)
_ = ax.plot(data[0],data[1])
ax.set_xscale('log')
ax.xaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:.2g}'))

# AXI6  (20200407): Convert log_10values on x axis to 10**(log_10 values)
@mpl.ticker.FuncFormatter
def little_formatter(x, pos):
    return "{:.2g}".format(10**x) if np.isclose(x %1,0) else ""

ax.xaxis.set_major_formatter(little_formatter)

# AXI7  (20200630): Nice scientific notation (base 10)
_ = ax.ticklabel_format(axis='y',style='sci',scilimits=(4,4))
_ = ax.yaxis.major.formatter._useMathText = True  # This chunk forces 10**x notation, as I recall

# AXI8  (20200701): Remove top and right spines
# Hide the right and top spines
_ = ax.spines['right'].set_visible(False)
_ = ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
_ = ax.yaxis.set_ticks_position('left')
_ = ax.xaxis.set_ticks_position('bottom')

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
