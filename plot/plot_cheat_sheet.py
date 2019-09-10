########################################
# AXISNESS
# AXI1  (20190611): Forandre fargene til sekundær axis
# AXI2  (20190711): Gridspec!
# AXI3  (20190823): Sekundær akse med forskjellige tick-merker
########################################
# BOXPLOTS
# BP1   (20190904): Multi boxplots, change colors, etc.
########################################
# COLOR MAPS
# COL1  (20190910): Get colormap values for an array
########################################
# LEGENDS
# LEG1  (20190604): Multiple axes, single legend
########################################
# SPACING
# SPA1  (20190723): fig.suptitle with plt.tight_layout()
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
ax1 = ax.twinx()
l1 = ax1.plot(x, y1, label="gamma")

ax1.axes.spines['right'].set_color(c2)
ax1.yaxis.label.set_color(c2)
ax1.tick_params(axis='y', colors=c2)


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

ls = l0+l1+l2
llabs = [l.get_label() for l in ls]
leg = ax.legend(ls, llabs, loc=0)

########################################
# SPACING
########################################

# SPA1  (20190723): fig.suptitle with plt.tight_layout()
# Just do dis!
fig.tight_layout(rect=[0, 0.03, 1, 0.95])

########################################
# TEXT
########################################

# TEX1  (20190904): Text annotation example
fontsize = 35
normtextpos = (0.02, 0.86)       # Normalized axis coordinates
junkera = ax.text(*normtextpos, 'a',
                  transform=ax.transAxes,
                  fontsize=fontsize)
