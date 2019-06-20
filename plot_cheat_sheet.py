########################################
# AXISNESS
# AXI1  (20160611): Forandre fargene til sekundær axis
########################################
# LEGENDS
# LEG1  (20190604): Multiple axes, single legend
########################################

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

# AXI1  (20160611): Forandre fargene til sekundær axis


sekundFarge = 'orange'

x, y0, y1, y2 = fabdata()

fig, ax = plt.subplots(1, 1, figsize=(16, 10))

ax.grid()

junk = ax.set_xlabel("Time (s)")

l0 = ax.plot(x, y0, label="normal")
ax1 = ax.twinx()
l1 = ax1.plot(x, y1, label="gamma")

ax1.axes.spines['right'].set_color(sekundFarge)
ax1.yaxis.label.set_color(sekundFarge)
ax1.tick_params(axis='y', colors=sekundFarge)

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
