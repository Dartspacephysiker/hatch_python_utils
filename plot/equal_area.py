# 2019/10/08
# def equal_area():
from hatch_python_utils.earth import coordinates as hCoord
from hatch_python_utils.plot import colormaps as hCM

import matplotlib.pyplot as plt

import numpy as np

from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Polygon


def draw_screen_poly(lats, lons, m, color='red', alpha=0.4):
    x, y = m(lons, lats)
    xy = list(zip(x, y))
    poly = Polygon(xy, facecolor=color, alpha=alpha)
    plt.gca().add_patch(poly)
    return poly


def draw_on_map(plotstat, m,
                ea=None,
                vmin=None,
                vmax=None,
                tryAattepunkter=True,
                alphaval=0.8,
                cmap=hCM.parula,
                verbose=False):

    if ea is None:
        print("Defaulting to NH!")
        ea = hCoord.EqualAreaBins(hemi='north')

    plotmax, plotmin = np.max(plotstat[np.isfinite(plotstat)]), np.min(
        plotstat[np.isfinite(plotstat)])

    if vmin is None:
        vmin = plotmin
    if vmax is None:
        vmax = plotmax

    cmapvals = cmap(plt.Normalize(vmin, vmax)(plotstat))

    polylist = []
    for i, stat in enumerate(plotstat):

        if ~np.isfinite(stat):
            continue

        if verbose:
            print(i, stat)

        tmpis = ea.ea.mini[i], ea.ea.maxi[i]
        tmpms = ea.ea.minm[i], ea.ea.maxm[i]
        tmpms = np.rad2deg(np.unwrap(np.deg2rad(np.array(tmpms)*15)))

        midlat = np.mean(tmpis)
        midlon = np.mean(tmpms)

        lowmidlat = np.mean([tmpis[0], midlat])
        upmidlat = np.mean([tmpis[1], midlat])

        lowmidlon = np.mean([tmpms[0], midlon])
        upmidlon = np.mean([tmpms[1], midlon])

        if tryAattepunkter:
            lats = [tmpis[0], lowmidlat, midlat, upmidlat, tmpis[1],  # LH side
                    tmpis[1], tmpis[1], tmpis[1], tmpis[1],         # Upper
                    upmidlat, midlat, lowmidlat, tmpis[0],          # Right
                    tmpis[0], tmpis[0], tmpis[0]]                  # Lower

            lons = [tmpms[0], tmpms[0], tmpms[0], tmpms[0], tmpms[0],  # LH side
                    lowmidlon, midlon, upmidlon, tmpms[1],           # Upper
                    tmpms[1], tmpms[1], tmpms[1], tmpms[1],          # Right
                    upmidlon, midlon, lowmidlon]                    # Lower
        else:
            lats = [tmpis[0], midlat, tmpis[1], tmpis[1],
                    tmpis[1], midlat, tmpis[0], tmpis[0]]
            lons = [tmpms[0], tmpms[0], tmpms[0], midlon,
                    tmpms[1], tmpms[1], tmpms[1], midlon]

        polylist.append(draw_screen_poly(lats, lons, m,
                                         color=cmapvals[i, :],
                                         alpha=alphaval))
    return polylist, cmapvals
