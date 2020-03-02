# 2019/10/08
# def equal_area():
from hatch_python_utils.earth import coordinates as hCoord
from hatch_python_utils.plot import colormaps as hCM

import matplotlib.pyplot as plt

import numpy as np

try:
    from mpl_toolkits.basemap import Basemap, cm
except:
    print("Can't load basemap in hatch_python_utils.plot.equal_area.py!")

from matplotlib.patches import Polygon


def PlotMap(hemi='north',
            fig=None,
            ax=None,
            isEqualArea=True,
            draw_coastlines=False,
            mirror_SH=False):

    if (fig is None) or (ax is None):
        figSize = (12, 12)
        fig = plt.figure(figsize=figSize)
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    littleAlpha = 0.5
    thickAlpha = 0.3

    ea = hCoord.EqualAreaBins(hemi=hemi)

    projection = 'cyl'
    projection = 'npstere'

    isSH = ea.hemi.lower() == 'south'
    isNH = ea.hemi.lower() == 'north'

    # eaParallels = np.array([48,60,72,84])
    eaParallels = np.array([48, 57, 66, 75, 84])
    eaboundLat = 57

    if isSH:
        if mirror_SH:
            if isEqualArea:
                parallels = eaParallels
                boundLat = eaboundLat
            else:
                parallels = np.arange(60., 86., 5.)
                boundLat = 60
            projection = 'npstere'
            latmax = np.max(parallels)
        else:
            if isEqualArea:
                parallels = np.flip(eaParallels*(-1.))
                boundLat = -1*eaboundLat
            else:
                parallels = np.arange(-85., -59., 5.)
                boundLat = -60
            projection = 'spstere'
            latmax = np.abs(np.min(parallels))
    elif isNH:
        if isEqualArea:
            parallels = eaParallels
            boundLat = eaboundLat
        else:
            parallels = np.arange(60., 86., 5.)
            boundLat = 60
        projection = 'npstere'
        latmax = np.max(parallels)

    if (projection == 'npstere') or (projection == 'spstere'):
        projOpts = dict(projection=projection,
                        boundinglat=boundLat,
                        lon_0=0,
                        resolution='l',
                        round=True)

    m = Basemap(ax=ax, **projOpts)

    if isEqualArea:
        meridians = np.arange(0, 360., 90.)
    else:
        meridians = np.arange(0, 360., 15.)

    latLabelLon = 45

    thinParDashes = [1, 1]
    thickParDashes = [1, 0]
    thickParWidth = 1.5

    def meridianFmt(mer):
        if mer in [-270., -90., 0., 90., 180., 270., 360.]:
            return "{:.0f} MLT".format(mer/15.)
        else:
            return ""

    def parallelFmt(par):
        return "{:.0f} MLat".format(par)

    junkMer = m.drawmeridians(meridians, labels=[True, True, True, True],
                              fmt=meridianFmt,
                              dashes=thinParDashes,
                              latmax=latmax)
    junkPar = m.drawparallels(parallels, labels=[True, True, True, True],
                              fmt=parallelFmt,
                              dashes=thinParDashes)

    if draw_coastlines:
        junk = m.drawcoastlines()

    return m


def draw_screen_poly(lats, lons, m, color='red', alpha=0.4):
    x, y = m(lons, lats)
    xy = list(zip(x, y))
    poly = Polygon(xy, facecolor=color, alpha=alpha)
    plt.gca().add_patch(poly)
    return poly


def make_2dgrid_object_from_bin_edges(xedges, yedges,
                                      flip_if_south_to_make_compatible_w_north=True):

    xmins, ymins = np.meshgrid(xedges[0:-1], yedges[0:-1])
    xmaxes, ymaxes = np.meshgrid(xedges[1:], yedges[1:])

    xmins = xmins.flatten()
    ymins = ymins.flatten()
    xmaxes = xmaxes.flatten()
    ymaxes = ymaxes.flatten()

    histothing = np.recarray(ymaxes.size,
                             dtype=np.dtype({'names': ['mini', 'maxi', 'minm', 'maxm'],
                                             'formats': [np.float32, np.float32, np.float32, np.float32]}))
    histothing.mini = xmins
    histothing.maxi = xmaxes

    histothing.minm = ymins
    histothing.maxm = ymaxes

    if flip_if_south_to_make_compatible_w_north:
        if np.where(histothing.mini < 0)[0].size == histothing.mini.size:
            print("This is southern grid! Flip to make compatible with northern grid ...")

            # tmpMinI = (-1.)*histothing.maxi
            # tmpMaxI = (-1.)*histothing.mini
            histothing.mini = np.flip(histothing.mini)
            histothing.maxi = np.flip(histothing.maxi)

    return histothing


def draw_on_map(plotstat, m,
                ea=None,
                vmin=None,
                vmax=None,
                tryAattepunkter=True,
                alphaval=0.8,
                cmap=hCM.parula,
                verbose=False):
    """
    'ea' can actually be whatever; it just needs to have members 'mini','maxi','minm','maxm' corresponding to each plotstat bin
    For example, use make_2dgrid_object_from_bin_edges(xedges,yedges) to make just such a thing
    """

    if ea is None:
        print("Defaulting to NH!")
        ea = hCoord.EqualAreaBins(hemi='north').ea
    elif hasattr(ea, 'ea'):
        ea = ea.ea

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

        tmpis = ea.mini[i], ea.maxi[i]
        tmpms = ea.minm[i], ea.maxm[i]
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
