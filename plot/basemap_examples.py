# 1  (20190910): Make stereographic map
# 2  (20190910): Draw polygons (equal-area binning)


from hatch_python_utils.earth import coordinates as hCoord
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import tkinter
mplBkgrnd = 'TkAgg'
mpl.use(mplBkgrnd)


def meridianFmt(mer):
    if mer in [-270., -90., 0., 90., 180., 270., 360.]:
        return "{:.0f} MLT".format(mer/15.)
    else:
        return ""


def parallelFmt(par):
    return "{:.0f} MLat".format(par)


########################################
# 1  (20190910): Make stereographic map

figSize = (12, 12)
fig = plt.figure(figsize=figSize)
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
projection = 'cyl'
projection = 'npstere'


ea = hCoord.EqualAreaBins()

isSH = ea.hemi.lower() == 'south'
isNH = ea.hemi.lower() == 'north'

if isSH:
    if mirror_SH:
        parallels = np.arange(60., 86., 5.)
        projection = 'npstere'
        boundLat = 60
        latmax = np.max(parallels)
    else:
        parallels = np.arange(-85., -59., 5.)
        projection = 'spstere'
        boundLat = -60
        latmax = np.abs(np.min(parallels))
elif isNH:
    parallels = np.arange(60., 86., 5.)
    projection = 'npstere'
    boundLat = 60
    latmax = np.max(parallels)


if (projection == 'npstere') or (projection == 'spstere'):
    projOpts = dict(projection=projection,
                    boundinglat=boundLat,
                    lon_0=0,
                    resolution='l',
                    round=True)
# lon_0,lat_0 = 0, 90

# if projection == 'nsper':
#     satellite_height = 400*1000
m = Basemap(**projOpts)
# elif projection == 'cyl':

#     m = Basemap(projection=projection,
#                 # lon_0=lon_0,lat_0=lat_0,
#                 satellite_height=satellite_height,
#                 resolution='l')

latLabelDelta = 10

littleAlpha = 0.5
thickAlpha = 0.3

# meridians = np.arange(-180.,181.,15.)
meridians = np.arange(0, 360., 15.)
latLabelLon = 45

thinParDashes = [1, 1]
thickParDashes = [1, 0]
thickParWidth = 1.5

junkMer = m.drawmeridians(meridians, labels=[True, True, True, True],
                          fmt=meridianFmt,
                          dashes=thinParDashes,
                          latmax=latmax)
junkPar = m.drawparallels(parallels, labels=[True, True, True, True],
                          fmt=parallelFmt,
                          dashes=thinParDashes)


########################################
# 2  (20190910): Draw polygons (equal-area binning)
# See journal__20190815__lek_med_Swarm_crosstrack.ipynb
