# 1  (20190910): Make stereographic map
# 2  (20190910): Draw polygons (equal-area binning)
# 3  (20191008): Apex grid lines
# 4  (20191016): Colorbar with equal-area kart

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
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Polygon

from hatch_python_utils.plot import equal_area as hEAPlot

# def draw_screen_poly( lats, lons, m,color='red',alpha=0.4):
#     x, y = m( lons, lats )
#     xy = list(zip(x,y))
#     poly = Polygon( xy, facecolor=color, alpha=alpha )
#     plt.gca().add_patch(poly)

figSize = (12,12)
fig = plt.figure(figsize=figSize)
ax = fig.add_axes([0.1,0.1,0.8,0.8])

ea = hCoord.EqualAreaBins()

projection = 'cyl'
projection = 'npstere'

isSH = ea.hemi.lower() == 'south'
isNH = ea.hemi.lower() == 'north'

isEqualArea = True
# eaParallels = np.array([48,60,72,84])
eaParallels = np.array([48,57,66,75,84])
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
# lon_0,lat_0 = 0, 90

# if projection == 'nsper':
#     satellite_height = 400*1000
m = Basemap(**projOpts)
# elif projection == 'cyl':
    
#     m = Basemap(projection=projection,
#                 # lon_0=lon_0,lat_0=lat_0,
#                 satellite_height=satellite_height,
#                 resolution='l')

# latLabelDelta = 10

littleAlpha = 0.5
thickAlpha = 0.3

# meridians = np.arange(-180.,181.,15.)
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

junk = m.drawcoastlines()

plotstat = np.random.normal(loc=0.0,scale=0.5,size=672)

polylist, cmapvals = hEAPlot.draw_on_map(plotstat, m,
                                         tryAattepunkter=True,
                                         alphaval=0.8,
                                         cmap=hCM.parula)

########################################
# 3  (20191008): Apex grid lines

from datetime import datetime
from apexpy import Apex

doApexGridLines = True

curTime = datetime(2015,1,1,hour=12)
apexRefHeight_km = 130

polcaplowlat = 70.
upperlat = 80.
deltalat = 10.
deltamlt = 3.

if isSH and not mirror_SH:
    lowlatlim__mlatglines = (-1.)*upperlat
    highlatlim__mlatglines = (-1.)*polcaplowlat

    lowlatlim__mlonglines = (-1.)*upperlat
    highlatlim__mlonglines = boundLat+10

else:
    lowlatlim__mlatglines = polcaplowlat
    highlatlim__mlatglines = upperlat

    lowlatlim__mlonglines = boundLat-10
    highlatlim__mlonglines = upperlat

if doApexGridLines:

    a = apexpy.Apex(curTime,refh=apexRefHeight_km)

    ##########
    # MLAT GRID LINES
    ##########
    mlonsForGrid = np.arange(-180.,181.,2)
    mlatsForGrid = np.arange(lowlatlim__mlatglines,highlatlim__mlatglines+1,deltalat)

    mlatGridLines = []
    for mlfg in mlatsForGrid:
        glatgrid,glongrid, gerrorgrid = a.apex2geo(np.array([mlfg]*len(mlonsForGrid)),
                                                   mlonsForGrid,
                                                   apexRefHeight_km)
        gridline = m.plot(glongrid,
                          glatgrid,
                          latlon=True,
                          linestyle='-',
                          color='black')
        
        mlatGridLines.append(gridline)

    ##########
    # MLON GRID LINES
    ##########

    # mlonsForGrid = np.arange(-180.,181.,30.)+lon_0
    mltsForGrid = np.arange(0,24,deltamlt)
    mlatsForGrid = np.arange(lowlatlim__mlonglines,highlatlim__mlonglines+1.,2.5)

    mlonsForGrid = []
    for mlt in mltsForGrid:
        mlonsForGrid.append(a.mlt2mlon(mlt,curTime))

    mlonGridLines = []
    for mlfg in mlonsForGrid:
        glatgrid,glongrid, gerrorgrid = a.apex2geo(mlatsForGrid,
                                                   np.array([mlfg]*len(mlatsForGrid)),
                                                   apexRefHeight_km)
        gridline = m.plot(glongrid,
                          glatgrid,
                          latlon=True,
                          linestyle=':',
                          color='black')
        
        mlonGridLines.append(gridline)
        # break

########################################
# 4  (20191016): Colorbar with equal-area kart

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
# norm(plotstat)
cax = fig.add_axes([0.87, 0.1, 0.03, 0.8])

cb1 = mpl.colorbar.ColorbarBase(cax,
                                cmap=hCM.parula,
                                norm=norm,
                                orientation='vertical')
yunk = cb1.set_label('Density')
