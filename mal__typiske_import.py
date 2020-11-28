########################################
# Klassikere
from datetime import datetime, timedelta
import importlib
import numpy as np
import pandas as pd

########################################
# MPL stuff
import matplotlib as mpl
import tkinter
mplBkgrnd = 'TkAgg'
mpl.use(mplBkgrnd)

from hatch_python_utils.plot import colormaps as hCM
# import colorcet as cc           # Perceptually uniform spacing
# cmap = cc.cm.colorwheel
cmap = hCM.parula

########################################
# Scipy stuff
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.signal import savgol_filter, argrelextrema
from scipy.interpolate import interp1d

########################################
# HPU
import hatch_python_utils
from hatch_python_utils.earth import coordinates as hCoord
from hatch_python_utils import arrays as hArr
junk = importlib.reload(hatch_python_utils)
junk = importlib.reload(hArr)
junk = importlib.reload(hCoord)

########################################
#Dirs
datadir = '/SPENCEdata/Research/database/Rockets/CAPER2/caper-2-ny-alesund-data/spencers-work/'
plotdir = '/SPENCEdata/Research/sandbox_and_journals/plots/Rockets/CAPER2/'

########################################
# Hvis
# import sys
# SwarmDir = '/SPENCEdata/Research/sandbox_and_journals/journals/Swarm/'
# if SwarmDir not in sys.path:
#     sys.path.insert(0,SwarmDir)

mpl.rcParams.update({'text.color': 'k'})
mpl.rcParams.update({'axes.labelcolor': 'k'})
mpl.rcParams.update({'xtick.color': 'k'})
mpl.rcParams.update({'ytick.color': 'k'})
mpl.rcParams.update({'font.size': 15})
mpl.rcParams.update({'font.family': 'sans-serif'})
mpl.rcParams.update({'font.sans-serif': 'Arial'})
mpl.rcParams.update({'text.usetex': False})

mpl.rcParams.update({'figure.figsize': [10.0, 8.0]})
# mpl.rcParams.update({'savefig.directory': plotdir})

import matplotlib.pyplot as plt
plt.ion()

doGenTOC = False
if doGenTOC:
    import hatch_python_utils.nb_tools
    importlib.reload(hatch_python_utils.nb_tools)
    from hatch_python_utils.nb_tools import gen_toc
    gen_toc()
