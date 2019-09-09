# 2019/07/03
from datetime import datetime,timedelta
import importlib
import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import interpolate
from scipy.signal import savgol_filter,argrelextrema
from scipy.interpolate import interp1d
import tkinter
mplBkgrnd = 'TkAgg'
mpl.use(mplBkgrnd)

import matplotlib.pyplot as plt

plt.ion()

import colorcet as cc           # Perceptually uniform spacing
cmap = cc.cm.colorwheel

mpl.rcParams.update({'text.color': 'k'})
mpl.rcParams.update({'axes.labelcolor': 'k'})
mpl.rcParams.update({'xtick.color': 'k'})
mpl.rcParams.update({'ytick.color': 'k'})
mpl.rcParams.update({'font.size': 15})
mpl.rcParams.update({'figure.figsize' : [10.0, 8.0]})

import hatch_python_utils
from hatch_python_utils import arrays as hArr
from hatch_python_utils.earth import coordinates as hCoord
junk = importlib.reload(hatch_python_utils)
junk = importlib.reload(hArr)
junk = importlib.reload(hCoord)

dataDir = '/SPENCEdata/Research/database/Rockets/CAPER2/caper-2-ny-alesund-data/spencers-work/'
plotDir = '/SPENCEdata/Research/sandbox_and_journals/plots/Rockets/CAPER2/'

#import CAPER2_tools
#junk = importlib.reload(CAPER2_tools)
#from CAPER2_tools import get_ELF_file, load_ELF_file, barplot_uniq_sRates, load_GPS_data,load_attitude_solution, despin_meas_using_Wallops_attitude, get_vXB_ENU, interp_measAB_to_common_tSeries, get_sun_spikes,load_magnetometer_file, get_magnetometer_file, get_interped_DCM_matrix
#
