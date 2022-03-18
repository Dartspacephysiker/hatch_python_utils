########################################
# Klassikere
from datetime import datetime, timedelta
from hatch_python_utils.time_tools import mjd2000_to_datetime as mjd2dt
from hatch_python_utils.time_tools import datetime_to_mjd2000 as dt2mjd

import importlib
import numpy as np
import pandas as pd

########################################
#Dirs
datadir = '/SPENCEdata/Research/database/'
plotdir = '/SPENCEdata/Research/sandbox_and_journals/plots/'

########################################
# MPL stuff
import matplotlib as mpl
import tkinter
mplBkgrnd = 'TkAgg'
mpl.use(mplBkgrnd)

########################################
# MPL opts

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

# To get pandas to automatically handle plot labels that are dates
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
