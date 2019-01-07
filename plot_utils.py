# 2019/01/07
import matplotlib.dates as mdates
import pandas as pd


class plotHelp(object):

    def __init__(self, df):
        self.df = df

    def tickFormat_quant(self, x, quant, pos=None):
        dsdf = pd.Timestamp(mdates.num2date(x))
        return "{:.1f}".format(self.df.iloc[
            self.df.index.get_loc(dsdf.tz_localize(None),
                                  method='nearest')][quant])

    def tickFormat_mlat(self, x, pos=None):
        dsdf = pd.Timestamp(mdates.num2date(x))
        return "{:.1f}".format(self.df.iloc[
            self.df.index.get_loc(dsdf.tz_localize(None),
                                  method='nearest')].mlat)

    def tickFormat_mlt(self, x, pos=None):
        dsdf = pd.Timestamp(mdates.num2date(x))
        return "{:.1f}".format(self.df.iloc[
            self.df.index.get_loc(dsdf.tz_localize(None),
                                  method='nearest')].mlt)
