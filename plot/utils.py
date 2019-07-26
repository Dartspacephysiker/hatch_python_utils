# 2019/01/07
import matplotlib.dates as mdates
import pandas as pd

# NOTE; if using tickFormat_ funcs as a lambda argument, you MUST include pos
# as an argument, e.g.,
# mltLabeler = lambda x, pos: hpu.tickFormat_quant(x, myDataFrame, 'mlt')
# mlatLabeler = lambda x, pos: hpu.tickFormat_quant(x, myDataFrame, 'mlat')


def tickFormat_quant_time(x, df, quant, pos=None):
    dsdf = pd.Timestamp(mdates.num2date(x))
    return "{:.1f}".format(df.iloc[
        df.index.get_loc(dsdf.tz_localize(None),
                         method='nearest')][quant])


def tickFormat_quant(x, df, quant, pos=None,
                     isTime=False,
                     timeFmtStr='%H:%M'):
    if isTime:
        return df.iloc[df.index.get_loc(x, method='nearest')][quant].strftime(timeFmtStr)
    else:
        return "{:.1f}".format(df.iloc[
            df.index.get_loc(x, method='nearest')][quant])


def tickFormat_mlt(x, df, pos=None):
    dsdf = pd.Timestamp(mdates.num2date(x))
    return "{:.1f}".format(df.iloc[
        df.index.get_loc(dsdf.tz_localize(None),
                         method='nearest')].mlt)


def tickFormat_mlat(x, df, pos=None):
    dsdf = pd.Timestamp(mdates.num2date(x))
    return "{:.1f}".format(df.iloc[
        df.index.get_loc(dsdf.tz_localize(None),
                         method='nearest')].mlat)

# class plotHelp(object):

#     def __init__(self, df):
#         self.df = df

#     def tickFormat_quant(self, x, quant, pos=None):
#         dsdf = pd.Timestamp(mdates.num2date(x))
#         # return "{:.1f}".format(self.df.iloc[
#         #     self.df.index.get_loc(dsdf.tz_localize(None),
#         #                           method='nearest')][quant])
#         that = self.df[quant]
#         return "{:.1f}".format(that.iloc[
#             that.index.get_loc(dsdf.tz_localize(None),
#                                method='nearest')])

#     def tickFormat_mlat(self, x, pos=None):
#         dsdf = pd.Timestamp(mdates.num2date(x))
#         # return "{:.1f}".format(self.df.iloc[
#         #     self.df.index.get_loc(dsdf.tz_localize(None),
#         #                           method='nearest')].mlat)
#         that = self.df['mlat']
#         return "{:.1f}".format(that.iloc[
#             that.index.get_loc(dsdf.tz_localize(None),
#                                method='nearest')])

#     def tickFormat_mlt(self, x, pos=None):
#         dsdf = pd.Timestamp(mdates.num2date(x))
#         return "{:.1f}".format(self.df.iloc[
#             self.df.index.get_loc(dsdf.tz_localize(None),
#                                   method='nearest')].mlt)
