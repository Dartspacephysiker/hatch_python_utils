# 2019/10/28
import numpy as np


def bin_median_getter(X, Y, binlines=None, statfunc=np.median):
    """
    bin_median_getter(X, Y, binlines=None, statfunc=np.median)
    """
    if binlines is None:
        binlines = np.arange(0, 70.1, 2.5)
    elif isinstance(binlines,int):
        binlines = np.linspace(np.min(X),np.max(X),binlines)

    binmid = np.diff(binlines)/2.+binlines[0:-1]
    vals = np.zeros(binmid.shape)*np.nan

    for i in range(len(binlines)-1):
        tmpbinedges = binlines[i], binlines[i+1]
        this = np.where((X >= tmpbinedges[0]) & (
            X <= tmpbinedges[1]) & np.isfinite(Y))[0]
        if len(this) > 0:
            vals[i] = statfunc(Y[this])
        # print(*tmpbinedges,vals[i])
    return binmid, vals

def bin_median_Q1_Q3_getter(X, Y, binlines=None):
    """
    binmidpts, median, Q1, Q3 = bin_median_Q1_Q3_getter(X, Y, binlines=None, statfunc=np.median)
    """

    Q1func = lambda x :np.quantile(x,0.25)
    Q3func = lambda x :np.quantile(x,0.75)

    midtpunkter, binVerdi = bin_median_getter(X, Y, binlines=binlines,statfunc=np.median)
    midtpunkter, Q1 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q1func)
    midtpunkter, Q3 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q3func)

    return midtpunkter, binVerdi, Q1, Q3

    # if binlines is None:
    #     binlines = np.arange(0, 70.1, 2.5)
    # elif isinstance(binlines,int):
    #     binlines = np.linspace(np.min(X),np.max(X),binlines)

    # binmid = np.diff(binlines)/2.+binlines[0:-1]
    # vals = np.zeros(binmid.shape)*np.nan

    # for i in range(len(binlines)-1):
    #     tmpbinedges = binlines[i], binlines[i+1]
    #     this = np.where((X >= tmpbinedges[0]) & (
    #         X <= tmpbinedges[1]) & np.isfinite(Y))[0]
    #     if len(this) > 0:
    #         vals[i] = statfunc(Y[this])
    #     # print(*tmpbinedges,vals[i])
    # return binmid, vals

def bin_mean_pmstddev_getter(X, Y, binlines=None,include_bincounts=True):
    """
    binmidpts, mean, mean-1stddev, mean+1stddev = bin_mean_pmstddev_getter(X, Y, binlines=None, include_bincounts=True)
    """

    midtpunkter, binVerdi = bin_median_getter(X, Y, binlines=binlines,statfunc=np.mean)
    midtpunkter, stdDev = bin_median_getter(X, Y, binlines=binlines,statfunc=np.std)

    Q1 = binVerdi - stdDev
    Q3 = binVerdi + stdDev
    # midtpunkter, Q1 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q1func)
    # midtpunkter, Q3 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q3func)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size)

        return midtpunkter, binVerdi, Q1, Q3, bincounts

    else:
        
        return midtpunkter, binVerdi, Q1, Q3

    # if binlines is None:
    #     binlines = np.arange(0, 70.1, 2.5)
    # elif isinstance(binlines,int):
    #     binlines = np.linspace(np.min(X),np.max(X),binlines)

    # binmid = np.diff(binlines)/2.+binlines[0:-1]
    # vals = np.zeros(binmid.shape)*np.nan

    # for i in range(len(binlines)-1):
    #     tmpbinedges = binlines[i], binlines[i+1]
    #     this = np.where((X >= tmpbinedges[0]) & (
    #         X <= tmpbinedges[1]) & np.isfinite(Y))[0]
    #     if len(this) > 0:
    #         vals[i] = statfunc(Y[this])
    #     # print(*tmpbinedges,vals[i])
    # return binmid, vals
