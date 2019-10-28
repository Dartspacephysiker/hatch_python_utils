# 2019/10/28
import numpy as np


def bin_median_getter(X, Y, binlines=None):

    if binlines is None:
        binlines = np.arange(0, 70.1, 2.5)

    binmid = np.diff(binlines)/2.+binlines[0:-1]
    vals = np.zeros(binmid.shape)*np.nan

    for i in range(len(binlines)-1):
        tmpbinedges = binlines[i], binlines[i+1]
        this = np.where((X >= tmpbinedges[0]) & (
            X <= tmpbinedges[1]) & np.isfinite(Y))[0]
        if len(this) > 0:
            vals[i] = np.median(Y[this])
        # print(*tmpbinedges,vals[i])
    return binmid, vals
