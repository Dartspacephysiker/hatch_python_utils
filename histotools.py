# 2019/10/28
import numpy as np
from scipy.stats import binom as binomDist

# THIS ARTICLE PROB DOES IT ALL
#Confidence interval for quantiles and percentiles (doi: 10.11613/BM.2019.010101) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6294150/

# Good discussion of median CI: https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median
def CI_95_median_rankvalue(N: int):
    """
    L, U = CI_95_median_rankvalue(N: int)
    Based on https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8
    """
    return np.array([np.int64(np.round(N/2-1.96*np.sqrt(N)/2)),
                     np.int64(np.round(1+N/2+1.96*np.sqrt(N)/2))])


def CI_95_quantile_rankvalue(N: int,q: float):
    """
    L, U = CI_95_quantile_rankvalue(N,q)

    Refs
    ===
    https://www-users.york.ac.uk/~mb55/intro/cicent.htm (who references Conover, W.J. (1980) Practical Nonparametric Statistics John Wiley and Sons, New York.)

    """
    return np.array([np.int64(np.round(N*q-1.96*np.sqrt(N*q*(1-q)))),
                     np.int64(np.round(N*q+1.96*np.sqrt(N*q*(1-q))))])


def CI_alpha_median_rankvalue(N: int,alpha: float):
    """
    L, U = CI_alpha_median_rankvalue(N,alpha)

    Refs
    ===
    https://stats.stackexchange.com/questions/122001/confidence-intervals-for-median

    """
    from scipy.stats import binom as binomDist

    p = 0.5                     # median
    L = binomDist.ppf(alpha/2,N,p)
    U = binomDist.ppf(1-alpha/2,N,p)

    return L,U

def CI_95_median(a):
    """
    Based on https://www.ucl.ac.uk/child-health/short-courses-events/about-statistical-courses/research-methods-and-statistics/chapter-8-content-8
    """
    return np.sort(a)[CI_95_median_rankvalue(a.size)]


def CI_95_quantile(a,q):
    """
    Refs
    ===
    https://www-users.york.ac.uk/~mb55/intro/cicent.htm (who references Conover, W.J. (1980) Practical Nonparametric Statistics John Wiley and Sons, New York.)
    """
    return np.sort(a)[CI_95_quantile_rankvalue(a.size,q)]


def bin_ind_getter(X, binlines=None,
                   reference_inds=None):
    
    """
    binned_indices = bin_ind_getter(X, binlines)

    Parameters
    ===
    X               : array-like, stuff to be binned

    Keywords
    ========
    binlines       (default None): boundary value of each bin (or integer)
    reference_inds (default None): array of indices (same size as X) into which you'd like binned_indices inserted
                                   This is useful if, for example, X is some array which you have already binned prior to entering this function
    """
    if binlines is None:
        binlines = np.arange(0, 70.1, 2.5)
    elif isinstance(binlines,int):
        binlines = np.linspace(np.min(X),np.max(X),binlines)

    if reference_inds is None:
        reference_inds = np.arange(X.size)
    else:
        assert reference_inds.size == X.size

    inds = []
    for i in range(len(binlines)-1):
        tmpbinedges = binlines[i], binlines[i+1]
        this = np.where((X >= tmpbinedges[0]) & (
            X <= tmpbinedges[1]))[0]
        inds.append(reference_inds[this])

    return inds

def bin_median_getter(X, Y, binlines=None, statfunc=np.median,include_CI95=False):
    """
    binmid, vals[, CI95_LOW, CI95_HIGH] = bin_median_getter(X, Y, binlines=None, statfunc=np.median[,include_CI95=True])
    """
    if binlines is None:
        binlines = np.arange(0, 70.1, 2.5)
    elif isinstance(binlines,int):
        binlines = np.linspace(np.min(X),np.max(X),binlines)

    binmid = np.diff(binlines)/2.+binlines[0:-1]
    vals = np.zeros(binmid.shape)*np.nan

    if include_CI95:
        CI95_LOW = np.zeros(binmid.shape)*np.nan
        CI95_HIGH = np.zeros(binmid.shape)*np.nan

    for i in range(len(binlines)-1):
        tmpbinedges = binlines[i], binlines[i+1]
        this = np.where((X >= tmpbinedges[0]) & (
            X <= tmpbinedges[1]))[0]
        if len(this) > 0:
            vals[i] = statfunc(Y[this])

            if include_CI95:
                CI95 = CI_95_median(Y[this])
                CI95_LOW[i] = CI95[0]
                CI95_HIGH[i] = CI95[1]

        # print(*tmpbinedges,vals[i])
    if include_CI95:
        return binmid, vals, CI95_LOW, CI95_HIGH
    else:
        return binmid, vals

def bin_median_CI95_getter(X, Y, binlines=None,include_bincounts=True):
    """
    binmidpts, median, CI95_LOW, CI95_HIGH[, bincounts] = bin_median_CI95_getter(X, Y, binlines=None, statfunc=np.median[,include_bincounts=True])
    """

    midtpunkter, binVerdi, CI95_LOW, CI95_HIGH = bin_median_getter(X, Y, binlines=binlines,statfunc=np.median,include_CI95=True)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size)

        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH, bincounts

    else:
        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH


def bin_median_Q1_Q3_getter(X, Y, binlines=None,include_bincounts=True):
    """
    binmidpts, median, Q1, Q3 = bin_median_Q1_Q3_getter(X, Y, binlines=None, statfunc=np.median)
    """

    Q1func = lambda x :np.quantile(x,0.25)
    Q3func = lambda x :np.quantile(x,0.75)

    midtpunkter, binVerdi = bin_median_getter(X, Y, binlines=binlines,statfunc=np.median)
    midtpunkter, Q1 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q1func)
    midtpunkter, Q3 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q3func)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size)

        return midtpunkter, binVerdi, Q1, Q3, bincounts

    else:
        return midtpunkter, binVerdi, Q1, Q3


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
