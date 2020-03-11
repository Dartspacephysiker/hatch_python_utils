# 2019/10/28
import numpy as np
from scipy.stats import binom as binomDist
from scipy.stats import median_absolute_deviation as MAD

def calc_binned_stats_dict(xvals,yvals,sjekkBins,
                           xlabel='',
                           doWeightedSum=False,
                           weights=None,
                           doMedian=True,
                           doMedian_CI95=True,
                           doMean=False,
                           dolog10=False,
                           doVariance=False,
                           doVariance_CI95=False,
                           doVariance__showStdDev=False,
                           treat_as_periodic=False,
                           periodic_Xbounds=None):
    """
    For use with a particular set of indices, xlabel, and ylabel, I have done sånt:
    xvals = dfUse[baseinds][xlabel]
    yvals = dfUse[baseinds][ylabel]
    """

    assert not doWeightedSum,"doWeightedSum kw not implemented!"

    if xlabel == 'scaledtidoffset':
        treat_as_periodic = True
        periodic_Xbounds = np.array([0,4])
    else:
        treat_as_periodic = False
        periodic_Xbounds = None
    
    if doWeightedSum:

        midAll, medAll, Q1All, Q3All, bincounts = bin_weighted_sum_getter(xvals,yvals,weights,
                                                                          binlines=sjekkBins,
                                                                          include_bincounts=True,
                                                                          treat_as_periodic=treat_as_periodic,
                                                                          periodic_Xbounds=periodic_Xbounds)

    elif doMedian:
        if doMedian_CI95:
            print("CI95")
            # midAll, medAll, CI95_LOW, CI95_HIGH, bincounts = bin_median_CI95_getter(xvals,
            midAll, medAll, Q1All, Q3All, bincounts = bin_median_CI95_getter(xvals,
                                                                             yvals,
                                                                             binlines=sjekkBins,
                                                                             include_bincounts=True,
                                                                             treat_as_periodic=treat_as_periodic,
                                                                             periodic_Xbounds=periodic_Xbounds)
        else:
            midAll, medAll, Q1All, Q3All, bincounts = bin_median_Q1_Q3_getter(xvals,
                                                                              yvals,
                                                                              binlines=sjekkBins,
                                                                              include_bincounts=True,
                                                                              treat_as_periodic=treat_as_periodic,
                                                                              periodic_Xbounds=periodic_Xbounds)
    elif doVariance:
        print("VARIANCE (maybelog10)")
        if doVariance_CI95:
            print("CI95")
            
            midAll, medAll, Q1All, Q3All, bincounts = bin_variance_CI95_getter(xvals,
                                                                               yvals,
                                                                               binlines=sjekkBins,
                                                                               include_bincounts=True,
                                                                               treat_as_periodic=treat_as_periodic,
                                                                               periodic_Xbounds=periodic_Xbounds)
        if doVariance__showStdDev:
            print("Convert to stdDev")
            statname = 'stddev'
    
            medAll = np.sqrt(medAll)
    
            Q1All = np.sqrt(Q1All)
            Q3All = np.sqrt(Q3All)
    
    elif doMean:
        print("Real Lindis")
        if dolog10:
            midAll, medAll, Q1All, Q3All,bincounts = bin_mean_pmstddev_getter(xvals,
                                                                              np.log10(yvals),
                                                                              binlines=sjekkBins,
                                                                              include_bincounts=True,
                                                                              treat_as_periodic=treat_as_periodic,
                                                                              periodic_Xbounds=periodic_Xbounds)
    
        else:
            midAll, medAll, Q1All, Q3All,bincounts = bin_mean_pmstddev_getter(xvals,
                                                                              yvals,
                                                                              binlines=sjekkBins,
                                                                              include_bincounts=True,
                                                                              treat_as_periodic=treat_as_periodic,
                                                                              periodic_Xbounds=periodic_Xbounds)
    elif doMAD:
        print("MAD")

        midAll, medAll, Q1All, Q3All, bincounts = bin_MAD_getter(xvals,
                                                                 yvals,
                                                                 binlines=sjekkBins,
                                                                 include_bincounts=True,
                                                                 treat_as_periodic=treat_as_periodic,
                                                                 periodic_Xbounds=periodic_Xbounds)
    curmean = np.mean(medAll)
    curstd = np.std(medAll)
    
    if doMean:
        errbars = np.vstack([(medAll-Q1All)/np.sqrt(bincounts),
                             (Q3All-medAll)/np.sqrt(bincounts)])
    else:
        errbars = np.vstack([(medAll-Q1All),
                             (Q3All-medAll)])
    
    return {'binCtr'                  :midAll,
            xlabel+'means'            :medAll,
            xlabel+'errbars'          :errbars,
            'meanOf'+xlabel+'means'   :curmean,
            'stddevOf'+xlabel+'means' :curstd,
            'count'                   :bincounts}


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


def CI_alpha_var(a, alpha,ddof=1):
    """
    # Get confidence intervals for population variance
    Based on https://www.thoughtco.com/interval-for-a-population-variance-3126221
    """
    from scipy.stats import chi2 as chi2Dist
    assert len(a.shape) == 1,"Support for multi-dim arrays not implemented!"

    n = a.size
    # df = n-1                    # degrees of freedom to use for chiSq dist
    sSq = np.var(a,ddof=ddof)
    
    A = chi2Dist.ppf(alpha/2,n-1)
    B = chi2Dist.ppf(1-alpha/2,n-1)

    return np.array([(n-1)*sSq/B,(n-1)*sSq/A])
    

# def CI_95_weightedsum(a,weights=None):

# SEE https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
# Apparently no standard definition of standard error of weighted mean
#     wgtsum = np.sum(weights)
#     wgtmean = np.sum(a*weights)/wgtsum

#     # Get variance assuming "frequency weights" and not "reliability weights" (see https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance)
#     stddev = np.sqrt(np.sum(weights*(a-wgtmean)**2)/(wgtsum-1))
    
#     return wgtmean,wgtmean


def CI_95_var(a,ddof=1):
    """
    # Get confidence intervals for population variance
    Based on https://www.thoughtco.com/interval-for-a-population-variance-3126221
    """

    return CI_alpha_var(a, 0.05,ddof=1)

    # from scipy.stats import chi2 as chi2Dist
    # assert len(a.shape) == 1,"Support for multi-dim arrays not implemented!"

    # n = a.size
    # # df = n-1                    # degrees of freedom to use for chiSq dist
    # sSq = np.var(a)
    
    # A = chi2Dist.ppf(alpha/2,n-1)
    # B = chi2Dist.ppf(1-alpha/2,n-1)

    # return np.array([(n-1)*sSq/B,(n-1)*sSq/A])


def CI_95_quantile(a,q):
    """
    Refs
    ===
    https://www-users.york.ac.uk/~mb55/intro/cicent.htm (who references Conover, W.J. (1980) Practical Nonparametric Statistics John Wiley and Sons, New York.)
    """
    return np.sort(a)[CI_95_quantile_rankvalue(a.size,q)]


def bin_ind_getter(X,
                   binlines=None,
                   reference_inds=None,
                   treat_as_periodic=False,
                   periodic_Xbounds=None):
    
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

    binlines_are_Nx2 = len(binlines.shape) == 2
    if binlines_are_Nx2:
        nBins = binlines.shape[0]
    else:
        nBins = len(binlines)-1

    for i in range(nBins):
        if binlines_are_Nx2:
            tmpbinedges = binlines[i,:]
        else:
            tmpbinedges = binlines[i], binlines[i+1]

        # PRE-20200217 OLD
        # this = np.where((X >= tmpbinedges[0]) & (
        #     X <= tmpbinedges[1]))[0]
        # 20200217 NYE
        this = get_periodic_bin_edge_indices(X,tmpbinedges,
                                             treat_as_periodic=treat_as_periodic,
                                             periodic_Xbounds=periodic_Xbounds)

        inds.append(reference_inds[this])

    return inds

def get_periodic_bin_edge_indices(X,binedges,treat_as_periodic=False,
                                  periodic_Xbounds=None):
    """
    indices = get_periodic_bin_edge_indices(X,binedges,Xbounds)
    Given two binedges and an array X, get the indices of X that fall within the two binedges
    """
    DEBUG = False

    if treat_as_periodic:
        assert binedges[0] < binedges[1]
        assert periodic_Xbounds[0] < periodic_Xbounds[1]

        extends_below = binedges[0] < periodic_Xbounds[0]
        extends_above = binedges[1] > periodic_Xbounds[1]
        extends_both_ways = extends_below and extends_above
        if extends_both_ways:
            print("PERIODIC BOTH SIDES")
            assert ~extends_both_ways,"NOT IMPLEMENTED"
        elif extends_below:

            extension_below = periodic_Xbounds[0]-binedges[0]  # Guaranteed a positive quantity
            assert extension_below >= 0,"This quantity skal være positiv!"

            both_extend_below = binedges[1] < periodic_Xbounds[0]

            if both_extend_below:
                binedges = binedges + (periodic_Xbounds[1]-periodic_Xbounds[0])

                indlets = np.where((X >= binedges[0]) & (X < binedges[1]))[0]

            else:
                indlets = X < binedges[1]
                indlets = np.where(indlets | (X >= (periodic_Xbounds[1] - extension_below)))[0]

                if DEBUG:
                    print("EXTENDS BELOW: {:.2f} - {:.2f}".format(periodic_Xbounds[1] - extension_below,binedges[1]))

        elif extends_above:

            extension_above = binedges[1]-periodic_Xbounds[1]  # Guaranteed a positive quantity
            assert extension_above >= 0,"This quantity skal være positiv!"

            both_extend_above = binedges[0] > periodic_Xbounds[1]

            if both_extend_above:
                binedges = binedges - (periodic_Xbounds[1]-periodic_Xbounds[0])
                indlets = np.where((X >= binedges[0]) & (X < binedges[1]))[0]
            else:
                indlets = X >= binedges[0]
                indlets = np.where(indlets | (X < (periodic_Xbounds[0] + extension_above)))[0]

                if DEBUG:
                    print("EXTENDS ABOVE: {:.2f} - {:.2f}".format(binedges[0],periodic_Xbounds[0] + extension_above))

        else:
            indlets = np.where((X >= binedges[0]) & (X < binedges[1]))[0]

    else:
        indlets = np.where((X >= binedges[0]) & (X < binedges[1]))[0]

    return indlets

def bin_median_getter(X, Y, binlines=None,
                      statfunc=np.median,
                      statfunc_kws=dict(),
                      include_CI95=False,
                      CI95func_kws=dict(),
                      variance__estimate_population_variance=False,
                      treat_as_periodic=False,
                      periodic_Xbounds=None):
    """
    binmid, vals[, CI95_LOW, CI95_HIGH] = bin_median_getter(X, Y, binlines=None, statfunc=np.median[,include_CI95=True])
    """

    isWeightedMean = False

    if include_CI95:

        if statfunc == np.median:
            print("95% CI for median") 
            CI_95_func = CI_95_median

        elif statfunc == np.var:
            print("95% CI for variance") 
            CI_95_func = CI_95_var

        elif statfunc == MAD:
            print("Can't do CI for MAD!")
            CI_95_func = CI_95_var

        elif statfunc == np.sum:

            assert 'weights' in CI95func_kws.keys(),"must provide weights!"

            isWeightedMean = True
            weights = CI95func_kws['weights'].copy()

            def weightedmean(Y,weights=None):
                return np.sum(Y*weights)/np.sum(weights)
            statfunc = weightedmean

            CI_95_func = CI_95_weightedsum

    # Are we doing sample variance or population variance?
    if statfunc == np.var:
        if variance__estimate_population_variance:
            print("Assuming you want an estimate of the population variance, not sample variance")
            statfunc = lambda x: np.var(x,ddof=1)  # ddof causes us to divide by N-1 instead of N


    if binlines is None:
        binlines = np.arange(0, 70.1, 2.5)
    elif isinstance(binlines,int):
        binlines = np.linspace(np.min(X),np.max(X),binlines)

    binlines_are_Nx2 = len(binlines.shape) == 2
    if binlines_are_Nx2:
        nBins = binlines.shape[0]
        binmid = np.mean(binlines,axis=1)
    else:
        nBins = len(binlines)-1
        binmid = np.diff(binlines)/2.+binlines[0:-1]

    vals = np.zeros(binmid.shape)*np.nan

    if include_CI95:
        CI95_LOW = np.zeros(binmid.shape)*np.nan
        CI95_HIGH = np.zeros(binmid.shape)*np.nan

    # for i in range(len(binlines)-1):
    for i in range(nBins):
        if binlines_are_Nx2:
            tmpbinedges = binlines[i,:]
        else:
            tmpbinedges = binlines[i], binlines[i+1]

        # this = np.where((X >= tmpbinedges[0]) & (
        #     X <= tmpbinedges[1]))[0]
        this = get_periodic_bin_edge_indices(X,tmpbinedges,
                                             treat_as_periodic=treat_as_periodic,
                                             periodic_Xbounds=periodic_Xbounds)

        if len(this) > 0:

            if isWeightedMean:
                statfunc_kws['weights'] = weights[this]
                CI95func_kws['weights'] = weights[this]

            vals[i] = statfunc(Y[this],**statfunc_kws)

            if include_CI95:
                CI95 = CI_95_func(Y[this],**CI95func_kws)
                CI95_LOW[i] = CI95[0]
                CI95_HIGH[i] = CI95[1]

        # print(*tmpbinedges,vals[i])
    if include_CI95:
        return binmid, vals, CI95_LOW, CI95_HIGH
    else:
        return binmid, vals

def bin_median_CI95_getter(X, Y, binlines=None,include_bincounts=True,
                           treat_as_periodic=False,
                           periodic_Xbounds=None):
    """
    binmidpts, median, CI95_LOW, CI95_HIGH[, bincounts] = bin_median_CI95_getter(X, Y, binlines=None[, include_bincounts=True])
    """

    midtpunkter, binVerdi, CI95_LOW, CI95_HIGH = bin_median_getter(X, Y, binlines=binlines,
                                                                   statfunc=np.median,
                                                                   include_CI95=True,
                                                                   treat_as_periodic=treat_as_periodic,
                                                                   periodic_Xbounds=periodic_Xbounds)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH, bincounts

    else:
        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH


def bin_weighted_sum_getter(X, Y, weights,
                            binlines=None,
                            include_bincounts=True,
                            treat_as_periodic=False,
                            periodic_Xbounds=None):
    """
    binmidpts, wgtedsum, NOTKNOW_LOW, NOTKNOW_HIGH[, bincounts] = bin_weighted_sum_getter(X, Y, weights, binlines=None[, include_bincounts=True])
    """

    midtpunkter, binVerdi, = bin_median_getter(X, Y,
                                               binlines=binlines,
                                               statfunc=np.sum,
                                               CI95func_kws=dict(weights=weights),
                                               include_CI95=False,
                                               treat_as_periodic=treat_as_periodic,
                                               periodic_Xbounds=periodic_Xbounds)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y,
                                                   binlines=binlines,
                                                   statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH, bincounts

    else:
        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH


def bin_MAD_getter(X, Y, binlines=None,include_bincounts=True,
                   treat_as_periodic=False,
                   periodic_Xbounds=None):
    """
    binmidpts, MAD, CI95_LOWBOGUS, CI95_HIGHBOGUS[, bincounts] = bin_MAD_getter(X, Y, binlines=None[, include_bincounts=True])
    """

    midtpunkter, binVerdi, CI95_LOW, CI95_HIGH = bin_median_getter(X, Y, binlines=binlines,
                                                                   statfunc=MAD,
                                                                   include_CI95=True,
                                                                   treat_as_periodic=treat_as_periodic,
                                                                   periodic_Xbounds=periodic_Xbounds)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH, bincounts

    else:
        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH


def bin_median_Q1_Q3_getter(X, Y, binlines=None,include_bincounts=True,
                            treat_as_periodic=False,
                            periodic_Xbounds=None):
    """
    binmidpts, median, Q1, Q3 = bin_median_Q1_Q3_getter(X, Y, binlines=None)
    """

    Q1func = lambda x :np.quantile(x,0.25)
    Q3func = lambda x :np.quantile(x,0.75)

    midtpunkter, binVerdi = bin_median_getter(X, Y, binlines=binlines,statfunc=np.median,
                                              treat_as_periodic=treat_as_periodic,
                                              periodic_Xbounds=periodic_Xbounds)
    midtpunkter, Q1 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q1func,
                                        treat_as_periodic=treat_as_periodic,
                                        periodic_Xbounds=periodic_Xbounds)
    midtpunkter, Q3 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q3func,
                                        treat_as_periodic=treat_as_periodic,
                                        periodic_Xbounds=periodic_Xbounds)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, Q1, Q3, bincounts

    else:
        return midtpunkter, binVerdi, Q1, Q3


def bin_mean_pmstddev_getter(X, Y, binlines=None,include_bincounts=True,
                             treat_as_periodic=False,
                             periodic_Xbounds=None):
    """
    binmidpts, mean, mean-1stddev, mean+1stddev = bin_mean_pmstddev_getter(X, Y, binlines=None, include_bincounts=True)
    """

    midtpunkter, binVerdi = bin_median_getter(X, Y, binlines=binlines,statfunc=np.mean,
                                              treat_as_periodic=treat_as_periodic,
                                              periodic_Xbounds=periodic_Xbounds)
    midtpunkter, stdDev = bin_median_getter(X, Y, binlines=binlines,statfunc=np.std,
                                            treat_as_periodic=treat_as_periodic,
                                            periodic_Xbounds=periodic_Xbounds)

    Q1 = binVerdi - stdDev
    Q3 = binVerdi + stdDev
    # midtpunkter, Q1 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q1func)
    # midtpunkter, Q3 = bin_median_getter(X, Y, binlines=binlines,statfunc=Q3func)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y, binlines=binlines,statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, Q1, Q3, bincounts

    else:
        
        return midtpunkter, binVerdi, Q1, Q3

def bin_variance_CI95_getter(X, Y,
                             binlines=None,
                             include_bincounts=True,
                             treat_as_periodic=False,
                             periodic_Xbounds=None):
    """
    binmidpts, variance, CI95_LOW, CI95_HIGH[, bincounts] = bin_variance_CI95_getter(X, Y, binlines=None[, include_bincounts=True])
    """

    midtpunkter, binVerdi, CI95_LOW, CI95_HIGH = bin_median_getter(X, Y,
                                                                   binlines=binlines,
                                                                   statfunc=np.var,
                                                                   include_CI95=True,
                                                                   variance__estimate_population_variance=True,
                                                                   treat_as_periodic=treat_as_periodic,
                                                                   periodic_Xbounds=periodic_Xbounds)

    if include_bincounts:
        midtpunkter, bincounts = bin_median_getter(X, Y,
                                                   binlines=binlines,
                                                   statfunc=np.size,
                                                   treat_as_periodic=treat_as_periodic,
                                                   periodic_Xbounds=periodic_Xbounds)

        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH, bincounts

    else:
        return midtpunkter, binVerdi, CI95_LOW, CI95_HIGH


