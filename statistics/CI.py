import numpy as np

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
    

def CI_95_weightedsum(a,weights=None):

    # SEE https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
    # Apparently no standard definition of standard error of weighted mean
    assert 2<0, "THIS DOESN'T WORK/ISN'T IMPLEMENTED"

    wgtsum = np.sum(weights)
    wgtmean = np.sum(a*weights)/wgtsum

    # Get variance assuming "frequency weights" and not "reliability weights" (see https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance)
    stddev = np.sqrt(np.sum(weights*(a-wgtmean)**2)/(wgtsum-1))
    
    return wgtmean,wgtmean


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


