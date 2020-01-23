import numpy as np

def Mu(x):
    """
    Maximum-likelihood estimator of mu parameter of log-normal distribution
    """
    return np.mean(np.log(x))

def Sigma(x):
    """
    Maximum-likelihood estimator of sigma parameter of log-normal distribution
    """
    return np.std(np.log(x)))

def Variance(x):
    """
    Maximum-likelihood estimator of variance of log-normal distribution
    """
    sigmaSq = sigma(x)**2
    return (np.exp(sigmaSq)-1)*np.exp(2*mu(x)+sigmaSq)

def Std(x):
    """
    Maximum-likelihood estimator of std deviation of log-normal distribution
    """
    return np.sqrt(variance(x))

def Median(x):
    """
    Maximum-likelihood estimator of median of log-normal distribution
    """
    return np.exp(mu(x))

def Mean(x):
    return np.exp(mu(x)+sigma(x)**2./2)
