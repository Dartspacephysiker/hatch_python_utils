import numpy as np
from scipy.stats import mode, gaussian_kde
from scipy.optimize import minimize, shgo
import warnings

def fastrobustmode(data):
    """
    From "On a fast, robust estimator of the mode: Comparisons to other robust estimators with applications"
    2006, Bickel and Frühwirth
    https://www.sciencedirect.com/science/article/pii/S0167947305001581
    """
    # data = np.random.normal(size=4)
    data = np.sort(data) 
    Nsamp = data.size

    # If there is only one number in the sample(n=1), then that number is the
    # estimate ofthe mode. Returnx1and stop.
    if data.size == 1:
        return data[0]

    # 2.  If there are two numbers in the sample(n=2), then the estimate of the mode
    # is themean of those numbers. Return(x1+x2)/2 and stop.
    elif data.size == 2:
        return np.mean(data)

    # 3.  If there are three numbers in the sample(n=3), then the estimate of the
    # mode is the meanof the two numbers that are closest together: ifx2−x1<x3−x2,
    # then return(x1+x2)/2and stop, but ifx2−x1<x3−x2, then return(x2+x3)/2 and
    # stop. Otherwise, ifn=3andx2−x1=x3−x2, then the mode is estimated to bex2;
    # returnx2and stop.
    elif data.size == 3:
        diffs = np.diff(data)
        if np.isclose(*diffs):
            return data[1]
        elif diffs[0] < diffs[1]:
            return np.mean(data[:2])
        elif diffs[0] > diffs[1]:
            return np.mean(data[1:])

    # 4.  If there are four or more numbers in the sample(n4), then the algorithm
    # will berecursively applied to about half of the sample. Set the minimum interval
    # width to therange of the sample:wmin=xn−x1. This will be replaced by the width
    # of the smallestinterval that contains approximately half of the values of the
    # sample. SetNequal to thesmallest integer greater than or equal ton/2. (This
    # algorithm could be generalized bysettingNequal to the smallest integer greater
    # than or equal to n, where 0< <1, butestimators have not been studied with
    # =1/2.) Fori=1,2,3,...,n−N+1, do thefollowing two steps:(a) Set the
    # widthw=xi+N−1−xi.(b) Ifw<wmin, then setwmin=wand setj=i.
    wmin = np.max(data)-np.min(data)
    N = np.int64(np.ceil(Nsamp/2))
    for i in range(N):
        w = data[N+i]-data[i]
        if w < wmin:
            wmin = w
            j = i
    return fastrobustmode(data[j:N+j])


def kde(array,
        cut_down=True,
        bw_method='scott'):
    if cut_down:
        bins, counts = np.unique(array, return_counts=True)
        if np.max(counts) == 1:
            warnings.warn("KDE: Cannot use 'cut_down' when all array values "
                          "appear only once (i.e., are unique)!")
        else:
            f_mean = counts.mean()
            f_above_mean = bins[counts > f_mean]
            bounds = [f_above_mean.min(), f_above_mean.max()]
            array = array[np.bitwise_and(bounds[0] < array, array < bounds[1])]

    return gaussian_kde(array, bw_method=bw_method)


# From "Modal Series #1 — How To Estimate The Global Mode From Data Sample"
# https://towardsdatascience.com/modal-series-1-how-to-estimate-the-global-mode-from-data-sample-e3e41380bfdb
def mode_estimation(array, cut_down=True, bw_method='scott'):
    kernel = kde(array, cut_down=cut_down, bw_method=bw_method)
    bounds = np.array([[array.min(), array.max()]])
    if len(array) <= 1000:
        nPoints = 100*len(array)
    elif len(array) <= 10000:
        nPoints = 10*len(array)
    else:
        nPoints = len(array)

    results = shgo(lambda x: -kernel(x)[0], bounds=bounds, n=nPoints)
    return results.x[0]


def mode_explicit(array, cut_down=True, bw_method='scott'):
    kernel = kde(array, cut_down=cut_down, bw_method=bw_method)
    height = kernel.pdf(array)
    return array[np.argmax(height)]


def refined_mode_estimation(array, cut_down=True, bw_method='scott'):
    kernel = kde(array, cut_down=cut_down, bw_method=bw_method)
    height = kernel.pdf(array)
    x0 = array[np.argmax(height)]
    span = array.max() - array.min()
    dx = span / 4
    bounds = np.array([[x0 - dx, x0 + dx]])
    linear_constraint = [{'type': 'ineq', 'fun': lambda x:  x - 0.5}]
    results = minimize(lambda x: -kernel(x)[0], x0=x0, bounds=bounds,
                       constraints=linear_constraint)
    return results.x[0]


def mode_estimation_suite(array,
                          do_rounding_trick=True,
                          kernelresampvals=None,
                          roundval=None,
                          **kws):
    """
    kernelresampvals: Where to resample the KDE estimate for obtaining the "mode_explicit" estimate
                      Can be one of [None,Integer,Array of quantiles between 0 and 1 to calculate]
    """

    verbose = kws['verbose'] if 'verbose' in kws else True

    if kernelresampvals is None:
        # print("No kernel resample values provided! Just doing a linspace for "
        #       "you")
        # kernelresampvals = np.linspace(array.min(), array.max(),
        #                                np.clip(10000, 0, array.size))
        print("No kernel resample values provided! Using quantiles")

        kernelresampvals = np.quantile(array,np.arange(0, 1.001, 0.001))

    elif np.issubdtype(type(kernelresampvals), np.integer):
        kernelresampvals = np.linspace(array.min(), array.max(),
                                       np.clip(kernelresampvals, 0, array.size))
    elif isinstance(kernelresampvals, np.ndarray):
        if verbose:
            print("Calculating kernelresampvals from provided quantile array ...")
        if (np.min(kernelresampvals) >= 0) & (np.min(kernelresampvals) <= 1):
            kernelresampvals = np.quantile(array,kernelresampvals)

    if verbose:
        import time
        fmtstr = "{:40s}: {:.5f}"
        t0 = time.time()

    if do_rounding_trick:
        log10minmax = np.log10([np.min(np.abs(array)),
                                np.max(np.abs(array))])
        if np.diff(log10minmax) >= 4:
            warnings.warn(">= 4 orders of magnitude difference between"
                          "min/max of provided array! Rounding trick for"
                          "'cut_down' may not work!")
        if roundval is None:
            roundval = -1*np.int(np.floor(np.log10(np.median(array))))+2
            print(f"Using roundval == {roundval}")
        narbef = array.size
        array = np.around(array, roundval)
        naraft = array.size
        if verbose:
            tround = time.time()
            print(fmtstr.format("Roundtrick execution time", tround-t0))
            tround = t0
            print("Shrunk array {:.3f}% ({:d} -> {:d})".format(
                (narbef-naraft)/narbef*100,narbef,naraft))

    # Scipy
    scipy = mode(array)[0][0]

    # Get KDE
    kernel = kde(array, **kws)
    if verbose:
        tkde = time.time()
        print(fmtstr.format("KDE execution time", tkde-t0))

    # From mode_estimation
    bounds = np.array([[array.min(), array.max()]])
    if len(array) <= 1000:
        nPoints = 100*len(array)
    elif len(array) <= 10000:
        nPoints = 10*len(array)
    else:
        nPoints = np.cut(10000, 0, len(array))

    mode_est = shgo(lambda x: -kernel(x)[0], bounds=bounds,
                    n=nPoints).x[0]
    # mode_est = 0
    if verbose:
        tmode = time.time()
        print(fmtstr.format("Mode (vanilla) execution time", tmode-tkde))

    # From mode_explicit
    height = kernel.pdf(kernelresampvals)
    x0 = array[np.argmax(height)]
    if verbose:
        texp = time.time()
        print(fmtstr.format("Mode (explicit) execution time", texp-tmode))

    # From refined_mode_estimation
    span = array.max() - array.min()
    dx = span / 4
    bounds = np.array([[x0 - dx, x0 + dx]])
    linear_constraint = [{'type': 'ineq', 'fun': lambda x:  x - 0.5}]
    refined_mode_est = minimize(lambda x: -kernel(x)[0], x0=x0, bounds=bounds,
                                constraints=linear_constraint).x[0]
    if verbose:
        tref = time.time()
        print(fmtstr.format("Mode (refined) execution time", tref-texp))

    suite = dict(scipy=scipy,
                 # estimate=mode_estimation(array, **kws),
                 # explicit=x0licit(array, **kws),
                 # refined=refined_mode_estimation(array, **kws))
                 estimate=mode_est,
                 explicit=x0,
                 refined=refined_mode_est)

    return suite


def mode_estimation_scott_vs_silverman(array,
                                       do_rounding_trick=True,
                                       kernelresampvals=None,
                                       roundval=None,
                                       cut_down=True,
                                       **kws):

    return dict(
        scott=mode_estimation_suite(array,
                                    cut_down=cut_down,
                                    do_rounding_trick=do_rounding_trick,
                                    kernelresampvals=kernelresampvals,
                                    roundval=roundval,
                                    bw_method='scott',
                                    **kws),
        silverman=mode_estimation_suite(array,
                                        cut_down=cut_down,
                                        do_rounding_trick=do_rounding_trick,
                                        kernelresampvals=kernelresampvals,
                                        roundval=roundval,
                                        bw_method='silverman',
                                        **kws))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    def generate_1D_data(size, mu=0.5, sigma=1, round=3):
        data = np.random.normal(mu, sigma, size)
        return np.around(data, round)

    mu = 0.5
    sigma = 0.2
    data = generate_1D_data(int(1e3), mu=mu, sigma=sigma)

    mode_ = mode(data)[0][0]
    _mode = mode_explicit(data)
    _mode_ = mode_estimation(data)
    _imode_ = refined_mode_estimation(data)
    
    print(f'actual mode is at {mu}')
    print(f'scipy mode: {mode_}')
    print(f'mode_explicit mode: {_mode}')
    print(f'mode_estimation mode: {_mode_}')
    print(f'refined_mode_estimation mode: {_imode_}')
    
    # Fig 1
    fig1,ax1 = plt.subplots(1,1)
    count, bins, ignored = plt.hist(data, 50, density=True, color='b', label='data')
    f_mean = count.mean()
    f_above_mean = bins[:-1][count > f_mean]
    bounds = [f_above_mean.min(), f_above_mean.max()]

    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), color='k', label='distribution')
    plt.axvline(bounds[0], label='lower bound', color='c')
    plt.axvline(bounds[1], label='upper bound', color='m')
    plt.axhline(f_mean , label='mean frequency', color='r')
    plt.title('normal distribution with range boundes defined by mean of frequencies')
    plt.xlabel('data point values')
    plt.ylabel('data point frequencies')
    plt.legend()

    # Fig 2
    fig2,ax2 = plt.subplots(1,1,figsize=(20, 10))

    bins = np.linspace(data.min(), data.max(), 100)
    plt.hist(data, bins=100, density=True, label='data', color='b')
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) ), linewidth=4, color='r', label='distribution')
    plt.axvline(mode_, linewidth=2, label='scipy mode', color='m')
    plt.axvline(_mode, linewidth=2, label='explicit mode', color='c')
    plt.axvline(_mode_, linewidth=2, label='KDE mode', color='g')
    plt.axvline(_imode_, linewidth=2, label='refined mode', color='y')
    plt.title('normal distribution with the estimated mode methods')
    plt.xlabel('data point values')
    plt.ylabel('data point frequencies')
    plt.grid()
    plt.legend()
    plt.show()
