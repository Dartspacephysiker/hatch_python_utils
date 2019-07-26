# 2019/07/26
# import colorednoise as cn
# Here, I completely rip off the colorednoise package to add an fmax argument
# As stated there, it's based on the algorithm in
# Timmer, J. and Koenig, M.: On generating power law noise. Astron. Astrophys. 300, 707-710 (1995)

from numpy import sqrt, newaxis
from numpy.fft import irfft, rfftfreq
from numpy.random import normal
from numpy import sum as npsum


def noisegen(color_or_beta, samples,
             sFreq=1,
             fmin=0,
             fmax=None,
             verbose=False):

    # haveBeta = False
    # haveColor = False

    try:
        _ = color_or_beta + 1
        beta = color_or_beta
        # haveBeta = True
    except:
        if not isinstance(color_or_beta, str):
            print("Must provide number or string to initialize noise!")
            return None

        color = color_or_beta
        # haveColor = True

        if color.lower() == 'pink':
            beta = 1
        elif color.lower() == 'brown':
            beta = 2
        elif color.lower() == 'white':
            beta = 0
        elif color.lower() == 'blue':
            beta = -1
        elif color.lower() == 'violet':
            beta = -2

        if verbose:
            print("Color: {:s} (beta = {:.2f})".format(color, beta))

    fmin2 = fmin/sFreq
    fmax2 = None
    if fmax is not None:
        fmax2 = fmax/sFreq
    else:
        fmax = 1
        fmax2 = 1

    if verbose:
        print("sFreq, fmin , fmax : {:.2f}, {:.2f}, {:.2f}".format(
            sFreq, fmin, fmax))
        print("sFreq, fmin2, fmax2: {:.2f}, {:.2f}, {:.2f}".format(
            sFreq, fmin2, fmax2))

    return powerlaw_psd_gaussian(beta,
                                 samples,
                                 fmin=fmin2,
                                 fmax=fmax2,
                                 verbose=verbose)


def powerlaw_psd_gaussian(exponent, size,
                          fmin=0,
                          fmax=None,
                          verbose=False):
    """Gaussian (1/f)**beta noise.

    Based on the algorithm in:
    Timmer, J. and Koenig, M.:
    On generating power law noise.
    Astron. Astrophys. 300, 707-710 (1995)

    Spencer's mod, which includes an fmax argument

    Normalised to unit variance

    Parameters:
    -----------

    exponent : float
        The power-spectrum of the generated noise is proportional to

        S(f) = (1 / f)**beta
        flicker / pink noise:   exponent beta = 1
        brown noise:            exponent beta = 2

        Furthermore, the autocorrelation decays proportional to lag**-gamma
        with gamma = 1 - beta for 0 < beta < 1.
        There may be finite-size issues for beta close to one.

    shape : int or iterable
        The output has the given shape, and the desired power spectrum in
        the last coordinate. That is, the last dimension is taken as time,
        and all other components are independent.

    fmin : float, optional
        Low-frequency cutoff.
        Default: 0 corresponds to original paper. It is not actually
        zero, but 1/samples.

    Returns
    -------
    out : array
        The samples.


    Examples:
    ---------

    # generate 1/f noise == pink noise == flicker noise
    >>> import colorednoise as cn
    >>> y = cn.powerlaw_psd_gaussian(1, 5)
    """

    # Make sure size is a list so we can iterate it and assign to it.
    try:
        size = list(size)
    except TypeError:
        size = [size]

    # The number of samples in each time series
    samples = size[-1]

    # Calculate Frequencies (we asume a sample rate of one)
    # Use fft functions for real output (-> hermitian spectrum)
    f = rfftfreq(samples)

    # Build scaling factors for all frequencies
    s_scale = f.copy()
    fmin = max(fmin, 1./samples)  # Low frequency cutoff
    ix = npsum(f < fmin)   # Index of the cutoff
    if ix and ix < len(s_scale):
        s_scale[:ix] = s_scale[ix]

    s_scale = s_scale**(-exponent/2.)

    # High-freq cutoff
    if fmax is not None:
        toKill = len(f)-npsum(f > fmax)  # Index of high-freq cutoff
        if verbose:
            print("Dropping {:d}/{:d} freqs above {:.3e} ...".format(
                npsum(f > fmax), len(f), fmax))
        if toKill <= samples:
            s_scale[toKill:] = 0

    # Calculate theoretical output standard deviation from scaling
    w = s_scale[1:].copy()
    w[-1] *= (1 + (samples % 2)) / 2.  # correct f = +-0.5
    sigma = 2 * sqrt(npsum(w**2)) / samples

    # Adjust size to generate one Fourier component per frequency
    size[-1] = len(f)

    # Add empty dimension(s) to broadcast s_scale along last
    # dimension of generated random power + phase (below)
    dims_to_add = len(size) - 1
    s_scale = s_scale[(newaxis,) * dims_to_add + (Ellipsis,)]

    # Generate scaled random power + phase
    sr = normal(scale=s_scale, size=size)
    si = normal(scale=s_scale, size=size)

    # If the signal length is even, frequencies +/- 0.5 are equal
    # so the coefficient must be real.
    if not (samples % 2):
        si[..., -1] = 0

    # Regardless of signal length, the DC component must be real
    si[..., 0] = 0

    # Combine power + corrected phase to Fourier components
    s = sr + 1J * si

    # Transform to real time series & scale to unit variance
    y = irfft(s, n=samples, axis=-1) / sigma

    return y
