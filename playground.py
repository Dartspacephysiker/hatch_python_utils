#2020/01/22
import numpy as np

def latmap(xs, coef, epsilon, x0, reflect = True):
    """
    xs      : latitudes (not colat)
    coef    : determines epsilon dependence
    epsilon : sw coupling
    x0      : low latitude boundary
    """
    xs = np.deg2rad(xs)
    sign = np.sign(xs)
    xs = np.abs(xs)

    # define constants
    A0 = np.e/4 * (np.pi - 2 * x0) # max amplitude
    A  = 1 - np.exp(-coef * epsilon)  # scaling factor
    a  = (np.pi / 2 - x0)**2

    x0 = x0 * np.pi/180 # convert to radian

    retvals = []
    # for i in range(len(epsilon)):
    xx = (xs - x0)
    f  = np.zeros_like(xs)
    iii = xs > x0
    f[iii] = -180/np.pi * A[iii] * A0 * (1 - a/xx[iii]**2) * np.exp(-a/xx[iii]**2)
    y =  f + xs* 180/np.pi

    # now find the values at x - possibly with reflection about y = x:
    # if reflect: # with reflection
    #     retvals.append(np.interp(x[i], y, xs) * 180 /np.pi)
    # else: # without reflection:
    #     x[i] = x[i] * np.pi / 180
    #     retvals.append(np.interp(x[i], xs, y))
    # return np.array(retvals) * sign 	
    return y

def latmapbegin(x, coef, epsilon, x0 = 75.):
    """
    x       : latitude (not colat)
    coef    : determines epsilon dependence
    epsilon : sw coupling
    x0      : low latitude boundary
    """

    x0 = x0 * np.pi/180 # conver to radian
    x  = x  * np.pi/180 # conver to radian

    assert (x >= 0).all()

    # define constants
    A0 = np.e/4 * (np.pi - 2 * x0) # max amplitude
    A  = 1 - np.exp(-coef * epsilon)  # scaling factor
    a  = (np.pi / 2 - x0)**2
    xx = (x - x0)
    f  = np.zeros_like(x)
    iii = x > -x0
    f[iii] = -180/np.pi * A[iii] * A0 * (1 - a/xx[iii]**2) * np.exp(-a/xx[iii]**2)
    return f + x* 180/np.pi 
