import numpy as np
from scipy.interpolate import BSpline,splev
import scipy.integrate as integrate


def get_bsplineobj_list(t,k,extrapolate=False):
    """Get only basis functions (represented as BSpline objects) that are non-zero
    somewhere within the internal knot interval (t[k], t[-k-1]). 

    """
    bsplineobjs = []
    bsplinefuncs = []

    internalknots = t[k:-k]
    nfuncs = len(t)-k-1
    for i in range(nfuncs):
        weights = np.zeros(nfuncs)
        weights[i] = 1

        # Too simple, leaves NaNs lying about
        # bsplineobjs.append(BSpline(t,weights,k,extrapolate=extrapolate))

        bsplobj = BSpline(t,weights,k,extrapolate=extrapolate)
        def bsplfunc(z,
                    zmin=internalknots[0],
                    zmax=internalknots[-1],
                    func=bsplobj):
    
            result = np.zeros_like(z)
            inrange = (z >= zmin) & (z <= zmax)
            # boverange = (z >= zmax)
            result[inrange] = func(z[inrange])
            # result[boverange] = maxval-zeroval
            return result

        bsplineobjs.append(bsplobj)
        bsplinefuncs.append(bsplfunc)

    return bsplineobjs,bsplinefuncs


def get_bsplineobjintegral_list(t,k,extrapolate=False):
    """Get only basis functions (represented as BSpline objects) that are non-zero
    somewhere within the internal knot interval (t[k], t[-k-1]). 

    """
    bsplineobjs = []

    internalknots = t[k:-k]

    nfuncs = len(t)-k-1
    for i in range(nfuncs):
        weights = np.zeros(nfuncs)
        weights[i] = 1
        extrapolate = False

        bsplobj = BSpline(t,weights,k,extrapolate=extrapolate)
        def antifunc(z,
                     zmin=internalknots[0],
                     zmax=internalknots[-1],
                     antiderfunc=bsplobj.antiderivative(),
                     zeroval=bsplobj.antiderivative()(internalknots[0]),
                     maxval=bsplobj.antiderivative()(internalknots[-1])):
    
            result = np.zeros_like(z)
            inrange = (z >= zmin) & (z <= zmax)
            boverange = (z >= zmax)
            result[inrange] = antiderfunc(z[inrange])-zeroval
            result[boverange] = maxval-zeroval
            return result
    
        bsplineobjs.append(antifunc)

    return bsplineobjs


def get_bsplinefunc_list(t,k,ext=1):
    """Get only basis functions that are non-zero somewhere within the internal knot
    interval (t[k], t[-k-1]). The sum of these basis functions (when they all
    have coeff 1) is 1.

    """
    bsplinefuncs = []

    nfuncs = len(t)-k-1
    for i in range(nfuncs):
        weights = np.zeros(nfuncs)
        weights[i] = 1
        bsplinefuncs.append(lambda h, t=t, weights=weights, k=k, ext=ext: splev(h, (t,weights,k), der=0, ext=ext))

    return bsplinefuncs
        

def get_bsplinefuncintegral_list(t,k,integlowbound,ext=1):
    """Get only integrals of basis functions that are non-zero somewhere within the internal knot
    interval (t[k], t[-k-1]). The sum of these basis functions (when they all
    have coeff 1) is 1.

    """
    bsplineintegfuncs = []

    nfuncs = len(t)-k-1
    for i in range(nfuncs):
        weights = np.zeros(nfuncs)
        weights[i] = 1
        tmpfunc = lambda h, t=t, weights=weights, k=k, ext=ext: splev(h, (t,weights,k), der=0, ext=ext)

        def tmpintegfunc(hmax, integlowbound=integlowbound, tmpfunc=tmpfunc):
            if hmax < integlowbound:
                return 0.

            return integrate.quad(tmpfunc, integlowbound, hmax, epsabs=1e-5)[0]

        # tmpintegfunc = lambda h, integlowbound=integlowbound: integrate.quad(tmpfunc, integlowbound, h,
        #                   epsabs=1e-4)[0]

        bsplineintegfuncs.append(tmpintegfunc)

    return bsplineintegfuncs
        


