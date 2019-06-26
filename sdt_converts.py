# 2019/06/25
# def sdt_converts():

# +
# FUNCTION:	n_3d(dat,ENERGY=en,ERANGE=er,EBINS=ebins,ANGLE=an,ARANGE=ar,BINS=bins)
# INPUT:
#	dat:	structure,	2d data structure filled by get_eesa_surv, get_eesa_burst, etc.
# KEYWORDS
#	ENERGY:	fltarr(2),	optional, min,max energy range for integration
#	ERANGE:	fltarr(2),	optional, min,max energy bin numbers for integration
#	EBINS:	bytarr(na),	optional, energy bins array for integration
#					0,1=exclude,include,
#					na = dat.nenergy
#	ANGLE:	fltarr(2,2),	optional, angle range for integration
#				theta min,max (0,0),(1,0) -90<theta<90
#				phi   min,max (0,1),(1,1)   0<phi<360
#	ARANGE:	fltarr(2),	optional, min,max angle bin numbers for integration
#	BINS:	bytarr(nb),	optional, angle bins array for integration
#					0,1=exclude,include,
#					nb = dat.ntheta
#	BINS:	bytarr(na,nb),	optional, energy/angle bins array for integration
#					0,1=exclude,include
# PURPOSE:
#	Returns the density, n, 1/cm^3
# NOTES:
#	Function normally called by "get_3dt" or "get_2dt" to
#	generate time series data for "tplot.pro".
#
# CREATED BY:
#	J.McFadden	95-7-27
# LAST MODIFICATION:
#	96-7-6		J.McFadden	added more keywords
# -
# def n_3d(dat2,ENERGY=en,ERANGE=er,EBINS=ebins,ANGLE=an,ARANGE=ar,BINS=bins)

import numpy as np


def moment_suite_3d(dat2, en=None, er=None, ebins=None, an=None, ar=None, bins=None):

    density = 0.
    flux3dx = 0.
    flux3dy = 0.
    flux3dz = 0.
    eflux3dx = 0.
    eflux3dy = 0.
    eflux3dz = 0.

    if not dat2.valid:
        print('Invalid Data')
        # return density

    if dat2.units_name.lower() != "eflux":
        print("Must have eflux units!")
        # return density

    # dat = conv_units(dat2,"eflux")		; Use Energy Flux
    dat = dat2
    na = dat.nenergy
    nb = dat.nbins

    ebins2 = np.ones(na, dtype=np.bool)
    # ebins2=replicate(1b,na)
    if en is not None:
        ebins2[:] = 0
        er2 = energy_to_ebin(dat, en)
        if er2[0] > er2[1]:
            er2 = np.flip(er2)
        ebins2[er2[0]:er2[1]] = 1
    if er is not None:
        ebins2[:] = 0
        er2 = er
        if er2[0] > er2[1]:
            er2 = np.flip(er2)
        ebins2[er2[0]:er2[1]] = 1

    if ebins is not None:
        ebins2 = ebins

    bins2 = np.ones(nb, dtype=np.bool)
    if an is not None:

        assert 2 < 0, "Haven't set this up yet! How to handle theta vs phi?"
        if len(an.shape) != 2:
            print('Error - angle keyword must be (2,2)')
        else:
            bins2 = angle_to_bins(dat, an)
    if ar is not None:
        bins2[:] = 0
        if ar[0] > ar[1]:
            bins2[ar[0]:nb-1] = 1
            bins2[0:ar[1]] = 1
        else:
            bins2[ar[0]:ar[1]] = 1

    if bins is not None:
        bins2 = bins

    if len(bins2.shape) != 2:
        bins2 = np.outer(ebins2, bins2)

    data = dat.data*bins2
    energy = dat.energy
    denergy = dat.denergy
    theta = np.deg2rad(dat.theta)
    phi = np.deg2rad(dat.phi)
    dtheta = np.deg2rad(dat.dtheta)
    dphi = np.deg2rad(dat.dphi)
    mass = dat.mass * 1.6e-22
    Const = np.sqrt(mass/(2.*1.6e-12))
    jConst = 1.
    jeConst = 1.6e-12

    if hasattr(dat, "domega"):
        if len(dat.domega.shape) == 1:
            domega = np.outer(np.ones(na), dat.domega)
    else:
        if len(dtheta.shape) == 1:
            dtheta = np.outer(np.ones(na), dtheta)
        if len(dphi.shape) == 1:
            dphi = np.outer(np.ones(na), dphi)
        domega = 2.*dphi*np.sin(theta)*np.sin(.5*dtheta)

    sumdata = np.sum(data*domega, axis=1)
    density = Const*np.sum((denergy[:, 0]*(energy[:, 0]**(-1.5)))*sumdata)
    # units are 1/cm^3
    # return density

    sumdataxj = np.sum(data*np.cos(phi)*domega*np.cos(theta), axis=1)
    sumdatayj = np.sum(data*np.sin(phi)*domega*np.cos(theta), axis=1)
    sumdatazj = np.sum(data*domega*np.sin(theta), axis=1)

    dnrgj = jConst*denergy*(energy**(-1))
    flux3dx = np.sum(dnrgj[:, 0]*sumdataxj)
    flux3dy = np.sum(dnrgj[:, 0]*sumdatayj)
    flux3dz = np.sum(dnrgj[:, 0]*sumdatazj)

    # Returns the flux, [Jx,Jy,Jz], 1/(cm^2-s)
    # return, [flux3dx,flux3dy,flux3dz]

    sumdataxje = np.sum(data*np.cos(phi)*domega*np.cos(theta), axis=1)
    sumdatayje = np.sum(data*np.sin(phi)*domega*np.cos(theta), axis=1)
    sumdatazje = np.sum(data*domega*np.sin(theta), axis=1)

    dnrgje = jeConst*denergy
    eflux3dx = np.sum(dnrgje[:, 0]*sumdataxje)
    eflux3dy = np.sum(dnrgje[:, 0]*sumdatayje)
    eflux3dz = np.sum(dnrgje[:, 0]*sumdatazje)

    # units are ergs/cm^2-sec
    # return, [eflux3dx,eflux3dy,eflux3dz]

    return {'density': density,
            'jx': flux3dx,
            'jy': flux3dy,
            'jz': flux3dz,
            'jex': eflux3dx,
            'jey': eflux3dy,
            'jez': eflux3dz}


def energy_to_ebin(dat, en,
                   bin2=None):

    if not dat.valid:
        print, 'Invalid Data'
        return np.nan

    if bin2 is not None:
        bin = bin2

    endim = np.size(en, 0)

    if endim == 0:
        if bin2 is None:
            bin = 0

        energy = np.squeeze(dat.energy[:, bin])

        ebin = np.argmin((energy-en)**2)

        return ebin

    else:

        ebin = np.zeros(endim, dtype=np.int64)

        if bin2 is None:
            bin = np.zeros(endim, dtype=np.int64)

        for a in range(endim):
            energy = np.squeeze(dat.energy[:, bin[a]])
            eb = np.argmin((energy-en[a])**2)
            ebin[a] = eb

        return ebin


def angle_to_bins(dat, an, ebin2=None):

    if not dat.valid:
        print('Invalid Data')
        return np.nan

    if (ebin2 is not None):
        if hasattr(ebin2, '__len__'):  # dimen(ebin2) ne 1 then begin
            print(' Error: ebin must be an integer!')
            return np.nan

    if ebin2 is not None:
        ebin = ebin2
    else:
        ebin = np.int64(dat.nenergy/2.)

    if len(np.shape(an)) <= 1:
        andim = np.size(an, 0)
        if andim != 2:
            print('Error in angle_to_bins: dimen1(an) must equal 2')
            return np.nan
        else:
            theta = np.squeeze(dat.theta[ebin, :])
            theta = 360.*(theta/360.-np.floor(theta/360.))
            th = 360.*(an/360.-np.floor(an/360.))
            if th[0] < th[1]:
                bins = (theta >= th[0]) & (theta < th[1])
            elif th[0] > th[1]:
                bins = (theta >= th[0]) | (theta < th[1])
            elif an[0] != an[1]:
                bins = bytarr(dat.nbins)
                bins[:] = 1
            else:
                bins = np.zeros(dat.nbins, dtype=np.bool)
                indmin = np.argmin(
                    np.abs(np.abs(np.abs(theta-th[0])-180.)-180.))
                bins[indmin] = 1

            return bins

    else:
        andim = dimen1(an)
        if andim != 2:
            print('Error in angle_to_bins: dimen1(an) must equal 2')
            return np.nan
        else:
            theta = np.squeeze(dat.theta[ebin, :])
            th = np.squeeze(an[:, 0])
            if (th[0] > th[1]) | (th[0] < -90.) | (th[1] > 90.):
                print('Error in angle_to_bins: -90. <= theta <= 90.')
                return np.nan

            phi = np.squeeze(dat.phi[ebin, :])
            phi = 360.*(phi/360.-np.floor(phi/360.))
            ph = np.squeeze(an[:, 1])
            ph = 360.*(ph/360.-np.floor(ph/360.))
            if th[0] != th[1]:
                if ph[0] < ph[1]:
                    bins = ((phi >= ph[0]) & (phi < ph[1])) & (
                        (theta >= th[0]) & (theta < th[1]))
                elif ph[0] > ph[1]:
                    bins = ((phi >= ph[0]) | (phi < ph[1])) & (
                        (theta >= th[0]) & (theta < th[1]))
                elif an[0, 1] != an[1, 1]:
                    bins = ((theta >= th[0]) & (theta < th[1]))
                else:
                    bins = np.zeros(dat.nbins, dtype=np.bool)
                    indmin = np.argmin(
                        np.abs(np.abs(np.abs(phi-ph[0])-180.)-180.))
                    phmin = phi[indmin]
                    bins = (phi == phmin) & (
                        (theta >= th[0]) & (theta < th[1]))

            else:
                indmin = np.argmin(np.abs(theta-th[0]))
                thmin = theta(indmin)
                if ph[0] < ph[1]:
                    bins = ((phi >= ph[0]) & (phi < ph[1])) & (theta == thmin)
                elif ph[0] > ph[1]:
                    bins = ((phi >= ph[0]) | (phi < ph[1])) & (theta == thmin)
                elif an[0, 1] != an[1, 1]:
                    bins = (theta == thmin)
                else:
                    # bins = bytarr(dat.nbins)
                    bins = np.zeros(dat.nbins, dtype=np.bool)
                    indmin = np.argmin(
                        (np.abs(np.abs(phi-ph[0])-180.)-180.)**2+(theta-th[0])**2)
                    bins[indmin] = 1

        return bins
