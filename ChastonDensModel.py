# 2019/03/13
# def ChastonDensModel():

import numpy as np


def dipolefield(alt):

    Bo = 5.40115e-5  # T at surface
    RE = 6370.e3
    fromcentre = alt+RE
    return Bo*(RE/fromcentre)**3.


def alfven_speed_mlt(alt, mlt,
                     doReturnFreqs=False):
    """
    Ripped off Chris Chaston's IDL Alfvén speed routine (including his implementations of dayside and nightside density models)

    Inputs
    ------
    alt : Geodetic altitude   (km      )
    mlt : Magnetic local time (0--24 hr)

    Output:
    if doReturnFreqs:
        Alfvén speed in km/s as well as omega_(p,hydro), omega_(p,oxy), omega_(c,hydro), omega_(c,oxy)
    """
    c = 3.0e5
    me = 9.1e-31  # electron mass in kg
    mp = 1.67e-27  # proton mass in kg
    mox = 16.*1.67e-27
    qe = 1.6e-19  # electron charge in C
    qp = 1.6e-19

    # density_func_night='density_model_night_dipole'
    # density_func_day='density_model_day_dipole'
    # dheight=call_function(density_func_night,alt)
    dheight = density_model_night_dipole(alt)
    bheight = dipolefield(alt*1000.0)
    hplasmaf = np.sqrt(
        (np.sum(dheight[1:])*1.0e6*(1.0*qe)**2)/(8.85e-12*1.0*mp))
    oplasmaf = np.sqrt((dheight[0]*1.0e6*(1.0*qe)**2)/(8.85e-12*16.0*mp))
    hcycf = (1.0*qe*bheight)/(1.0*mp)
    ocycf = (1.0*qe*bheight)/(16.0*mp)
    va = c*1./np.sqrt((hplasmaf/hcycf)**2+(oplasmaf/ocycf)**2)
    va_km_per_sec_night = va/np.sqrt(1+va**2/c**2)

    # dheight=call_function(density_func_day,alt)
    dheight = density_model_day_dipole(alt)
    hplasmaf = np.sqrt(
        (np.sum(dheight[1:])*1.0e6*(1.0*qe)**2)/(8.85e-12*1.0*mp))
    oplasmaf = np.sqrt((dheight[0]*1.0e6*(1.0*qe)**2)/(8.85e-12*16.0*mp))
    hcycf = (1.0*qe*bheight)/(1.0*mp)
    ocycf = (1.0*qe*bheight)/(16.0*mp)
    va = c*1./np.sqrt((hplasmaf/hcycf)**2+(oplasmaf/ocycf)**2)
    va_km_per_sec_day = va/np.sqrt(1+va**2/c**2)

    va_mlt = va_km_per_sec_night * \
        np.abs(mlt-12.)/12.+va_km_per_sec_day*(1.-np.abs(mlt-12.)/12.)

    if doReturnFreqs:
        return va_mlt, hplasmaf, oplasmaf, hcycf, ocycf
    else:
        return va_mlt

# ;;2016/05/15
# ;;Here are the variables that differ between dayside and nightside
# ;;Night
# ;;  nh=1000.0
# ;;  no_FF=8.0e4
# ;;  scale_upper_h=70.0
# ;;  scale_upper_o_e=70.0
# ;;  scale_upper_o_f=70.0
# ;;  n_out_o=6000.0
# ;;  scale_upper_out_o=400.
#
# ;;day
# ;;  nh=2000.0
# ;;  no_FF=1.0e5
# ;;  scale_upper_h=140.0
# ;;  scale_upper_o_e=140.0
# ;;  scale_upper_o_f=210.0
# ;;  n_out_o=4000.0
# ;;  scale_upper_out_o=600.
#


def density_model_night_dipole(r):
    # used this function for dayside density profile using Chapman profiles
    Re = 6370.0  # km
    altitude = r  # (r-Re)
    nh = 1000.0
    no_E = 1.0e4
    no_FF = 8.0e4
    E_alt = 110.0
    F_alt = 200.0
    scale_lower_h = 30.0
    scale_upper_h = 70.0

    scale_lower_o_e = 5.0
    scale_upper_o_e = 70.0

    scale_lower_o_f = 10.0
    scale_upper_o_f = 70.0

    scale_lower_ps = 20000.0
    scale_upper_ps = 20000.0

    power_law_h = -2.0
    power_law_o = -2.0
    power_law_ps = -2

    n_out_h = 1000.0
    out_alt_h = 500.0
    scale_lower_out_h = 1000.
    scale_upper_out_h = 800.

    n_out_o = 6000.0
    out_alt_o = 500.0
    scale_lower_out_o = 1000.
    scale_upper_out_o = 400.

    cluster_alt = 5.4*Re
    nh_f = 1.0

    # HYDROGEN
    # Hydrogen ionospheric
    if altitude >= F_alt:
        nh_ionos = nh*np.exp(1-np.abs(altitude-F_alt)/scale_upper_h -
                             np.exp(-np.abs(altitude-F_alt)/scale_upper_h))
    else:
        nh_ionos = nh*np.exp(1-np.abs(altitude-F_alt)/scale_lower_h -
                             np.exp(-np.abs(altitude-F_alt)/scale_lower_h))

    # hydrogen outflow
    if altitude >= out_alt_h:
        nh_ionos_out = n_out_h*np.exp(1-np.abs(altitude-out_alt_h)/scale_upper_out_h -
                                      np.exp(-np.abs(altitude-out_alt_h)/scale_upper_out_h))
    else:
        nh_ionos_out = n_out_h*np.exp(1-np.abs(altitude-out_alt_h)/scale_lower_out_h -
                                      np.exp(-np.abs(altitude-out_alt_h)/scale_lower_out_h))

    nh_ionos = nh_ionos+nh_ionos_out

    # hydrogen plasmasheet
    # nh_ps=1.0*tanh((altitude-100.)/(2.*6371.2))
    # if r GE sp_r-6370.0 then nh_ps=nh_f*(r/sp_r)**power_law_ps else nh_ps=nh_f*((sp_r-6370.0)/sp_r)**(power_law_ps)
    # if altitude GE cluster_alt then begin
    #   nh_ps=nh_f*np.exp(1-np.abs(altitude-cluster_alt)/scale_upper_ps-np.exp(-np.abs(altitude-cluster_alt)/scale_upper_ps))
    # endif else begin
    #   nh_ps=nh_f*np.exp(1-np.abs(altitude-cluster_alt)/scale_lower_ps-np.exp(-np.abs(altitude-cluster_alt)/scale_lower_ps))

    nh_ps = 1.0

    # OXYGEN
    # E-region
    if altitude >= E_alt:
        no_ionos = no_E*np.exp(1-np.abs(altitude-E_alt)/scale_upper_o_e -
                               np.exp(-np.abs(altitude-E_alt)/scale_upper_o_e))
    else:
        no_ionos = no_E*np.exp(1-np.abs(altitude-E_alt)/scale_lower_o_e -
                               np.exp(-np.abs(altitude-E_alt)/scale_lower_o_e))

    # F-region
    if altitude >= F_alt:
        no_ionos = no_ionos+no_FF * \
            np.exp(1-np.abs(altitude-F_alt)/scale_upper_o_f -
                   np.exp(-np.abs(altitude-F_alt)/scale_upper_o_f))
    else:
        no_ionos = no_ionos+no_FF * \
            np.exp(1-np.abs(altitude-F_alt)/scale_lower_o_f -
                   np.exp(-np.abs(altitude-F_alt)/scale_lower_o_f))

    # ionospheric outflow
    if altitude >= out_alt_o:
        no_ionos_out = n_out_o*np.exp(1-np.abs(altitude-out_alt_o)/scale_upper_out_o -
                                      np.exp(-np.abs(altitude-out_alt_o)/scale_upper_out_o))
    else:
        no_ionos_out = n_out_o*np.exp(1-np.abs(altitude-out_alt_o)/scale_lower_out_o -
                                      np.exp(-np.abs(altitude-out_alt_o)/scale_lower_out_o))

    no_ionos = no_ionos_out+no_ionos

    return np.array([no_ionos, nh_ps, nh_ionos])

#
# ##2016/05/15
# ##Here are the variables that differ between dayside and nightside
# ##Night
# ##  nh=1000.0
# ##  no_FF=8.0e4
# ##  scale_upper_h=70.0
# ##  scale_upper_o_e=70.0
# ##  scale_upper_o_f=70.0
# ##  n_out_o=6000.0
# ##  scale_upper_out_o=400.
#
# ##day
# ##  nh=2000.0
# ##  no_FF=1.0e5
# ##  scale_upper_h=140.0
# ##  scale_upper_o_e=140.0
# ##  scale_upper_o_f=210.0
# ##  n_out_o=4000.0
# ##  scale_upper_out_o=600.
#


def density_model_day_dipole(r):
    # used this function for dayside density profile using Chapman profiles

    Re = 6370.0  # km
    altitude = r  # (r-Re)
    nh = 2000.0
    no_E = 1.0e4
    no_FF = 1.0e5
    E_alt = 110.0
    F_alt = 200.0
    scale_lower_h = 30.0
    scale_upper_h = 140.0

    scale_lower_o_e = 5.0
    scale_upper_o_e = 140.0

    scale_lower_o_f = 10.0
    scale_upper_o_f = 210.0

    scale_lower_ps = 20000.0
    scale_upper_ps = 20000.0

    power_law_h = -2.0
    power_law_o = -2.0
    power_law_ps = -2

    n_out_h = 1000.0
    out_alt_h = 500.0
    scale_lower_out_h = 1000.
    scale_upper_out_h = 800.

    n_out_o = 4000.0
    out_alt_o = 500.0
    scale_lower_out_o = 1000.
    scale_upper_out_o = 600.

    cluster_alt = 5.4*Re
    nh_f = 1.0

    # HYDROGEN
    # Hydrogen ionospheric
    if altitude >= F_alt:
        nh_ionos = nh*np.exp(1-np.abs(altitude-F_alt)/scale_upper_h -
                             np.exp(-np.abs(altitude-F_alt)/scale_upper_h))
    else:
        nh_ionos = nh*np.exp(1-np.abs(altitude-F_alt)/scale_lower_h -
                             np.exp(-np.abs(altitude-F_alt)/scale_lower_h))

    # hydrogen outflow

    if altitude >= out_alt_h:
        nh_ionos_out = n_out_h*np.exp(1-np.abs(altitude-out_alt_h)/scale_upper_out_h -
                                      np.exp(-np.abs(altitude-out_alt_h)/scale_upper_out_h))
    else:
        nh_ionos_out = n_out_h*np.exp(1-np.abs(altitude-out_alt_h)/scale_lower_out_h -
                                      np.exp(-np.abs(altitude-out_alt_h)/scale_lower_out_h))

    nh_ionos = nh_ionos+nh_ionos_out

    # hydrogen plasmasheet
    # nh_ps=1.0*tanh((altitude-100.)/(2.*6371.2))
    # if r >= sp_r-6370.0 then nh_ps=nh_f*(r/sp_r)**power_law_ps else nh_ps=nh_f*((sp_r-6370.0)/sp_r)**(power_law_ps)
    # if altitude >= cluster_alt:
    #   nh_ps=nh_f*np.exp(1-np.abs(altitude-cluster_alt)/scale_upper_ps-np.exp(-np.abs(altitude-cluster_alt)/scale_upper_ps))
    # else:
    #   nh_ps=nh_f*np.exp(1-np.abs(altitude-cluster_alt)/scale_lower_ps-np.exp(-np.abs(altitude-cluster_alt)/scale_lower_ps))

    nh_ps = 1.0

    # OXYGEN
    # E-region
    if altitude >= E_alt:
        no_ionos = no_E*np.exp(1-np.abs(altitude-E_alt)/scale_upper_o_e -
                               np.exp(-np.abs(altitude-E_alt)/scale_upper_o_e))
    else:
        no_ionos = no_E*np.exp(1-np.abs(altitude-E_alt)/scale_lower_o_e -
                               np.exp(-np.abs(altitude-E_alt)/scale_lower_o_e))

    # print,no_ionos
    # F-region
    if altitude >= F_alt:
        no_ionos = no_ionos+no_FF * \
            np.exp(1-np.abs(altitude-F_alt)/scale_upper_o_f -
                   np.exp(-np.abs(altitude-F_alt)/scale_upper_o_f))
    else:
        no_ionos = no_ionos+no_FF * \
            np.exp(1-np.abs(altitude-F_alt)/scale_lower_o_f -
                   np.exp(-np.abs(altitude-F_alt)/scale_lower_o_f))

    # ionospheric ooutflow
    if altitude >= out_alt_o:
        no_ionos_out = n_out_o*np.exp(1-np.abs(altitude-out_alt_o)/scale_upper_out_o -
                                      np.exp(-np.abs(altitude-out_alt_o)/scale_upper_out_o))
    else:
        no_ionos_out = n_out_o*np.exp(1-np.abs(altitude-out_alt_o)/scale_lower_out_o -
                                      np.exp(-np.abs(altitude-out_alt_o)/scale_lower_out_o))

    no_ionos = no_ionos+no_ionos_out

    return np.array([no_ionos, nh_ps, nh_ionos])
