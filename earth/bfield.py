import numpy as np
from datetime import datetime
import ppigrf
from pysymmetry import geodesy
import apexpy

def geoc2apex_dB_vector(gclat,glon,r_km,apex_date,
                        B_NEC_vec,
                        igrf_NEC_vec=None,
                        refh=110,
                        verbose=True):
    """
    INPUTS
    ======
    gclat        : Geocentric latitude [deg]
    glon         : Longitude           [deg]
    r_km         : radius              [km]
    apex_date    : date-time like   
    B_NEC_vec    : Vector of B-field measurements, NEC coords [nT]

    OPTIONAL
    ========
    igrf_NEC_vec : Vector of IGRF main field, NEC coords [nT]
    refh         : Reference height for Modified Apex coordinate system [km]


    NOTE: B_NEC_vec and igrf_NEC_vec must have shape (3,N), where N is the number of measurements.

    RETURNS
    =======
    dBe1         : Field perturbation in Modified Apex e1-direction (approx. E-W)
    dBe2         : Field perturbation in Modified Apex e2-direction (approx. N-S)

    EXAMPLE
    =======
    #Swarm A B-field measurement from 2017-04-18 00:00:00.019

    import numpy as np
    from datetime import datetime

    gclat,glon,r_km = -83.04955,-36.1731,6827.577
    B_NEC_vec = np.array([1.411483e4, 7.766765e2, -3.739525e4])
    B_NEC_vec = B_NEC_vec[:,np.newaxis] #Force shape to be (3,1)
    apex_date = datetime(2017,4,18,0,0,0,19000)

    dBe1, dBe2 = geoc2apex_dB_vector(gclat,glon,r_km,apex_date,B_NEC_vec)
    """

    # 
    # assert len(B_NEC_vec.shape) == len(igrf_NEC_vec.shape) == 2,"B_NEC_vec and igrf_NEC_vec must both have shape (3,N)!"
    # assert B_NEC_vec.shape[0] == igrf_NEC_vec.shape[0] == 3,"B_NEC_vec and igrf_NEC_vec must both have shape (3,N)!"
    assert (len(B_NEC_vec.shape) == 2) and (B_NEC_vec.shape[0] == 3),"B_NEC_vec must have shape (3,N)!"

    # Get IGRF if not provided by user
    if igrf_NEC_vec is not None:
        assert (len(igrf_NEC_vec.shape) == 2) and (igrf_NEC_vec.shape[0] == 3),"igrf_NEC_vec must have shape (3,N)!"
        igrfN,igrfE,igrfC = igrf_NEC_vec
    else:
        if verbose:
            print(f"Getting IGRF from ppigrf for {apex_date} ...")

        igrfr, igrftheta, igrfphi = ppigrf.igrf_gc(r_km, 90.-gclat, glon, apex_date) 
        igrfN,igrfE,igrfC = -igrftheta,igrfphi,-igrfr

    B_N,B_E,B_C = B_NEC_vec

    # Get ΔB in geodetic coordinates
    gdlat, alt, DX, DZ = geodesy.geoc2geod(90.-gclat,
                                           r_km,
                                           -(B_N-igrfN),  # colatitudinal component
                                           -(B_C-igrfC))  # radial component
        
    # Create apex object, get basevectors
    a = apexpy.Apex(apex_date, refh=refh)

    f1, f2, f3, g1, g2, g3, d1, d2, d3, e1, e2, e3 = a.basevectors_apex(
        gdlat, glon, alt, coords='geo')

    dB_E = B_E-igrfE
    dB_NGeod = DX
    dB_UGeod = -DZ
        
    # Convert ΔB in geodetic coordinates to Modified-Apex(refh) coordinates
    dBe1 = dB_E*e1[0,] \
        + dB_NGeod*e1[1,] \
        + dB_UGeod*e1[2,]
        
    dBe2 = dB_E*e2[0,] \
        + dB_NGeod*e2[1,] \
        + dB_UGeod*e2[2,]
        
    # dBpar = (B_E*igrfENorm \
    #          + B_NGeod*igrfNGeodNorm \
    #          + B_UGeod*igrfUGeodNorm) - igrfMag

    return dBe1, dBe2


def field_geometry_lonlat(lon,lat,date,h_km=0):
    """
    Calculate inclination, declination, and dip latitude 
    for given location(s) using IGRF.

    Example:
    --------
    inc, dec, diplat = field_geometry_lonlat()

    Parameters
    ----------
    lon : array
        longitude [deg], postiive east, of IGRF calculation
    lat : array
        geodetic latitude [deg] of IGRF calculation
    date : date(s)
        one or more dates to evaluate IGRF coefficients

    Keywords
    ----------
    h : array
        height [km] above ellipsoid for IGRF calculation

    Returns
    -------
    inc    : array
        Geomagnetic field inclination [degrees]
    dec    : array
        Geomagnetic field declination [degrees]
    diplat : array
        Geomagnetic dip latitude      [degrees]
    

    Definitions
    -----------
    Inclination : Angle made with horizontal by Earth's magnetic field lines
    Declination : Angle between geodetic north and magnetic north
    Dip latitude: Latitude of given dip angle (related to magnetic latitude?)
    S. M. Hatch
    Mar 2022
    
    """

    Be, Bn, Bu = ppigrf.igrf(lon, lat, h_km, date) # returns east, north, up

    Z = -Bu
    H = np.sqrt(Be**2+Bn**2)

    inc    = np.arctan2(Z,H)
    dec    = np.arccos(Bn/H)
    diplat = np.arctan2(Z,2*H)
    
    inc,dec,diplat = map(lambda x: np.ravel(np.rad2deg(x)), [inc,dec,diplat])

    return inc, dec, diplat

if __name__ == '__main__':
    
    loc = 'Layton, UT'
    lat,lon = 41.066, -111.961
    EW = 'E' if lon > 0 else 'W'
    date = datetime(2012,8,5,0)

    inc,dec,diplat = field_geometry_lonlat(lon,lat,date)
    inc,dec,diplat = map(lambda x: x[0], [inc,dec,diplat])

    print(f"{loc} ({lat:.2f} N, {np.abs(lon):.2f} {EW})")
    print(date)
    print("==============================")
    print("")
    print("{:20s} : {:.2f}".format("Inclination  [deg]",inc))
    print("{:20s} : {:.2f}".format("Declination  [deg]",dec))
    print("{:20s} : {:.2f}".format("Dip latitude [deg]",diplat))

    print("")
    print("Declination from NCEI geomagnetic calculator (https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml):")
    print('"2012-08-05	12.10° E changing by  0.12° W per year"')
