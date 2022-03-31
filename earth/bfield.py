import numpy as np
from datetime import datetime
import ppigrf

def field_geometry_lonlat(lon,lat,date,h_km=0):
    """
    Calculate inclination, declination, and dip latitude 
    for given location(s) using IGRF.

    Example:
    --------
    inc, dec = field_geometry_lonlat()

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
    Mar 2020
    
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
