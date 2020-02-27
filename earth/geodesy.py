import numpy as np

def get_max_sza(z,
                R=6371.):
    """
    z is altitude in km
    # R = 6371. # Earth radius (default)
    """

    # assert not hasattr(R,'__iter__')
    zIsArray = hasattr(z,'size')

    if not zIsArray:
        z = np.array([z])

    h = z + R
    max_sza = 180.-np.rad2deg(np.arctan2(np.sqrt(1-(R/h)**2.), (h/R-R/h)))

    # if zIsArray:
        
    max_sza[np.isclose(z,0)] = 90.
        
    max_sza[z < 0] = np.nan
        
    # else:
        
    #     if np.isclose(z,0):
    #         max_sza = 90.
    #     elif (z < 0):
    #         max_sza = np.nan

    if zIsArray:
        return max_sza
    else:
        return max_sza[0]

def get_h2d_bin_areas(minlats, maxlats, minlons, maxlons, haversine=True,
                      rearth=6.370949e3,
                      altitude=0,
                      spherical_rectangle=True,
                      do_extra_width_calc=True):
    # +
    # NAME:                        GET_H2D_BIN_AREAS
    #
    # PURPOSE:                     And now you wanna know the area, too?!?!?!
    #                              Well, you can have it--in km^2
    # CATEGORY:
    #
    # MODIFICATION HISTORY:      2016/03/12 Barnebarn
    #
    # Notáte bien—se dan las áreas en kilometros cuadrados
    # -

    ## BINEDGE1=Binedge1, BINEDGE2=Binedge2

    # breakpoint()

    # rearth = 6.370949e3+altitude
    # print("rearth:", rearth)

    # print("Want more precision? Try https://math.stackexchange.com/questions/1205927/how-to-calculate-the-area-covered-by-any-spherical-rectangle")
    areas = latlong_pair_area(
        minlons, minlats, maxlons, maxlats, rearth=rearth, haversine=haversine,
        altitude=altitude,
        spherical_rectangle=spherical_rectangle,
        do_extra_width_calc=do_extra_width_calc)  # use radius at 100 km

    return areas


def latlong_pair_area(lon1, lat1,
                      lon2, lat2,
                      rearth=6.370949e3,
                      altitude=0,
                      haversine=True,
                      spherical_rectangle=True,
                      do_extra_width_calc=True):
    """
    Get area (in km^2) of a rectangle on a sphere defined by two lat/lon pairs
    """
    # 2016/03/12 Spence needs this
    # Example: Surf. Area of Utah in two squares
    # ##:Upper square:
    # lonLat1  = [-114.043579, 41.951320]
    # lonLat2  = [-111.055298, 41.021355]
    # upArea   = latlong_pair_area(lonLat1[0],lonLat1[1],lonLat2[0],lonLat2[1])
    #
    # ##:Lower square:
    # lonLat1  = [-114.043579, 41.021355]
    # lonLat2  = [-109.046173, 37.009133]
    # downArea = latlong_pair_area(lonLat1[0],lonLat1[1],lonLat2[0],lonLat2[1])
    #
    # uparea+downarea = 212565.86 km^2
    # actual area of Utah = 84899 sq. mi. * (1.60934 km/mi)^2 = 219886.3 sq. km.
    # Error = 3.33% (7320 sq. km.). Not bad!

    # rearth += altitude

    height = geo_dist(lon1, lat1,
                      lon1, lat2,
                      rearth=rearth,
                      altitude=altitude,
                      haversine=haversine)

    if do_extra_width_calc:
        width1 = geo_dist(lon1, lat1,
                          lon2, lat1,
                          rearth=rearth,
                          altitude=altitude,
                          haversine=haversine)
        width2 = geo_dist(lon1, lat2,
                          lon2, lat2,
                          rearth=rearth,
                          altitude=altitude,
                          haversine=haversine)

        # width = np.max([width1, width2])

        # width = np.mean([width1, width2])
        width = np.mean(np.vstack([width1,width2]),axis=0)
    else:
        width = geo_dist(lon1, lat1,
                         lon2, lat1,
                         rearth=rearth,
                         altitude=altitude,
                         haversine=haversine)

    if spherical_rectangle:

        return 4*rearth**2.*np.arcsin(np.tan(width/rearth/2.)*np.tan(height/rearth/2.))
    else:
        # print("Vanlig")
        return height*width


def geo_dist(lon1, lat1,
             lon2, lat2,
             deg=True,
             rearth=6.370949e3,
             altitude=0,
             haversine=True):  # default to earth radius in km
    # MODIFICATION HISTORY:
    # 	Written by:	Daithi A. Stone (stoned@uvic.ca), 2000-06-29.
    #	Modified:	DAS, 2000-07-06 (removed LENGTH.pro, added
    #			DIMENSION.pro).
    #	Modified:	DAS, 2000-07-24 (added Degrad constant).
    #	Modified:	DAS, 2000-08-17 (coverted from Cosine to Haversine
    #			Haversine formula, added SINGLE keyword).
    #	Modified:	DAS, 2002-08-12 (complies with CONSTANTS.pro revision)
    #   Modified:       SMH, 2016-03-12 Made my own stuff
    #   Modified:       SMH, 2019-11-18 From IDL code

    rearth += altitude

    if deg:
        lat1, lat2, lon1, lon2 = map(np.deg2rad, (lat1, lat2, lon1, lon2))

    # Difference coordinates
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    # breakpoint()

    if haversine:

        # ***********************************************************************
        # Haversine Formula

        # Main calculation
        a = (np.sin(dlat / 2.))**2 + np.cos(lat1) * \
            np.cos(lat2) * (np.sin(dlon / 2.))**2
        a = np.sqrt(a)

        # A fix if a>1
        id = np.where(a > 1)[0]

        if id.size > 0:
            a[id] = 1

        # Convert to distance
        dist = rearth * 2 * np.arcsin(a)

    else:
        # ***********************************************************************
        # Cosine Formula (not used, but I thought I would leave it in)

        # Convert to Cartesian coordinates
        x1 = np.cos(lon1) * np.cos(lat1)
        y1 = np.sin(lon1) * np.cos(lat1)
        z1 = np.sin(lat1)
        x2 = np.cos(lon2) * np.cos(lat2)
        y2 = np.sin(lon2) * np.cos(lat2)
        z2 = np.sin(lat2)

        # Direction cosine
        dx = x1*x2+y1*y2+z1*z2

        # A fix if the |cosine| > 1
        id = np.where(dx > 1)[0]
        if id.size > 0:
            dx[id] = 1.

        #Output (distance)
        dist = rearth*np.arccos(dx)

    return dist
