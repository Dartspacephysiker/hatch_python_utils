def sphDist(lat1, mlt1, lat2, mlt2):
    """
    Great-circle distance
    lat1: Latitude of 1st point in degrees
    lat2: Latitude of 2nd point in degrees
    mlt1: Magnetic local time of 1st point in hours (so between 0 and 24)
    mlt2: Magnetic local time of 2nd point in hours
    """
    lat1R = np.deg2rad(lat1)
    mlt1R = np.deg2rad(mlt1*15.)
    lat2R = np.deg2rad(lat2)
    mlt2R = np.deg2rad(mlt2*15.)

    return np.rad2deg(np.arccos(np.sin(lat1R)*np.sin(lat2R)+np.cos(lat1R)*np.cos(lat2R)*np.cos(np.abs(mlt1R-mlt2R))))
