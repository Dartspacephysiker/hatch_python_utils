# 2019/11/05
import datetime
from operator import attrgetter


def doy_single(dt):
    tt = dt.timetuple()
    doy = tt.tm_yday+tt.tm_hour/24.+tt.tm_min/1440.+tt.tm_sec/86400.
    return doy


def doy_from_tuple(timetuple):
    """
    Tuple skal være ('tm_yday','tm_hour','tm_min','tm_sec')
    """
    return timetuple[0]+timetuple[1]/24.+timetuple[2]/1440.+timetuple[3]/86400.


def datetime_to_doy(dts):
    tider = attrgetter('tm_yday', 'tm_hour', 'tm_min', 'tm_sec')
    return [doy_from_tuple(tider(dt.timetuple())) for dt in dts]


def jd_to_datetime(jds):
    noonJan1200UT__JD = 2451545.
    tdeltas = jds-noonJan1200UT__JD
    return [datetime.datetime(2000, 1, 1, 12)+datetime.timedelta(days=tdiff) for tdiff in tdeltas]
