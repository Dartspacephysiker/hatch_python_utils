# 2019/11/05
import datetime
from operator import attrgetter


def get_UT(dts):
    try:
        UT = dts.hour+dts.minute/60+dts.second/3600
    except:
        try:
            UT = [dtid.hour+dtid.minute/60+dtid.second/3600 for dtid in dts]
        except:
            print("Couldn't convert dts to UT!")
            return

    return UT

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


def datetime_to_yearfrac(dts):
    """
    Convert list of datetime-like objects to year fractions.
    For example, 2001-06-30 becomes 2001.4945364574537

    2020-04-17    SMH 
    """
    return [dt.year+doy_single(dt)/doy_single(datetime.datetime(dt.year,12,31,23,59)) for dt in dts]


def jd_to_datetime(jds):
    noonJan1200UT__JD = 2451545.
    tdeltas = jds-noonJan1200UT__JD
    return [datetime.datetime(2000, 1, 1, 12)+datetime.timedelta(days=tdiff) for tdiff in tdeltas]

def mjd2000_to_datetime(mjd2000s):
    jds = mjd2000s + 2400000.5 + 51544.0
    return jd_to_datetime(jds)
