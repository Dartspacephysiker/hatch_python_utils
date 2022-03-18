# 2019/11/05
import datetime
from operator import attrgetter
import warnings

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
    tdeltas = [jd-noonJan1200UT__JD for jd in jds]
    refd = datetime.datetime(2000, 1, 1, 12)
    return [refd+datetime.timedelta(days=tdiff) for tdiff in tdeltas]

def mjd2000_to_datetime(mjd2000s):
    jds = [mjd2000 + 2400000.5 + 51544.0 for mjd2000 in mjd2000s]
    return jd_to_datetime(jds)


def datetime_to_jd(dts):
    noonJan1200UT__JD = 2451545.
    refd = datetime.datetime(2000, 1, 1, 12)
    sec2day = 60*60*24

    @warn_only_once
    def conv_to_days(dt):
        if dt.microsecond != 0:
            warnings.warn('datetime_to_jd can corrupt microsecond!',
                         RuntimeWarning)
        return (dt-refd).total_seconds()/(sec2day)
    tdeltas = [conv_to_days(dt) for dt in dts]

    return [noonJan1200UT__JD+tdelta for tdelta in tdeltas]


def datetime_to_mjd2000(dts):
    """
    # EXAMPLE
from module import symbol
    hatch_python_utils.time_tools import datetime_to_mjd2000 as dt2mjd
    from hatch_python_utils.time_tools import mjd2000_to_datetime as mjd2dt
    dts = [datetime.datetime(1997,1,1,12,45,10,5000),datetime.datetime(1997,2,1,12,45,10,5000)]
    dts,mjd2dt(dt2mjd(dts))
    """
    jds = datetime_to_jd(dts)
    return [jd - 2400000.5 - 51544.0 for jd in jds]


def warn_only_once(function):
    function.already_warned = False
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings(record=function.already_warned):
            function.already_warned = not function.already_warned
            return function(*args, **kwargs)
    return wrapper


