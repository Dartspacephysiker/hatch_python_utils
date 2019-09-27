# 2018/10/31
# Velg denne dag
import random
import time
from datetime import datetime


def toYearFraction(date):
    """
    From https://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years
    """
    def sinceEpoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    s = sinceEpoch

    year = date.year
    startOfThisYear = datetime(year=year, month=1, day=1)
    startOfNextYear = datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def strTimeProp(start, end, format, prop):
    """Get a time at a proportion of a range of two formatted times.

    av https://stackoverflow.com/questions/553303/generate-a-random-date-between-two-other-dates

    start and end should be strings specifying times formated in the
    given format (strftime-style), giving an interval [start, end].
    prop specifies how a proportion of the interval to be taken after
    start.  The returned time will be in the specified format.

    Example:
    =======
    start = '2015-01-01T00:00:00'
    stop  = '2015-12-31T11:59:59'
    fmt = '%Y-%m-%dT%H:%M:%S'
    this = np.random.uniform()
    print(this)
    strTimeProp(start,stop,fmt,this)

    """

    stime = time.mktime(time.strptime(start, format))
    etime = time.mktime(time.strptime(end, format))

    ptime = stime + prop * (etime - stime)

    return time.strftime(format, time.localtime(ptime))


def cdf_epoch_to_utc(cdf_epochs):
    return (cdf_epochs-62167219200000.0000)/1000.0


def cdf_epoch_to_datetime(cdf_epochs):
    return [datetime.utcfromtimestamp(utc) for utc in cdf_epoch_to_utc(cdf_epochs)]
