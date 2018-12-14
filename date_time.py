# 2018/10/31
# Velg denne dag
import random
import time


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
