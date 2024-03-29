# 2018/10/31
import numpy as np

r2d = 180/np.pi


def isSorted(array):
    return np.array_equal(array, np.sort(array))


def isSorted_alt(x, key=lambda x: x):
    """
    Just another way to do it
    """
    return all([key(x[i]) <= key(x[i + 1]) for i in range(len(x) - 1)])


def group_consecutives(vals, maxDiff=1,
                       min_streak=None,
                       do_absDiff=False,
                       print_summary=False,
                       print__maxVal=None):
    """Return list of consecutive lists of numbers from vals (number list).

    av https://stackoverflow.com/questions/7352684/
    how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy

    En dag du skulle legge til en parameter som tar hensyn til streak-størrelse
    """

    # assert np.issubdtype(maxDiff,np.integer),"maxDiff must be of type integer!"
    if maxDiff <= 0:
        print("maxDiff must be >= 0!")
        return np.array([], dtype=np.int64)

    if do_absDiff:
        this = np.split(vals, np.where(np.abs(np.diff(vals)) > maxDiff)[0]+1)
    else:
        this = np.split(vals, np.where(np.diff(vals) > maxDiff)[0]+1)

    # if min_streak is None:
        # return this
    # else:
    if min_streak is not None:
        keep = []
        for liszt in this:
            # print(liszt.size)
            if liszt.size >= min_streak:
                keep.append(liszt)

        # return keep
        this = keep

    if print_summary:
        titleStr = "{:2s} {:{width}s}   {:{width}s} - {:{width}s} (/{:{width}s})"
        rowStr = titleStr.replace("s", "d")

        # rowStr = "{:2d} {:{width}d}   {:{width}d} - {:{width}d} (/{:{width}d})"

        nDigs = str(np.int64(np.log10(np.max(vals))))

        if print__maxVal is None:
            print__maxVal = np.max(vals)
        print(titleStr.format("i",
                              "N",
                              "strt",
                              "stop",
                              "tot",
                              width=nDigs))
        for i, dudG in enumerate(this):
            print(rowStr.format(i,
                                len(dudG),
                                dudG[0],
                                dudG[-1],
                                print__maxVal,
                                width=nDigs))

    return this


def find_nearest(array, value,OLDWAY=False):
    """
    Bruk heller get_closest, which is multidim
    """

    if not OLDWAY:

        return get_closest(value,array)

    else:

        array = np.asarray(array)
        if hasattr(value, '__len__'):
            if array.ndim > 1:
                print("Don't know what to do with multidimensional array and search value!")
                return -1
            idx = (np.abs(array.reshape(array.size, 1) - value)).argmin(axis=0)
        else:
            idx = (np.abs(array - value)).argmin()
        return idx


def get_closest(x,bins):
    """
    Here is the NWO in terms of speed. Use this'n instead of find_nearest.
    """
    assert np.allclose(np.sort(bins)-bins,0)  # If this isn't true, algorithm doesn't work

    nbins = len(bins)
    res = np.digitize(x,bins)-1

    res[res < 0] = 0
    res[res > (nbins-1)] = nbins-1
    # altres = np.clip(res+1,0,nbins-1)
    # res = res + np.argmin([np.abs(x-bins[res]),np.abs(x-bins[altres])],axis=0)
    res = res + np.argmin([np.abs(x-bins[res]),np.abs(x-bins[np.clip(res+1,0,nbins-1)])],axis=0)
    
    return res


def get_streaks(array, minNStreak=None, maxGap=1):

    start_i = -1
    stop_i = -1
    streakLens = -1

    if not np.array_equal(array, np.sort(array)):
        print("Bogus array! You need to sort")
        return start_i, stop_i, streakLens

    aSize = array.size
    if not isinstance(aSize, int):
        if len(aSize) > 1:
            print("Multi-d array! Don't know what to do ...")
            return start_i, stop_i, streakLens

    diff = np.diff(array)
    splitPoints = np.where(diff >= maxGap)[0]
    start_i = np.append(np.array([0]), splitPoints+1)
    stop_i = np.append(splitPoints, np.array([array.size-1]))

    streakLens = stop_i - start_i

    if minNStreak is not None:
        keep_ii = np.where(streakLens >= minNStreak)[0]

        start_i = start_i[keep_ii]
        stop_i = stop_i[keep_ii]
        streakLens = streakLens[keep_ii]

    return start_i, stop_i, streakLens

# def group_consecutives(vals, stepsize=1, min_streak=None):
#     """Return list of consecutive lists of numbers from vals (number list).

#     av https://stackoverflow.com/questions/7352684/
#     how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy

#     En dag du skulle legge til en parameter som tar hensyn til streak-størrelse
#     """

#     this = np.split(vals, np.where(np.diff(vals) != stepsize)[0]+1)

#     if min_streak is None:
#         return this
#     else:
#         keep = []
#         for liszt in this:
#             # print(liszt.size)
#             if liszt.size >= min_streak:
#                 keep.append(liszt)
#         return keep


# # def group_consecutives(vals, step=1):
# #     run = []
# #     result = [run]
# #     expect = None
# #     for v in vals:
# #         if (v == expect) or (expect is None):
# #             run.append(v)
# #         else:
# #             run = [v]
# #             result.append(run)
# #         expect = v + step
# #     return result

# getIndices = lambda seq, vals: np.array([(np.abs(seq-val)).argmin() for val in vals])

def getIndices(seq, vals):
    return np.array([(np.abs(seq-val)).argmin() for val in vals])


def value_locate(seq, vals):
    """
    Et svar på IDL'S VALUE_LOCATE
    """
    return getIndices(seq, vals)


def windowmaker(winsize, winslide, totSec,
                sample_rate_Hz=1,
                verbose=False):

    nslides = np.ceil((totSec) / (winslide/sample_rate_Hz))+1

    startSecs = []
    stopSecs = []
    frontDupes = []
    backDupes = []

    winslideSec = winslide/sample_rate_Hz
    halfwinSec = winsize/2/sample_rate_Hz

    for i in range(int(nslides)):
        # startSec = i*winslide - halfwin
        startSec = i*winslideSec - halfwinSec

        toStartSec = 0
        toFrontOfStopSec = 0

        if startSec < 0:
            toFrontOfStopSec = -startSec
            startSec = 0
        # stopSec = i*winslide + halfwin + toFrontOfStopSec
        stopSec = i*winslideSec + halfwinSec + toFrontOfStopSec

        if stopSec > totSec:
            assert toFrontOfStopSec == 0, "Bad!"
            # print("Bad!")
            # break
            # else:
            toStartSec = totSec-stopSec
            stopSec = totSec
            startSec += toStartSec

        startSecs.append(startSec)
        stopSecs.append(stopSec)

        if verbose:
            print("{:3d} {:7.2f} {:7.2f} {:7.2f}".format(
                i, startSec, stopSec, startSec+stopSec))

    startSecs, stopSecs = np.array(startSecs), np.array(stopSecs)

    dupeFront = startSecs == startSecs[0]
    dupeBack = stopSecs == stopSecs[-1]
    dupeFront[0] = False
    dupeBack[-1] = False

    startSecs = startSecs[(~dupeFront) & (~dupeBack)]
    stopSecs = stopSecs[(~dupeFront) & (~dupeBack)]

    return startSecs, stopSecs


def bruteNMin(this, nMinsWant,
              quiet=True):
    """
    Accomplish exactly the same thing with this:
    np.argpartition(this,range(nMinsWant))

    Just try it!

    from hatch_python_utils import arrays as hArr
    this = np.random.uniform(0,1,10)
    nMinsWant = 4
    print(hArr.bruteNMin(this,nMinsWant))
    print(np.argpartition(this,range(nMinsWant)))
    """
    count = 0
    mask = np.ones(len(this), dtype=bool)
    minInds = np.zeros(nMinsWant, dtype=np.int64)
    while count < nMinsWant:
        tmpInd = this[mask].argmin()

        realInd = np.where((np.cumsum(mask)-1) == tmpInd)[0][0]

        if not quiet:
            print("tmpInd, realInd: {:4d}, {:4d}".format(tmpInd, realInd))

        minInds[count] = realInd
        mask[realInd] = False

        count += 1

    return minInds


def get_rollangle_for_min_anglerange(x,
                                     dangle=5,
                                     minangle=0,
                                     maxangle=180,
                                     DEBUG=False):
    """
    When using scipy.stat's binned_statistic with an independent variable that is periodic
    (say, angles over 0–180 deg), it sometimes happens that more bins are created than are 
    actually necessary. This function finds the angle by which to roll the x axis that 
    will minimize the spread of the independent variable

    EXAMPLE
    =======
    import numpy as np
    from hatch_python_utils.arrays import get_rollangle_for_min_anglerange
    x = np.array([5,20,40,150,160,180])
    print(x.max()-x.min())
    best = get_rollangle_for_min_anglerange(x)
    xtmp = x - best
    xtmp[xtmp < 0] += 180
    print(xtmp.max()-xtmp.min())
    """

    if np.where((x < minangle) | (x > maxangle))[0].size > 0:
        assert 2<0,"Gonna break your algorithm! x is incompatible with given minangle/maxangle!"

    diffs = np.arange(minangle,maxangle+dangle/2,dangle,
                      dtype=np.int64 if np.isclose(np.int64(dangle),dangle) else np.float64)
    spread = np.zeros(len(diffs))*np.nan    

    for idd,diff in enumerate(diffs):
        tmp = x - diff
        tmp[tmp < minangle] += maxangle
        if DEBUG:
            print(idd,diff,tmp)
        spread[idd] = tmp.max() - tmp.min()
    
    winna = np.argmin(spread)
    return diffs[winna]

def carVec_to_sphVec(theta, phi,
                     vecx, vecy, vecz,
                     deg=False):

    if deg == False:
        conv = 1.
    else:
        conv = r2d

    rComp = vecx * np.sin(theta * conv) * np.cos(phi * conv) + \
        vecy * np.sin(theta * conv) * np.sin(phi * conv) + \
        vecz * np.cos(theta * conv)

    tComp = vecx * np.cos(theta * conv) * np.cos(phi * conv) + \
        vecy * np.cos(theta * conv) * np.sin(phi * conv) + \
        vecz * (-1.) * np.sin(theta * conv)

    pComp = vecx * (-1.) * np.sin(phi * conv) + \
        vecy * np.cos(phi * conv)

    return np.vstack((rComp, tComp, pComp))

    # return np.vstack((r * np.sin(theta * conv * conv) * np.cos(phi * conv),
    #                   r * np.sin(theta * conv * conv) * np.sin(phi * conv),
    #                   r * np.cos(theta * conv)))


def sphVec_to_carVec(theta, phi,
                     vecR, vecT, vecP,
                     deg=False):

    if deg == False:
        conv = 1.
    else:
        conv = r2d

    xComp = vecR * np.sin(theta * conv) * np.cos(phi * conv) + \
        vecT * np.cos(theta * conv) * np.cos(phi * conv) - \
        vecP * np.sin(phi * conv)

    yComp = vecR * np.sin(theta * conv) * np.sin(phi * conv) + \
        vecT * np.cos(theta * conv) * np.sin(phi * conv) + \
        vecP * np.cos(phi * conv)

    zComp = vecR * np.cos(theta * conv) - \
        vecT * np.sin(theta * conv)

    return np.vstack((xComp, yComp, zComp))

    # return np.vstack((r * np.sin(theta * conv * conv) * np.cos(phi * conv),
    #                   r * np.sin(theta * conv * conv) * np.sin(phi * conv),
    #                   r * np.cos(theta * conv)))
