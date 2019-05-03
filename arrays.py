# 2018/10/31
import numpy as np

r2d = 180/np.pi


def group_consecutives(vals, maxDiff=1, min_streak=None):
    """Return list of consecutive lists of numbers from vals (number list).

    av https://stackoverflow.com/questions/7352684/
    how-to-find-the-groups-of-consecutive-elements-from-an-array-in-numpy

    En dag du skulle legge til en parameter som tar hensyn til streak-størrelse
    """

    this = np.split(vals, np.where(np.diff(vals) > maxDiff)[0]+1)

    if min_streak is None:
        return this
    else:
        keep = []
        for liszt in this:
            # print(liszt.size)
            if liszt.size >= min_streak:
                keep.append(liszt)
        return keep


def find_nearest(array, value):
    array = np.asarray(array)
    if hasattr(value, '__len__'):
        if array.ndim > 1:
            print("Don't know what to do with multidimensional array and search value!")
            return -1
        idx = (np.abs(array.reshape(array.size, 1) - value)).argmin(axis=0)
    else:
        idx = (np.abs(array - value)).argmin()
    return idx


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
