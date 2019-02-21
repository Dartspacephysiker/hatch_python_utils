# 2018/10/31
import numpy as np


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
