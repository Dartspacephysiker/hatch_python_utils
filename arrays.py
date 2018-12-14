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
