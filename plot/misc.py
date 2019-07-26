# 2018/12/28
import numpy as np
import matplotlib.pyplot as plt


def mystep(x, y, ax=None, giveXandY=False, **kwargs):
    # x = np.array(x)
    # y = np.array(y)
    # x2 = x*0.
    X = np.c_[x[:-1], x[1:], x[1:]]
    # Y = np.c_[y[:-1],y[:-1],np.zeros_like(x2[:-1])*np.nan]
    Y = np.c_[y[:-1], y[:-1], np.zeros(x.shape[0]-1)*np.nan]
    if not ax:
        ax = plt.gca()
    # return ax.plot(X.flatten(), Y.flatten(), **kwargs)
    if giveXandY:
        return (X, Y, ax.plot(X.flatten(), Y.flatten(), **kwargs))
    else:
        return ax.plot(X.flatten(), Y.flatten(), **kwargs)
