import numpy as np

def median(x,weights):
    """
    # Based on Wikipedia article (https://en.wikipedia.org/wiki/Weighted_median)

    # Test 
    x,weights = np.arange(1,6),np.array([0.15,0.1,0.2,0.3,0.25])
    print(np.median(x),weighted.median(x,weights))

    # Test special case where weights land you right in between 0.5 (see Wikipedia article on weighted median)
    x,weights = np.arange(1,5),np.array([0.25,0.25,0.25,0.25])
    print(np.median(x),weighted.median(x,weights))
    """
    return quantile(x,weights,0.5)

def quantile(x,weights,q):

    assert (q >= 0) and (q <= 1)

    sortie = np.argsort(x)
    x = x[sortie]
    weights = weights[sortie] / np.sum(weights)

    weightcum = np.cumsum(weights)

    ind = np.argmin(np.abs(weightcum-q))
    if np.isclose(weightcum[ind],q) and np.isclose(1-weightcum[ind+1]+weights[ind+1],1-q):
        return np.mean([x[ind],x[ind+1]])
    else:
        ind += int(weightcum[ind] < q)  # Add 1 if we're below the requested quantile

    return x[ind]
	

