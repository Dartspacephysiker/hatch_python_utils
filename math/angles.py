def in_angle_interval(x, a, b,degrees=True):
    """
    Find out if x is in angle interval [a,b]
    From https://stackoverflow.com/questions/42133572/create-a-circular-list-by-using-a-range-of-angles-python
    """
    if degrees:
        maxAngle = 360
    else:
        maxAngle = 2*np.pi

    return (x - a) % maxAngle <= (b - a) % maxAngle
