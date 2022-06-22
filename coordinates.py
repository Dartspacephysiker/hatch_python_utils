"""
As originally employed in journal__20220613__does_rotating_into_centered_pole_system_remove_multivaluednessQQ_plus_test_for_multivaluedness.py
"""

import numpy as np


def sph2cart_coords(theta,phi,r=1,deg=True):

    if deg:
        tr = np.deg2rad(theta)
        pr = np.deg2rad(phi)
    else:
        tr, pr = theta, phi

    x = r * np.sin(tr) * cos(pr)
    y = r * np.sin(tr) * sin(pr)
    z = r * np.cos(tr)

    return (x,y,z)


def cart2sph_coords(x,y,z,deg=True):

    r = np.sqrt(x**2+y**2+z**2)
    H = np.sqrt(x**2+y**2)
    theta = np.arctan2(H,z)
    phi = np.arctan2(y,x)

    if deg:
        theta = np.rad2deg(theta)
        phi = np.rad2deg(phi)

    return (r,theta,phi)
