import os
import numpy as np
import h5py

from hatch_python_utils.models.gemini.helpers_for_cartesian_simgrid import reshape_rthetaphi_to_x1x2x3

RE = 6370e3

def HHMMSS_to_s(h,m,s):
    return h*3600+m*60+s

def s_to_HHMMSS(x):
    h = x//3600
    m = (x%3600)//60
    s = x%60

    return h,m,s


def load_bfield_stuff(direc,bfieldpoints_fname,bfield_fname,thetactr,phictr,iscartesiangrid=True,
                      simple_decimate_factor=None):

    bfpoints_fullfn = os.path.join(direc, bfieldpoints_fname)
    bf_fullfn = os.path.join(direc, bfield_fname)

    if iscartesiangrid:
        print("Assuming GEMINI grid for this run is local Cartesian (UEN)")

    if simple_decimate_factor is None:
        decimator = slice(None,None,None)
    else:
        decimator = slice(None,None,simple_decimate_factor)

    mag = {}
    # magfull = {}
    keys = ('r','theta','phi','gridsize','lpoints')
    with h5py.File(bfpoints_fullfn) as maggie:
        gridsize = maggie['gridsize'][:]

        # Should we reshape these guys?
        do_reshape = len(gridsize) == 3
        if do_reshape:          # Passed first check for reshaping, now second check
            do_reshape = ( gridsize[1] != -1 ) and ( gridsize[2] != -1 )

            if simple_decimate_factor is not None:
                print("Can't 'simple decimate' B-field measurements since they are in a grid (or at least I haven't figured that out)!")

        for key in keys:
            if key in ('r','theta','phi'):
                # mag[key] = maggie[key][:]
                # magfull[key] = maggie[key][:]

                if do_reshape:
                    mag[key] = reshape_rthetaphi_to_x1x2x3(maggie[key][:],gridsize)
                else:
                    mag[key] = maggie[key][:][decimator]
    
            else:
                mag[key] = np.array(maggie[key])
    
    mag['h_km'] = (mag['r']-RE)/1000
    
    keys = ('Br','Btheta','Bphi')
    with h5py.File(bf_fullfn) as maggie:
        for key in keys:
            # mag[key] = maggie['magfields'][key][:]
            # magfull[key] = maggie['magfields'][key][:]

            if do_reshape:
                mag[key] = reshape_rthetaphi_to_x1x2x3(maggie['magfields'][key][:],gridsize)
            else:
                mag[key] = maggie['magfields'][key][:][decimator]

    # Reorder gridsize after running reshape_rthetaphi_to_x1x2x3
    if do_reshape:
        gridsize = gridsize[[0,2,1]]
    
    # Where we have B-field perturbations in geomagnetic ECEF
    mag['x'] = mag['r']*np.sin(mag['theta'])*np.cos(mag['phi'])
    mag['y'] = mag['r']*np.sin(mag['theta'])*np.sin(mag['phi'])
    mag['z'] = mag['r']*np.cos(mag['theta'])
    
    # Get B-field perturbation locations in x1, x2, x3 coords, doing the same warping crap as Matt
    if iscartesiangrid:
        mag['gamma2'] = thetactr-mag['theta']
        mag['gamma1'] = mag['phi']-phictr
        mag['x1'] = mag['r']-RE
        mag['x2'] = RE*np.sin(thetactr)*mag['gamma1']
        mag['x3'] = RE*mag['gamma2']
    
    return mag

