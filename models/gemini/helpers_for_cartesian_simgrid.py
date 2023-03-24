"""Functions that are appropriate for use with GEMNI runs based on a local
Cartesian grid (Up-East-North) with a background magnetic field 
B = B0 <0, 0, -1>

October 2022
SMH

"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator as rginterp

RE = 6370e3                     # Should be value from GEMINI, but I haven't checked

def get_ped_hall(dat,xg,
                 refalt_for_efield=500.*1000.,
                 return_enu_coordinates=False):
    """
    NOTE: I think this ONLY works for a background B-field that is constant everywhere and only in the Cartesian DOWN direction

    Also, arrays in xg are given in geomagnetic Cartesian coordinates so that order is (x_mag, y_mag, z_mag)

    If you want to do something in "model" coordinates (x1, x2, x3), you need to use model unit vectors, e.g., xg['e1']
    """

    print("Getting hall and pedersen currents and stuff in geomagnetic XYZ coords ...")
    print("# Code for plotting")
    print("#    * Histograms of Pedersen and Hall ")
    print("#    * Current density lines of Pedersen, Hall, and total current vectors")
    print("#    * Conductivity profiles")
    print("#    * Average current density profiles and cumulative height-integrated current density profiles")
    print("#    * Angles between full current vector and Hall/Pedersen currents")
    print("#    * How does E-field direction (or rather, v_i x B) change with altitude?")
    print("# can be found in in journal__20221025__trace_gemini_hall_and_pedersen_currents.py")


    # Get velocity vector and current density vector in geomagnetic XYZ coordinates
    vvec = np.transpose(dat['v1'].values[:,:,:,np.newaxis]*xg['e1']+\
                        dat['v2'].values[:,:,:,np.newaxis]*xg['e2']+\
                        dat['v3'].values[:,:,:,np.newaxis]*xg['e3'],axes=[3,0,1,2])

    Jvec = np.transpose(dat['J1'].values[:,:,:,np.newaxis]*xg['e1']+\
                        dat['J2'].values[:,:,:,np.newaxis]*xg['e2']+\
                        dat['J3'].values[:,:,:,np.newaxis]*xg['e3'],axes=[3,0,1,2])

    # bhat = np.stack([-np.ones_like(dat['v1']),
    #                  np.zeros_like(dat['v1']),
    #                  np.zeros_like(dat['v1'])])

    bhat = -np.transpose(xg['e1'],axes=[3,0,1,2])

    # Get E-field in V/m
    E = -np.cross(vvec,
                  np.abs(np.median(xg['Bmag']))*bhat,
                  axis=0)

    # But we don't want just any E-field! We want the E-field above where things are super collisional
    zind = np.argmin(np.abs(refalt_for_efield-xg['x1'][2:-2]))
    EUEN = np.stack([np.sum(np.transpose(xg['e1'],axes=[3,0,1,2])*E,axis=0),
                     np.sum(np.transpose(xg['e2'],axes=[3,0,1,2])*E,axis=0),
                     np.sum(np.transpose(xg['e3'],axes=[3,0,1,2])*E,axis=0)])
    refEs_UEN = EUEN[:,zind,:,:]
    Emapped_UEN = np.zeros_like(EUEN)
    Emapped_UEN[:,:,:,:] = refEs_UEN[:,np.newaxis,:,:]
    Emapped = np.transpose(Emapped_UEN[0,...][...,np.newaxis]*xg['e1']+\
                           Emapped_UEN[1,...][...,np.newaxis]*xg['e2']+\
                           Emapped_UEN[2,...][...,np.newaxis]*xg['e3'],
                           axes=[3,0,1,2])


    bhatcrossE = np.cross(bhat,Emapped,axis=0)

    Emag = np.sqrt(np.sum(Emapped**2,axis=0))
    ehat = Emapped/Emag
    bhatcrossE_hat = bhatcrossE/Emag

    # Pedersen stuff
    Jped = np.sum(Jvec*ehat,axis=0)
    Sigmaped = Jped/Emag
    Jpedvec = Jped*ehat

    # Hall stuff
    Jhall = np.sum(Jvec*bhatcrossE_hat,axis=0)
    Sigmahall = Jhall/Emag
    Jhallvec = Jhall*bhatcrossE_hat

    # samre dei opp
    derivs = dict(Jp=Jped,
                  Jh=Jhall,
                  Jpvec=Jpedvec,
                  Jhvec=Jhallvec,
                  SigmaP=Sigmaped,
                  SigmaH=Sigmahall,
                  E=Emapped)

    return derivs


def get_ped_hall_ENU(dat,xg,
                     refalt_for_efield=500.*1000.):
    """
    NOTE: I think this ONLY works for a background B-field that is constant everywhere and only in the Cartesian DOWN direction

    Also, arrays in xg are given in geomagnetic Cartesian coordinates so that order is (x_mag, y_mag, z_mag)

    If you want to do something in "model" coordinates (x1, x2, x3), you need to use model unit vectors, e.g., xg['e1']
    """

    print("Getting hall and pedersen currents and stuff in ENU coords ...")
    print("# Code for plotting")
    print("#    * Histograms of Pedersen and Hall ")
    print("#    * Current density lines of Pedersen, Hall, and total current vectors")
    print("#    * Conductivity profiles")
    print("#    * Average current density profiles and cumulative height-integrated current density profiles")
    print("#    * Angles between full current vector and Hall/Pedersen currents")
    print("#    * How does E-field direction (or rather, v_i x B) change with altitude?")
    print("# can be found in in journal__20221025__trace_gemini_hall_and_pedersen_currents.py")

    # Get velocity vector and current density vector in geomagnetic XYZ coordinates
    vvec = np.stack([dat['v2'].values,dat['v3'].values,dat['v1'].values])
    Jvec = np.stack([dat['J2'].values,dat['J3'].values,dat['J1'].values])

    bhat = np.zeros_like(vvec)
    bhat[2,...] = -1.

    # Get E-field in V/m
    E = -np.cross(vvec,
                  np.abs(np.median(xg['Bmag']))*bhat,
                  axis=0)

    # But we don't want just any E-field! We want the E-field above where things are super collisional
    zind = np.argmin(np.abs(refalt_for_efield-xg['x1'][2:-2]))
    refEs = E[:,zind,:,:]
    Emapped = np.zeros_like(E)
    Emapped[:,:,:,:] = refEs[:,np.newaxis,:,:]

    bhatcrossE = np.cross(bhat,Emapped,axis=0)

    Emag = np.sqrt(np.sum(Emapped**2,axis=0))
    # bhatcrossEmag = np.sqrt(np.sum(bhatcrossE**2,axis=0))
    ehat = Emapped/Emag
    bhatcrossE_hat = bhatcrossE/Emag
    # bhatcrossE_hat = bhatcrossE/bhatcrossEmag

    # Pedersen stuff
    Jped = np.sum(Jvec*ehat,axis=0)
    Sigmaped = Jped/Emag
    Jpedvec = Jped*ehat

    # Hall stuff
    Jhall = np.sum(Jvec*bhatcrossE_hat,axis=0)
    Sigmahall = Jhall/Emag
    Jhallvec = Jhall*bhatcrossE_hat

    # Calculate conductivity gradients (x1, x2, x3 correspond to U,E,N respectively)
    gradSigmaP = np.gradient(Sigmaped,xg['x1'][2:-2],xg['x2'][2:-2],xg['x3'][2:-2])
    gradSigmaH = np.gradient(Sigmahall,xg['x1'][2:-2],xg['x2'][2:-2],xg['x3'][2:-2])

    # Set in ENU order
    gradSigmaP = np.stack([gradSigmaP[1],gradSigmaP[2],gradSigmaP[0]])
    gradSigmaH = np.stack([gradSigmaH[1],gradSigmaH[2],gradSigmaH[0]])

    # samre dei opp
    derivs = dict(Jp=Jped,
                  Jh=Jhall,
                  Jpvec=Jpedvec,
                  Jhvec=Jhallvec,
                  SigmaP=Sigmaped,
                  SigmaH=Sigmahall,
                  gradSigmaP=gradSigmaP,
                  gradSigmaH=gradSigmaH,
                  E=Emapped)

    return derivs


def get_UEN_indices_of_point(point,xg):
    """
    Index by [2:-2] to drop ghost cells
    point in km, whereas xg['x1'] and friends are in m
    """
    zind = np.argmin(np.abs(point[0]*1000.-xg['x1'][2:-2]))
    xind = np.argmin(np.abs(point[1]*1000.-xg['x2'][2:-2]))
    yind = np.argmin(np.abs(point[2]*1000.-xg['x3'][2:-2]))
    
    return (zind,yind,xind)


def get_interp_x1x2x3_funcs_from_magCartesian(vec,xg):
    """Assumes vec has shape (3, x1, x2, x3) and that vec components are given in magnetospheric Cartesian coordinates
    """
    v1 = np.sum(np.transpose(xg['e1'],axes=[3,0,1,2])*vec,axis=0)
    v2 = np.sum(np.transpose(xg['e2'],axes=[3,0,1,2])*vec,axis=0)
    v3 = np.sum(np.transpose(xg['e3'],axes=[3,0,1,2])*vec,axis=0)

    v1interp = rginterp((xg['x1'][2:-2],xg['x2'][2:-2],xg['x3'][2:-2]),v1, method='linear', bounds_error=False, fill_value=np.nan)
    v2interp = rginterp((xg['x1'][2:-2],xg['x2'][2:-2],xg['x3'][2:-2]),v2, method='linear', bounds_error=False, fill_value=np.nan)
    v3interp = rginterp((xg['x1'][2:-2],xg['x2'][2:-2],xg['x3'][2:-2]),v3, method='linear', bounds_error=False, fill_value=np.nan)

    # return v1interp,v2interp,v3interp

    v1interp2 = lambda x,y,z: v1interp((x,y,z))
    v2interp2 = lambda x,y,z: v2interp((x,y,z))
    v3interp2 = lambda x,y,z: v3interp((x,y,z))

    return v1interp2,v2interp2,v3interp2



def get_interp_x1x2x3_funcs_from_UEN(vec,x1,x2,x3):
    """Assumes vec has shape (3, x1, x2, x3) and that vec components are given in (x1, x2, x3) coordinates (LOCAL UEN)
    """
    # v1 = np.sum(np.transpose(xg['e1'],axes=[3,0,1,2])*vec,axis=0)
    # v2 = np.sum(np.transpose(xg['e2'],axes=[3,0,1,2])*vec,axis=0)
    # v3 = np.sum(np.transpose(xg['e3'],axes=[3,0,1,2])*vec,axis=0)

    v1 = vec[0,:]
    v2 = vec[1,:]
    v3 = vec[2,:]

    v1interp = rginterp((x1,x2,x3),v1, method='linear', bounds_error=False, fill_value=np.nan)
    v2interp = rginterp((x1,x2,x3),v2, method='linear', bounds_error=False, fill_value=np.nan)
    v3interp = rginterp((x1,x2,x3),v3, method='linear', bounds_error=False, fill_value=np.nan)

    # return v1interp,v2interp,v3interp

    v1interp2 = lambda x1,x2,x3: v1interp((x1,x2,x3))
    v2interp2 = lambda x1,x2,x3: v2interp((x1,x2,x3))
    v3interp2 = lambda x1,x2,x3: v3interp((x1,x2,x3))

    return v1interp2,v2interp2,v3interp2


def magnetic_spherical_vector_to_local_UEN(Br,Bt,Bp,t,p,tc,pc):
    """
    Br: (spherical magnetic) radial vector component
    Bt: (spherical magnetic) colatitudinal vector component
    Bp: (spherical magnetic) azimuthal vector component
    t : theta [rad]
    p : phi   [rad]
    tc: center theta for UEN coordinate sys [rad]
    pc: center phi   for UEN coordinate sys [rad]
    """

    stc,ctc = np.sin(tc),np.cos(tc)
    spc,cpc = np.sin(pc),np.cos(pc)

    st,ct = np.sin(t),np.cos(t)
    sp,cp = np.sin(p),np.cos(p)

    a00 =  st *stc*cp *cpc + st *stc*sp *spc + ct *ctc
    a01 = -st *cp *stc     + st *sp *cpc
    a02 = -st *ctc*cp *cpc - st *ctc*sp *spc + ct *stc

    a10 =  ct *cp *stc*cpc + ct *sp *stc*spc - st *ctc
    a11 = -ct *cp *spc     + ct *sp *cpc
    a12 = -ct *cp *ctc*cpc - ct *sp *ctc*spc - st *stc 

    a20 = -sp *stc*cpc     + cp *stc*spc
    a21 =  sp *spc         + cp *cpc
    a22 =  sp *ctc*cpc     - cp *ctc*spc

    B1 = Br*a00 + Bt*a10 + Bp*a20
    B2 = Br*a01 + Bt*a11 + Bp*a21
    B3 = Br*a02 + Bt*a12 + Bp*a22

    return B1,B2,B3


def reshape_rthetaphi_to_x1x2x3(a,gridsize):
    """
    Reshaper to go from arrays organized by (r, theta, phi)
    to (x1, x2, x3), where 
    x1hat = zhat ≈ rhat
    x2hat = xhat ≈ phihat
    x3hat = yhat ≈ -thetahat
    """
    return np.flip(np.transpose(a.reshape(gridsize),axes=[0,2,1]),axis=2)


def get_cube(X1,X2,X3,cubedim: dict()):
    """
    Given simulation coordinate matrices X1, X2, X3 (which all have three dimensions) and cubedim dictionary,
    return indices, submatrices, etc., corresponding to cube dimensions

    cubedim dict is something like the following
    x1cubectr,x2cubectr,x3cubectr = 130*1e3,0.,0.  # Define box center at (130 km, 0 km, 0 km)
    dx1cube,dx2cube,dx3cube = 100*1e3,60*1e3,60*1e3  # Define box extent (±50 km, ±30 km, ±30 km)
    x1cubemin,x1cubemax = x1cubectr-dx1cube/2,x1cubectr+dx1cube/2
    x2cubemin,x2cubemax = x2cubectr-dx2cube/2,x2cubectr+dx2cube/2
    x3cubemin,x3cubemax = x3cubectr-dx3cube/2,x3cubectr+dx3cube/2
    
    cubedim = dict()            # Make "cube dimensions" dictionary
    cubedim['x1min'] = x1cubemin
    cubedim['x1max'] = x1cubemax
    cubedim['x2min'] = x2cubemin
    cubedim['x2max'] = x2cubemax
    cubedim['x3min'] = x3cubemin
    cubedim['x3max'] = x3cubemax

    """

    cube = dict()
        
    if X1.ndim == 3 and X2.ndim == 3 and X3.ndim == 3:

        cube['keepinds'] = (X1 >= cubedim['x1min']) & (X1 <= cubedim['x1max']) & \
        (X2 >= cubedim['x2min']) & (X2 <= cubedim['x2max']) & \
        (X3 >= cubedim['x3min']) & (X3 <= cubedim['x3max'])
    
        cube['x1unik'] = X1[:,0,0]              # z (egentlig r)
        cube['x2unik'] = X2[0,:,0]              # x (egentlig phi)
        cube['x3unik'] = X3[0,0,:]              # y (egentlig theta)
        
        cube['X1inds'] = (cubedim['x1min'] <= X1 ) & (X1 <= cubedim['x1max'])
        cube['X2inds'] = (cubedim['x2min'] <= X2 ) & (X2 <= cubedim['x2max'])
        cube['X3inds'] = (cubedim['x3min'] <= X3 ) & (X3 <= cubedim['x3max'])
        
        cube['x1inds'] = (cubedim['x1min'] <= cube['x1unik'] ) & (cube['x1unik'] <= cubedim['x1max'])
        cube['x2inds'] = (cubedim['x2min'] <= cube['x2unik'] ) & (cube['x2unik'] <= cubedim['x2max'])
        cube['x3inds'] = (cubedim['x3min'] <= cube['x3unik'] ) & (cube['x3unik'] <= cubedim['x3max'])
        cube['Nx1'],cube['Nx2'],cube['Nx3'] = np.sum(cube['x1inds']),np.sum(cube['x2inds']),np.sum(cube['x3inds'])
        cube['Nmeas'] = cube['Nx1']*cube['Nx2']*cube['Nx3']
        
        cube['shape'] = (cube['Nx1'],cube['Nx2'],cube['Nx3'])
    
        # NOTE: If a variable ends with "cube," it means it only includes measurements within a cube, and it is ordered by (x1,x2,x3) (or equivalently, UEN) even if raveled
        cube['X1'] = X1[cube['keepinds']].reshape(cube['Nx1'],cube['Nx2'],cube['Nx3'])
        cube['X2'] = X2[cube['keepinds']].reshape(cube['Nx1'],cube['Nx2'],cube['Nx3'])
        cube['X3'] = X3[cube['keepinds']].reshape(cube['Nx1'],cube['Nx2'],cube['Nx3'])
        
    else:
        
        cube['keepinds'] = (X1 >= cubedim['x1min']) & (X1 <= cubedim['x1max']) & \
            (X2 >= cubedim['x2min']) & (X2 <= cubedim['x2max']) & \
            (X3 >= cubedim['x3min']) & (X3 <= cubedim['x3max'])
    
        # cube['x1unik'] = X1[:,0,0]              # z (egentlig r)
        # cube['x2unik'] = X2[0,:,0]              # x (egentlig phi)
        # cube['x3unik'] = X3[0,0,:]              # y (egentlig theta)
        
        cube['X1inds'] = (cubedim['x1min'] <= X1 ) & (X1 <= cubedim['x1max'])
        cube['X2inds'] = (cubedim['x2min'] <= X2 ) & (X2 <= cubedim['x2max'])
        cube['X3inds'] = (cubedim['x3min'] <= X3 ) & (X3 <= cubedim['x3max'])
        
        cube['x1inds'] = (cubedim['x1min'] <= X1 ) & (X1 <= cubedim['x1max'])
        cube['x2inds'] = (cubedim['x2min'] <= X2 ) & (X2 <= cubedim['x2max'])
        cube['x3inds'] = (cubedim['x3min'] <= X3 ) & (X3 <= cubedim['x3max'])

        cube['Nx1'],cube['Nx2'],cube['Nx3'] = np.sum(cube['x1inds']),np.sum(cube['x2inds']),np.sum(cube['x3inds'])
        cube['Nmeas'] = np.sum(cube['x1inds'] & cube['x2inds'] & cube['x3inds'])
        
        cube['shape'] = (cube['Nmeas'],)
    
        # NOTE: If a variable ends with "cube," it means it only includes measurements within a cube, and it is ordered by (x1,x2,x3) (or equivalently, UEN) even if raveled
        cube['X1'] = X1[cube['keepinds']].reshape(cube['shape'])
        cube['X2'] = X2[cube['keepinds']].reshape(cube['shape'])
        cube['X3'] = X3[cube['keepinds']].reshape(cube['shape'])
        

    return cube

