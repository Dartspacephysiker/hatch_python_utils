"""Continuous CECS utils

NOTE: Here we use a Cartesian coordinate system x,y,z such that z is UPWARD
(and, say, x is east and y is north). This is DIFFERENT from the Cartesian
coordinate system that Vanhamäki (2007) and Vanhamäki, Viljanen, and Amm (2005)
use to define CECS. For them, x is north, y is east, and z is DOWN.

Why does it matter? Because if you compare expressions for J and B in this
implementation with their expressions, you will notice sign differences.

Also, some conventions (established by references, not me):
•A positive-amplitude CECS CF system corresponds to a DOWNWARD field-aligned
 current and a radially OUTWARD current sheet.
•A positive-amplitude CECS DF system corresponds to a current that circulates in
 the CLOCKWISE direction when seen from above (i.e., the -phihat direction).


REFERENCES
==========

Vanhamäki, H. (2007) ‘Theoretical modeling of ionospheric electrodynamics
including induction effects’, Finnish Meteorological Institute Contributions.

Vanhamäki, H., Viljanen, A. and Amm, O. (2005) ‘Induction effects on ionospheric
electric and magnetic fields’, Annales Geophysicae, 23(5), pp. 1735–1746. doi:
10.5194/angeo-23-1735-2005.

Spencer M. Hatch
October 2022
"""

from scipy.interpolate import BSpline
from scipy.interpolate import splev
from scipy.special import logit
import scipy.integrate as integrate
from scipy.interpolate import RegularGridInterpolator as rginterp

import numpy as np

from hatch_python_utils.math.bspline import get_bsplineobj_list, get_bsplineobjintegral_list

d2r = np.pi/180
MU0 = 4 * np.pi * 1e-7



class BsplineCECS(object):

    def __init__(self,x_cecs,y_cecs,knotvector,
                 polyorder=3,
                 current_type='divergence_free',
                 grid_type='regular',
                 exclude_laplacian_edge=None):
        """
        x_cecs: x coordinates of continuous CECS (
        y_cecs: y coordinates of continuous CECS 

        grid_type: Either "regular" or "irregular". 

                   If "regular," x_cecs and y_cecs are EITHER 2D arrays with the 
                   same shape, or 1D arrays that will be broadcast to create a 
                   2D array.

        exclude_laplacian_edge: Whether or not to use outer edge of cCECS grid in
                        evaluating model B-field and J vectors. What value this
                        keyword takes on depends on whether grid_type is regular
                        or irregular.

                        If grid_type is "regular", exclude_laplacian_edge may be set to
                        True.  In this case the edge rows and columns will be
                        taken to be the laplacian cCECS.

                        If grid_type is "irregular", the user must provide an
                        array of booleans to indicate which cCECS comprise the
                        laplacian edge. 

        current_type: string, optional
            The type of CECS function. This must be either 
            'divergence_free' (default): divergence-free basis functions
            'curl_free': curl-free basis functions

        """

        assert grid_type in ['regular','irregular'],"'grid_type' must be either 'regular' or 'irregular'! See doc string"

        self.grid_type = grid_type

        # Determine how to set up grid; also handle laplacian edge, if any
        if self.grid_type == 'regular':
            if x_cecs.ndim == 1:
                assert y_cecs.ndim == 1,"If x_cecs is 1D, y_cecs must also be 1D!"
                assert (len(np.unique(x_cecs)) == len(x_cecs)) and \
                    (len(np.unique(y_cecs)) == len(y_cecs)),"For regular grid and 1D x_cecs and y_cecs, elements of x_cecs and y_cecs must be unique!"
                
                x_cecs,y_cecs = np.meshgrid(x_cecs,y_cecs,indexing='ij')

            if x_cecs.ndim == 2:
                assert x_cecs.shape == y_cecs.shape,"When providing 2D array for cCECS x and y locations, x_cecs and y_cecs must have same shape!"

            # Note, for following code to work, it is required that x_cecs.ndim == 2
            self.laplacian_edge = np.zeros_like(x_cecs,dtype=bool)
            if exclude_laplacian_edge is not None:
                assert exclude_laplacian_edge in [True,False],"For regular grid, exclude_laplacian_edge must be True or False"
                if exclude_laplacian_edge:
                    self.laplacian_edge[[0,-1],:] = True
                    self.laplacian_edge[:,[0,-1]] = True
            self.laplacian_edge = self.laplacian_edge.ravel()

        elif self.grid_type == 'irregular':
                
            assert x_cecs.size == y_cecs.size,f"For irregular grid, x_cecs and y_cecs must have the same size! You've provided "\
                f"x_cecs and y_cecs of size {x_cecs.size} and {y_cecs.size}, respectively"

            self.laplacian_edge = np.zeros_like(x_cecs,dtype=bool)
            if exclude_laplacian_edge is not None:
                assert exclude_laplacian_edge.size == x_cecs.size,f"laplacian_edge array must have same size as x_cecs and y_cecs! You've provided "\
                    f"laplacian_edge and x_cecs of size {laplacian_edge.size} and {x_cecs.size}, respectively"
                assert np.issubdtype(exclude_laplacian_edge.dtype,bool),"For irregular grid, laplacian_edge must have dtype=bool"
                self.laplacian_edge = exclude_laplacian_edge


        # CECS locations
        self.xc = x_cecs.ravel()
        self.yc = y_cecs.ravel()

        self.Ncoeffs = len(knotvector)-polyorder-1         # Number of B-spline coefficients
        self.NCECS = len(self.xc)                   # Number of cCECS locations
        self.Nlaplacian_edge = np.sum(self.laplacian_edge)  # Number of cCECS comprising "laplacian edge"

        self.knotvec = knotvector
        self.polyorder = polyorder
        self.current_type = current_type
        self.Bspl_func_peakheights = self.knotvec[self.polyorder-1:-self.polyorder+1]
        self.internalknotinterval = (self.knotvec[self.polyorder],self.knotvec[-self.polyorder-1])
        self.zmin = self.internalknotinterval[0]
        self.zmax = self.internalknotinterval[1]
        # Create basis functions
        # Since we assume every continuous CECS uses the same knot vector and polyorder, we only need to do this once
        # self.bsplinefuncs = get_bsplinefunc_list(self.knotvec,self.polyorder,ext=1)
        # self.bsplineintegfuncs = get_bsplinefuncintegral_list(self.knotvec,self.polyorder,
        #                                                       self.zmin,
        #                                                       ext=1)
        self.bsplineobjs,self.bsplinefuncs = get_bsplineobj_list(self.knotvec,self.polyorder,extrapolate=False)
        # self.bsplinefuncs = get_bsplineobj_list(self.knotvec,self.polyorder,extrapolate=False)
        self.bsplineintegfuncs = get_bsplineobjintegral_list(self.knotvec,self.polyorder,
                                                             extrapolate=False)

        # Should have size self.NCECS*self.Ncoeffs
        # self.m = np.zeros((self.NCECS, self.Ncoeffs))
        self.m = np.zeros(self.NCECS*self.Ncoeffs)

    def get_J_G_matrices(self,x,y,z,
                         constant = 1./(2*np.pi), 
    ):
        """ Calculate matrices Gx and Gy that relate CECS amplitudes to current density 
            vector components.
    
        Parameters
        ----------
        x,y,z: array-like
            Array of Cartesian coordinate of evaluation points [m]
            Flattened arrays must all have the same size
        {x,y,z}_cecs: array-like
            Array of CECS pole Cartesian coordinates [m]
            Flattened arrays must all have the same size
        constant: float, optional
            The CECS functions are scaled by the factor 1/(2pi), which is
            the default value of 'constant'.
    
        Returns
        -------
        If current_type is 'divergence_free' or 'curl_free':
        Gx: 2D array
            2D array with shape (x.size, x_cecs.size), relating CECS amplitudes
            m to the x-direction current densities at (x, y, z) via 'jx = Gx.dot(m)'
        Gy: 2D array
            2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
            m to the y-direction current densities at (x, y, z) via 'jy = Gy.dot(m)'
        # Gz: 2D array
        #     2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
        #     m to the z-direction current densities at (x, y, z) via 'jz = Gz.dot(m)'
        If current_type is 'potential' or 'scalar':
        G: 2D array
            2D array with shape (lat.size, lat_cecs.size), relating amplitudes m
            to scalar field magnitude at (lat, lon) via 'z = G.dot(m)'
    
        """
        
        if self.current_type not in ['divergence_free','curl_free']:
            assert 2<0,"types 'potential' and 'scalar' not yet implemented!"
    
        xp, yp = get_horizontal_prime_coordinates(x,y,self.xc,self.yc)
    
        rho = get_rho(x,y,self.xc,self.yc)
    
        # Get B-spline function values for each z
        # bsplvals = np.stack([self.bsplinefuncs[i](z) for i in range(self.Ncoeffs)]).T
        bsplvals = np.vstack([self.bsplinefuncs[i](z) for i in range(self.Ncoeffs)]).T

        if self.current_type == 'divergence_free':
            unit_vec = (-1.)*get_phihat(x,y,self.xc,self.yc)
        elif self.current_type == 'curl_free':
            unit_vec = get_rhohat(x,y,self.xc,self.yc)
        elif self.current_type in ['potential', 'scalar']:
            unit_vec = 1
        else:
            raise Exception('type must be "divergence_free", "curl_free", "potential", or "scalar"')
    
        # get the scalar part of Amm's divergence-free CECS:    
        if self.current_type in ['divergence_free', 'curl_free']:
            coeff = constant / rho[:,:,np.newaxis]*bsplvals[:,np.newaxis,:]
    
            # G matrices
            Gx = coeff * unit_vec[:, :, np.newaxis, 0]
            Gy = coeff * unit_vec[:, :, np.newaxis, 1]

            Gx = Gx.reshape((Gx.shape[0],np.prod(Gx.shape[1:]))) 
            Gy = Gy.reshape((Gy.shape[0],np.prod(Gy.shape[1:]))) 

            return Gx, Gy
        else: # self.current_type is 'potential' or 'scalar'
            # NOT SURE WHAT THESE SHOULD BE
            # if self.current_type == 'potential':
            #     return -2*constant*np.log(np.sin(theta/2)).T
            # elif self.current_type == 'scalar':
            #     return    constant      / np.tan(theta/2).T
    
            return np.nan, np.nan
    

    def get_B_G_matrices(self,x,y,z,
                         constant = MU0/(4.*np.pi), 
                         df__only_z=False,
                         round_rho_dec_place=None,
                         verbose = False,
                         debug = False,
    ):
        """Calculate matrices Gx and Gy that relate CECS amplitudes to B-field
            perturbation vector components.
    
        Parameters
        ----------
        x,y,z: array-like
            Array of Cartesian coordinate of evaluation points [m]
            Flattened arrays must all have the same size
        {x,y,z}_cecs: array-like
            Array of CECS pole Cartesian coordinates [m]
            Flattened arrays must all have the same size
        constant: float, optional
            The CECS functions are scaled by the factor 1/(2pi), which is
            the default value of 'constant'.
    
        df__only_z: Only calculate z component of G matrix for divergence-free CECS
        round_rho_dec_place : This can save you LOTS of time! 

        Returns
        -------
        If current_type is 'divergence_free' or 'curl_free':
        Gx: 2D array
            2D array with shape (x.size, x_cecs.size), relating CECS amplitudes
            m to the x-direction B-field perturbations at (x, y, z) via 'Bx = Gx.dot(m)'
        Gy: 2D array
            2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
            m to the y-direction B-field perturbations at (x, y, z) via 'By = Gy.dot(m)'
        # Gz: 2D array
        #     2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
        #     m to the z-direction B-field perturbations at (x, y, z) via 'Bz = Gz.dot(m)'
        If current_type is 'potential' or 'scalar':
        G: 2D array
            2D array with shape (lat.size, lat_cecs.size), relating amplitudes m
            to scalar field magnitude at (lat, lon) via 'z = G.dot(m)'

        NOTES
        =====
        In this method I use the following notation for the various indices
        k -> Number of measurements (i.e., length of x, y, or z)
        i -> Number of CECS lines (i.e., self.NCECS)
        j -> Number of B-spline coefficients per CECS line (i.e., self.Ncoeffs)
        """
        
        if self.current_type not in ['divergence_free','curl_free']:
            assert 2<0,"types 'potential' and 'scalar' not yet implemented!"
    
        xp, yp = get_horizontal_prime_coordinates(x,y,self.xc,self.yc)
    
        rho = get_rho(x,y,self.xc,self.yc)  # Shape is (Nmeas, self.NCECS)

        Nmeas = len(z)
        NCECS = self.NCECS
        Ncoeff = self.Ncoeffs

        if self.current_type == 'divergence_free':
            rhohat = get_rhohat(x,y,self.xc,self.yc)
            zhat = np.zeros_like(rhohat)
            zhat[...,2] = 1.

            # OLD, don't consider degenerateness of integrals
            # k indexes measurements
            # i indexes (x,y) positions of BsplineCECS lines
            # j indexes B-spline functions
            
            # Cij = np.zeros((Nmeas,NCECS,Ncoeff))
            # Dij = np.zeros((Nmeas,NCECS,Ncoeff))
            # zeval = np.linspace(*self.internalknotinterval,501)
            # for k in range(Nmeas):
            #     if (k % 100) == 0:
            #         print("{:04d}/{:04d}".format(k,Nmeas))
            #     for i in range(NCECS):
            #         for j in range(Ncoeff):

            #             Cij[k,i,j] = Bdf_rhointeg_trapz(self.bsplinefuncs[j],rho[k,i],z[k],zeval)


            #             Dij[k,i,j] = Bdf_zinteg_trapz(self.bsplinefuncs[j],rho[k,i],z[k],zeval)


            # NEW, consider degenerateness of z and rho
            # k indexes unique z
            # i indexes unique rho
            # j indexes B-spline functions
            zunik,z_k = np.unique(z,return_inverse=True)
            if round_rho_dec_place is not None:
                rhounik,rho_i = np.unique(np.round(rho,round_rho_dec_place),return_inverse=True)
            else:
                rhounik,rho_i = np.unique(rho,return_inverse=True)
            rho_i = rho_i.reshape(rho.shape)

            Nzunik = len(zunik)
            Nrhounik = len(rhounik)

            Cunik = np.zeros((Nzunik,Nrhounik,Ncoeff))
            Dunik = np.zeros((Nzunik,Nrhounik,Ncoeff))

            # Points for evaluating integrands
            zeval = np.linspace(*self.internalknotinterval,1001)

            if df__only_z:
                if verbose:
                    print("Calculating coefficients for B_df z-component G matrix ...")
                    print(f"Total number of elements of G matrix for Bdf_z: {Nmeas*NCECS*Ncoeff}")
                    print(f"Total number of unique integrals to perform: {Nzunik*Nrhounik*Ncoeff}")

                for k in range(Nzunik):
                    if (k % 5) == 0:
                        print("{:04d}/{:04d}".format(k,Nzunik))
                    for i in range(Nrhounik):
                        for j in range(Ncoeff):

                            # Trapezoidal integration
                            Dunik[k,i,j] = Bdf_zinteg_trapz(self.bsplinefuncs[j],rhounik[i],zunik[k],zeval)

            else:
                if verbose:
                    print("Calculating coefficients for B_df matrices ...")
                    print(f"Total number of elements from combination of Bdf_rho and Bdf_z: {Nmeas*NCECS*Ncoeff*2}")
                    print(f"Total number of unique integrals to perform: {Nzunik*Nrhounik*Ncoeff*2}")

                for k in range(Nzunik):
                    if (k % 5) == 0:
                        print("{:04d}/{:04d}".format(k,Nzunik))
                    for i in range(Nrhounik):
                        for j in range(Ncoeff):
                
                            # Trapezoidal integration
                            Cunik[k,i,j] = Bdf_rhointeg_trapz(self.bsplinefuncs[j],rhounik[i],zunik[k],zeval)
                            Dunik[k,i,j] = Bdf_zinteg_trapz(self.bsplinefuncs[j],rhounik[i],zunik[k],zeval)
                
                            # Scipy's quad function (takes longer for some reason)
                            # Cunik[k,i,j] = Bdf_rhointeg_quad(self.bsplinefuncs[j],rhounik[i],zunik[k],
                            #                             self.internalknotinterval)
                            # Dunik[k,i,j] = Bdf_zinteg_quad(self.bsplinefuncs[j],rhounik[i],zunik[k],
                            #                           self.internalknotinterval)


            # Get indices into these wonder machines
            # Ultimately want Cij to have shape (N meas (Nmeas), N CECS lines, N B-spline functions)
            Cij = np.zeros((Nmeas,NCECS,Ncoeff))
            Dij = np.zeros((Nmeas,NCECS,Ncoeff))

            if verbose:
                print("Populating coefficient matrices from matrices of unique integrals  ...")

            for k in range(Nmeas):
                if (k % 1000) == 0 and debug:
                    print("{:05d}/{:05d}".format(k,Nmeas))
                    showme = True
                    showj = int(np.floor(np.random.uniform(0,Ncoeff)))  # Show a random j
                    showi = int(np.floor(np.random.uniform(0,NCECS)))  # Show a random i
                for i in range(NCECS):
                    
                    Cij[k,i,:] = Cunik[z_k[k],rho_i[k,i],:]
                    Dij[k,i,:] = Dunik[z_k[k],rho_i[k,i],:]

                    if (k % 1000) == 0 and debug:
                        # if showme and (i == showi) and (j == showj):
                        if showme and (i == showi):
                            Cquad = Bdf_rhointeg_quad(self.bsplinefuncs[showj],rho[k,i],z[k],
                                                     self.internalknotinterval)
                            Ctrapz = Bdf_rhointeg_trapz(self.bsplinefuncs[showj],rho[k,i],z[k],zeval)
                            Dquad = Bdf_zinteg_quad(self.bsplinefuncs[showj],rho[k,i],z[k],
                                                   self.internalknotinterval)
                            Dtrapz = Bdf_zinteg_trapz(self.bsplinefuncs[showj],rho[k,i],z[k],zeval)
                            
                            print(f"k, i, j: {k:4d}, {i:4d}, {showj:2d}")
                            print(f"Cunik, Cquad, Ctrapz: {Cij[k,i,showj]:14.7g}, {Cquad:14.7g}, {Ctrapz:14.7g}")
                            print(f"Dunik, Dquad, Dtrapz: {Dij[k,i,showj]:14.7g}, {Dquad:14.7g}, {Dtrapz:14.7g}")
                            print("")
                            showme = False

            # OLD, before trickery
            Ccoeff = constant / rho[:,:,np.newaxis]*Cij
            Dcoeff = constant / rho[:,:,np.newaxis]*Dij
    
            Ccoeff = constant / rho[:,:,np.newaxis]*Cij
            Dcoeff = constant / rho[:,:,np.newaxis]*Dij

            # G matrices
            Gx = Ccoeff * rhohat[:, :, np.newaxis, 0]
            Gy = Ccoeff * rhohat[:, :, np.newaxis, 1]

            Gz = Dcoeff * zhat[:, :, np.newaxis, 2] 

            # NEW trickery
            # Cij = np.transpose(Cij,axes=[0,2,1])
            # Dij = np.transpose(Dij,axes=[0,2,1])

            # Ccoeff = constant / rho[:,np.newaxis,:]*Cij
            # Dcoeff = constant / rho[:,np.newaxis,:]*Dij

            # # G matrices
            # Gx = Ccoeff * rhohat[:, np.newaxis, :, 0]
            # Gy = Ccoeff * rhohat[:, np.newaxis, :, 1]

            # Gz = Dcoeff * zhat[:, np.newaxis, :, 2] 

            Gx = Gx.reshape((Gx.shape[0],np.prod(Gx.shape[1:]))) 
            Gy = Gy.reshape((Gy.shape[0],np.prod(Gy.shape[1:]))) 
            Gz = Gz.reshape((Gz.shape[0],np.prod(Gz.shape[1:]))) 

            return Gx, Gy, Gz

        elif self.current_type == 'curl_free':
            # raise Exception("'curl_free' not implemented")

            phihat = get_phihat(x,y,self.xc,self.yc)  # shape (Nmeas, self.NCECS, 3), or (k, i, 3)
            coeff = -constant*2. / rho # Shape is (Nmeas, self.NCECS)
            

            # First, evaluate integral of B-spline functions at every z
            zunik,z_k = np.unique(z,return_inverse=True)
            nz = len(zunik)
            zmin = self.internalknotinterval[0]
            zLTzmin = zunik < zmin

            Cunik = np.array([self.bsplineintegfuncs[j](zunik) for j in range(self.Ncoeffs)]).T # Resulting shape after transpose is (nz, self.Ncoeffs)

            if np.sum(zLTzmin) > 0:
                # print("Make sure that this is right. Do you want to zero out integrals where z is less than zmin? Maybe have a look at coefficients for these 'bad' zs: Cunik[zLTzmin,:]")
                Cunik[zLTzmin,:] = 0.

            # Map back into z
            Ckj = Cunik[z_k,:]    # Shape is (Nmeas, self.Ncoeffs), or (k, j)

            # Now form G matrices
            Gx = Ckj[:,np.newaxis,:] * coeff[:,:,np.newaxis] * phihat[:,:,np.newaxis,0]
            Gy = Ckj[:,np.newaxis,:] * coeff[:,:,np.newaxis] * phihat[:,:,np.newaxis,1]
            Gz = Ckj[:,np.newaxis,:] * coeff[:,:,np.newaxis] * phihat[:,:,np.newaxis,2]

            Gx = Gx.reshape((Gx.shape[0],np.prod(Gx.shape[1:]))) 
            Gy = Gy.reshape((Gy.shape[0],np.prod(Gy.shape[1:]))) 
            Gz = Gz.reshape((Gz.shape[0],np.prod(Gz.shape[1:]))) 

            return Gx, Gy, Gz

        elif self.current_type in ['potential', 'scalar']:
            raise Exception("'potential' and 'scalar' not implemented")
            unit_vec = 1
        else:
            raise Exception('type must be "divergence_free", "curl_free", "potential", or "scalar"')
    

    def set_model_coeffs(self,m):
        if m.ndim > 1:
            assert self.m.shape == m.shape,f"Model coefficient vector must have shape (N CECS lines(={self.NCECS}), Ncoeffs(={self.Ncoeffs})) "
            m = m.ravel()
            
        self.m = m

    def get_J(self,x,y,z,
              jpar_area=None,
              constant = 1./(2*np.pi)):
        GJx, GJy = self.get_J_G_matrices(x,y,z,constant=constant)

        # Take care of laplacian, if we have any
        Nlaplacian = self.Nlaplacian_edge

        muse = self.m.reshape(self.NCECS,self.Ncoeffs)[~self.laplacian_edge,:].ravel()
        # muse = self.m[~self.laplacian_edge]
        initshape = (GJx.shape[0],self.NCECS,self.Ncoeffs)
        finalshape = (GJx.shape[0],(self.NCECS-Nlaplacian)*self.Ncoeffs)
        GJx = GJx.reshape(initshape)[:,~self.laplacian_edge,:].reshape(finalshape)
        GJy = GJy.reshape(initshape)[:,~self.laplacian_edge,:].reshape(finalshape)
        # GJx,GJy = GJx[:,~self.laplacian_edge],GJy[:,~self.laplacian_edge]
        Jx, Jy = GJx@muse,GJy@muse

        # Calculate current density in z direction
        if self.current_type == 'curl_free' and jpar_area is not None and self.grid_type == 'regular':
            
            xunik = np.unique(self.xc[~self.laplacian_edge])
            yunik = np.unique(self.yc[~self.laplacian_edge])

            # First, evaluate integral of B-spline functions at every z
            zunik = np.unique(z)
            nz = len(zunik)

            integrals_by_zunik = np.array([self.bsplineintegfuncs[j](zunik) for j in range(self.Ncoeffs)]).T # Resulting shape after transpose is (nz, self.Ncoeffs)
            
            breakpoint()
            
            # Now calculate j_z at all zunik for each CECS line
            Iz_est_atCECS = muse.reshape(self.NCECS-Nlaplacian,self.Ncoeffs)[:,np.newaxis,:]*integrals_by_zunik[np.newaxis,:,:]  # resulting shape is (NcCECS,nz,Ncoeffs)
            Jz_est_atCECS = -np.sum(Iz_est_atCECS,axis=2)/jpar_area  # Sum over coefficients, resulting shape is (NcCECS,nz)

            # … And now make an interpolating function so that we can get j_z at
            #   points that do not lie on CECS line
            # Jzinterp = rginterp((self.xc, self.yc, zunik), Jz_est_atCECS, method='linear', bounds_error=False, fill_value=np.nan)
            Jzinterp = rginterp((xunik, yunik, zunik),
                                Jz_est_atCECS.reshape(np.unique(self.xc).size,np.unique(self.yc).size,nz),
                                method='linear', bounds_error=False, fill_value=np.nan)
            
            # Finally, use interpolating function to get j_z at requested locations
            Jz = Jzinterp(np.vstack([x,y,z]).T)

        else:
            Jz = np.zeros_like(Jx)

        return np.stack([Jx,Jy,Jz]).T

    def get_B(self,x,y,z,
              constant = MU0/(4*np.pi),
              round_rho_dec_place=None):
        GBx, GBy, GBz = self.get_B_G_matrices(x,y,z,constant=constant,
                                              round_rho_dec_place=round_rho_dec_place)

        muse = self.m.copy()

        if self.Nlaplacian_edge > 0:
            muse = muse.reshape(self.NCECS,self.Ncoeffs)[~self.laplacian_edge,:].ravel()
            initshape = (GBx.shape[0],self.NCECS,self.Ncoeffs)
            finalshape = (GBx.shape[0],(self.NCECS-Nlaplacian)*self.Ncoeffs)
            GBx = GBx.reshape(initshape)[:,~self.laplacian_edge,:].reshape(finalshape)
            GBy = GBy.reshape(initshape)[:,~self.laplacian_edge,:].reshape(finalshape)
            GBz = GBz.reshape(initshape)[:,~self.laplacian_edge,:].reshape(finalshape)

        return np.stack([GBx@muse,GBy@muse,GBz@muse]).T

    def splev(h, weights, der=0, ext=1):
        """
        Evaluate spline function with particular set of weights
        """
        return splev(h, (self.knotvec,weights,self.polyorder), der=der, ext=ext)


    def splinteg_quad(hmax, weights, ext=1, epsabs=1e-5):
        """
        Integrate spline function (which represents dI/dz) up to a height hmax
        """
        tmpfunc = lambda h, t=self.knotvec, weights=weights, k=self.polyorder, ext=ext: splev(h, (t,weights,k), der=0, ext=ext)

        def tmpintegfunc(hmax, integlowbound=integlowbound, tmpfunc=tmpfunc):
            if hmax < integlowbound:
                return 0.

            return integrate.quad(tmpfunc, integlowbound, hmax, epsabs=epsabs)[0]


    def splinteg_trapz(zmax, weights, der=0, ext=1, npts=1001):

        if zmax <= self.zmin:
            return 0.

        zc = np.linspace(self.zmin,zmax,npts)

        # tmpfunc = lambda h, t=self.knotvec, weights=weights, k=self.polyorder, ext=ext: splev(h, (t,weights,k), der=0, ext=ext)
        # y = tmpfunc(zc)

        y = splev(zc, (self.knotvec,weights,self.polyorder), der=der, ext=ext)

        return integrate.trapz(y,zc)


def Bdf_rhointeg_quad(Bsplfunc_j,rho_i,z,integbounds):
    func = lambda zc, z=z, rho_i=rho_i: Bsplfunc_j(zc)*Bdfrho(rho_i,z,zc)

    # return integrate.quad(func, *integbounds,
    #                       epsabs=1e-4)[0]

    zc = np.linspace(*integbounds,1001)
    y = func(zc)

    return integrate.trapz(y,zc)


def Bdf_rhointeg_trapz(Bsplfunc_j,rho_i,z,zcs):
    func = lambda zc, z=z, rho_i=rho_i: Bsplfunc_j(zc)*Bdfrho(rho_i,z,zc)

    # return integrate.quad(func, *integbounds,
    #                       epsabs=1e-4)[0]

    y = func(zcs)

    return integrate.trapz(y,zcs)


def Bdf_zinteg_quad(Bsplfunc_j,rho_i,z,integbounds):
    return integrate.quad(lambda zc, z=z, rho_i=rho_i: Bsplfunc_j(zc)*Bdfz(rho_i,z,zc), *integbounds,
                          epsabs=1e-4)[0]


def Bdf_zinteg_trapz(Bsplfunc_j,rho_i,z,zcs):

    func = lambda zc, z=z, rho_i=rho_i: Bsplfunc_j(zc)*Bdfz(rho_i,z,zc)

    y = func(zcs)

    return integrate.trapz(y,zcs)


def Bdfrho(rho,z,zc):
    """
    """
    return (1. - np.abs(z-zc)/np.sqrt(rho**2+(z-zc)**2))*np.sign(zc-z)


def Bdfz(rho,z,zc):
    """
    """
    return -rho/np.sqrt(rho**2+(z-zc)**2)


def get_horizontal_prime_coordinates(x, y, x_cecs, y_cecs):
    
    x   = np.array(x).flatten()[:, np.newaxis]
    y   = np.array(y).flatten()[:, np.newaxis]
    x_c = np.array(x_cecs).flatten()[np.newaxis, :]
    y_c = np.array(y_cecs).flatten()[np.newaxis, :]

    xp = x-x_c
    yp = y-y_c

    return xp,yp


def get_prime_coordinates(x,y,z, x_cecs, y_cecs, z_cecs):
    
    x   = np.array(x).flatten()[:, np.newaxis]
    y   = np.array(y).flatten()[:, np.newaxis]
    z   = np.array(z).flatten()[:, np.newaxis]
    x_c = np.array(x_cecs).flatten()[np.newaxis, :]
    y_c = np.array(y_cecs).flatten()[np.newaxis, :]
    z_c = np.array(z_cecs).flatten()[np.newaxis, :]

    xp = x-x_c
    yp = y-y_c
    zp = z-z_c

    return xp,yp,zp


def get_distance(x, y, z, x_cecs, y_cecs, z_cecs):
    """" calculate distance between data point and cecs node.

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size
    return_degrees: bool, optional
        Set to True if you want output in degrees. Default is False (radians)

    Returns
    -------
    dist: 2D array (x.size, x_cecs.size)
        Array of distances between the points
        described by (x, y, z) and the points described by 
        (x_cecs, y_cecs, z_cecs).
    """

    # reshape
    xp, yp, zp = get_prime_coordinates(x,y,z,x_cecs,y_cecs,z_cecs)

    dist = np.sqrt(xp**2+yp**2+zp**2)

    return dist


def get_rho(x, y, x_cecs, y_cecs):
    """" calculate rho (horizontal distance) between data points and cecs nodes.

    Parameters
    ----------
    x,y: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size

    Returns
    -------
    rho: 2D array (x.size, x_cecs.size)
        Array of horizontal distances between data points
        described by (x, y, z) and cecs nodes at (x_cecs, y_cecs, z_cecs).
    """

    xp, yp = get_horizontal_prime_coordinates(x, y, x_cecs, y_cecs)

    rho = np.sqrt( xp**2+yp**2)

    return rho


def get_phi(x, y, x_cecs, y_cecs, return_degrees = False):
    """" calculate polar (horizontal) angle of data point relative to cecs node.

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size
    return_degrees: bool, optional
        Set to True if you want output in degrees. Default is False (radians)

    Returns
    -------
    phi: 2D array (x.size, x_cecs.size)
        Array of polar (horizontal) angles of data points
        described by (x, y, z) relative the points described by 
        (x_cecs, y_cecs, z_cecs). Unit in radians unless return_degrees is set
        to True
    """

    xp, yp = get_horizontal_prime_coordinates(x,y,x_cecs,y_cecs)

    phi = np.arctan2(yp,xp)

    if return_degrees:
        phi = phi/d2r

    return phi


def get_phihat(x,y, x_cecs, y_cecs):
    """" calculate polar (horizontal) unit vector for data point relative to cecs node.

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size

    Returns
    -------
    phihat: 3D array (x.size, x_cecs.size, 3)
        Array of polar (horizontal) unit vectors in Cartesian coordinates of data points
        described by (x, y, z) relative the points described by 
        (x_cecs, y_cecs, z_cecs).
    """
    
    phiprime = get_phi(x, y, x_cecs, y_cecs)
    
    phihat = np.transpose(np.stack([-np.sin(phiprime),np.cos(phiprime),np.zeros_like(phiprime)]),axes=[1,2,0])

    return phihat


def get_rhohat(x, y, x_cecs, y_cecs):
    """" calculate rho unit vector from cecs node to data point.

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size
    lon_cecs: array-like
        Array of CECS pole longitudes [deg]
        Flattened array must havef same size as lat_cecs
        Output will be a 2D array with shape (mlat.size, mlat_cecs.size)

    Returns
    -------
    rhohat: 3D array (3, x.size, x_cecs.size)
        Array of rho unit vectors in Cartesian coordinates pointing from
        cecs nodes to data points.
    """
    
    phiprime = get_phi(x, y, x_cecs, y_cecs)
    
    rhohat = np.transpose(np.stack([np.cos(phiprime),np.sin(phiprime),np.zeros_like(phiprime)]),axes=[1,2,0])

    return rhohat


def get_CECS_J_G_matrices(x,y,z, x_cecs, y_cecs, z_cecs, 
                          current_type = 'divergence_free', constant = 1./(2*np.pi), 
                          dz_tolerance=10.,
):
    """ Calculate matrices Gx and Gy which relate CECS amplitudes to current density 
        vector components.

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size
    current_type: string, optional
        The type of CECS function. This must be either 
        'divergence_free' (default): divergence-free basis functions
        'curl_free': curl-free basis functions
        'potential': scalar field whose negative gradient is curl-free CECS
        'scalar': 
    constant: float, optional
        The CECS functions are scaled by the factor 1/(2pi), which is
        the default value of 'constant'. Change if you want something 
        different.
    dz_tolerance: Since CECS were not originally designed to be used in 3D, we have to make sure that
                      in evaluating the current at a particular altitude, the CECS poles at higher or lower
                      altitudes (within tolerance given by dz_tolerance) do not contribute to the
                      evaluated current

    Returns
    -------
    If current_type is 'divergence_free' or 'curl_free':
    Gx: 2D array
        2D array with shape (x.size, x_cecs.size), relating CECS amplitudes
        m to the x-direction current densities at (x, y, z) via 'jx = Gx.dot(m)'
    Gy: 2D array
        2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
        m to the y-direction current densities at (x, y, z) via 'jy = Gy.dot(m)'
    # Gz: 2D array
    #     2D array with shape (x.size, y_cecs.size), relating CECS amplitudes
    #     m to the z-direction current densities at (x, y, z) via 'jz = Gz.dot(m)'
    If current_type is 'potential' or 'scalar':
    G: 2D array
        2D array with shape (lat.size, lat_cecs.size), relating amplitudes m
        to scalar field magnitude at (lat, lon) via 'z = G.dot(m)'

    """
    
    if current_type not in ['divergence_free','curl_free']:
        assert 2<0,"types 'potential' and 'scalar' not yet implemented!"

    xp, yp, zp = get_prime_coordinates(x,y,z,x_cecs,y_cecs,z_cecs)

    rho = get_rho(x,y,x_cecs,y_cecs)

    if current_type == 'divergence_free':
        unit_vec = (-1.)*get_phihat(x,y,x_cecs,y_cecs)
    elif current_type == 'curl_free':
        unit_vec = get_rhohat(x,y,x_cecs,y_cecs)
    elif current_type in ['potential', 'scalar']:
        unit_vec = 1
    else:
        raise Exception('type must be "divergence_free", "curl_free", "potential", or "scalar"')

    # get the scalar part of Amm's divergence-free CECS:    
    if current_type in ['divergence_free', 'curl_free']:
        coeff = constant / rho

        # G matrices
        Gx = coeff * unit_vec[:, :, 0]
        Gy = coeff * unit_vec[:, :, 1]
    
        shouldbezero = np.abs(zp) > dz_tolerance
        Gx[shouldbezero] = 0.
        Gy[shouldbezero] = 0.

        return Gx, Gy
    else: # current_type is 'potential' or 'scalar'
        # NOT SURE WHAT THESE SHOULD BE
        # if current_type == 'potential':
        #     return -2*constant*np.log(np.sin(theta/2)).T
        # elif current_type == 'scalar':
        #     return    constant      / np.tan(theta/2).T

        return np.nan, np.nan


def get_CECS_B_G_matrices(x, y, z, x_cecs, y_cecs, z_cecs,
                          current_type = 'divergence_free',
                          constant = MU0/(4.*np.pi), 
):
    """ Calculate matrices Gx, Gy, and Gz that relate CECS amplitudes to magnetic field

    Based on equations (A.10) and (A.11) of Vanhamäki (2007).

    Parameters
    ----------
    x,y,z: array-like
        Array of Cartesian coordinate of evaluation points [m]
        Flattened arrays must all have the same size
    {x,y,z}_cecs: array-like
        Array of CECS pole Cartesian coordinates [m]
        Flattened arrays must all have the same size
    current_type: string, optional
        The type of CECS function. This must be either 
        'divergence_free' (default): divergence-free basis functions
        'curl_free': curl-free basis functions
    constant: float, optional
        The B^cf CECS function is scaled by the factor MU0/(2pi), while
        B^df CECS function is scaled by MU0/(4pi)

    Returns
    -------
    Gx: 2D array
        2D array with shape (x.size, x_cecs.size), relating CECS amplitudes
        m to the x-component magnetic field at (x, y, z) via 'Bx = Gx.dot(m)'
    Gy: 2D array
        2D array with shape (y.size, y_cecs.size), relating CECS amplitudes
        m to the y-component magnetic field at (x, y, z) via 'By = Gy.dot(m)'
    Gz: 2D array
        2D array with shape (z.size, z_cecs.size), relating CECS amplitudes
        m to the z-component magnetic field at (x, y, z) via 'Bz = Gz.dot(m)'


    Notes
    ------
    Variables with 'p' in name signify a quantity in the CECS node's local (i.e., "p"rimed) coordinate system

    2022/10/15 SMH
    """
    
    # reshape z:
    if np.array(z).size == 1:
        z = np.ones_like(x) * z
    else:
        z = np.array(z).flatten()[:, np.newaxis]

    xp, yp, zp = get_prime_coordinates(x,y,z,x_cecs,y_cecs,z_cecs)

    rho = get_rho(x,y,z,x_cecs,y_cecs,z_cecs)

    coeff = constant / rho

    # G matrix scale factors
    if current_type == 'divergence_free':

        rhohat = get_rhohat(x,y,z,x_cecs,y_cecs,z_cecs)
        zhat   = np.zeros_like(rhohat)
        zhat[:,:,2] = 1.

        Grhoprime = -coeff * (1 - np.abs(zp) / np.sqrt(rho**2+zp**2)) * np.sign(zp)

        Gx = Grhoprime * rhohat[:,:,0]
        Gy = Grhoprime * rhohat[:,:,1]

        Gz = -coeff * rho/np.sqrt(rho**2+zp**2)


    elif current_type == 'curl_free':

        phihat = get_phihat(x,y,z,x_cecs,y_cecs,z_cecs)

        heaviside = np.zeros_like(zp,dtype=zp.dtype)
        heaviside[((zp) > -1e8) | (np.isclose(0.,zp))] = 1.

        Gphiprime = -2. * coeff * heaviside

        Gx = Gphiprime * phihat[:,:,0]
        Gy = Gphiprime * phihat[:,:,1]
        Gz = np.zeros_like(phihat[:,:,2])

    return Gx, Gy, Gz


if __name__ == '__main__':
    import hatch_python_utils.math.cecs_utils
    importlib.reload(hatch_python_utils.math.cecs_utils)
    from hatch_python_utils.math import cecs_utils
    import importlib 
    importlib.reload(cecs_utils)
    print("TEST OUT CECS")
    # x,y,z = np.array([[1,0,0],[2,0,0],[3,0,0]]).T
    # xc,yc,zc = np.array([[0,0,-1],[0,0,0],[0,0,1]]).T
    xyz = np.array([[1,0,0],[1,1,0],[-1,1,1]]).T
    xyzc = np.array([[0,0,-1],[0,0,0],[0,0,1],[1,0,1],[-5,0,1]]).T
    # cecs_utils.get_phi(x,y,z,xc,yc,zc)
    phi = cecs_utils.get_phi(*xyz,*xyzc,return_degrees=True)
    phihat = cecs_utils.get_phihat(*xyz,*xyzc)
    rhohat = cecs_utils.get_rhohat(*xyz,*xyzc)
    Gx, Gy = cecs_utils.get_CECS_J_G_matrices(*xyz,*xyzc)

    GBx, GBy, GBz = cecs_utils.get_CECS_B_G_matrices(*xyz,*xyzc,current_type='divergence_free')

    print(phi.shape)
    print(phi)


