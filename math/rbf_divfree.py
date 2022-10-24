import numpy as np
from hatch_python_utils.math.gauss3d import gaussian3d,gauss3d_gradient_of_laplacian_ith

MU0 = 1.25663706212e-06

class RBFs(object):

    def __init__(self,nodevectors,
                 sigmavectors,
                 cvectors=None):
        """
        nodevectors  : positions of divergence-free RBF nodes in Cartesian coords

        cvectors     : weighting coefficients for each RBF node

        sigmavectors : Vector of sigma values

        For nodevectors, cvectors, and sigmavectors:
                        shape is [M,3], where
                        M is number of RBF nodes
                        '3' says these vectors are in 3D Cartesian coordinates

        SMH
        2022/06/02
        """

        self.nodevectors = nodevectors
        self.nuvectors = 1/(np.sqrt(2)*sigmavectors)

        if cvectors is not None:
            self.set_cvectors(cvectors)
        else:
            self.set_cvectors(np.zeros(self.nodevectors.shape))

        self.NRBFs = self.nodevectors.shape[0]

    def set_cvectors(self,cvectors):
        self.cvectors = cvectors

    # def get_curl(self,posvectors):
    #     """
    #     Jvec = rbfs.get_curl(posvectors)
    #     Units of Jvec are [unit of cvectors]/(mu0*[unit of position vectors])

    #     For example, if c vectors (which are magnetic field) are in nT and position is in km,
    #     Jvec has units nT/mu0/km
    #     """

    #     relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

    #     return J_gauss3d_vec(relposvec,self.cvectors,self.nuvectors)


    def get_curl(self,posvectors):
        """
        Jvec = rbfs.get_curl(posvectors)
        Units of Jvec are [unit of cvectors]/(mu0*[unit of position vectors])

        For example, if c vectors (which are magnetic field) are in nT and position is in km,
        Jvec has units nT/mu0/km
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        GJ = self.get_curl_Gmatrix(posvectors)

        return (GJ@(self.cvectors.T.ravel())).reshape(posvectors.T.shape).T 


    def get_field(self,posvectors):
        """
        Bvec = rbfs.get_field(posvectors)
        Units of B are given by [unit of cvectors]

        For example, if c vectors are in nT, so is Bvec
        """

        # OLD
        # relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]
        # return B_gauss3d_vec(relposvec,self.cvectors,self.nuvectors)

        # MATRIX WAY
        G = self.get_field_Gmatrix(posvectors)
        return (G@(self.cvectors.T.ravel())).reshape(posvectors.T.shape).T


    def get_divergence(self,posvectors,delta=0.01):
        """
        div_B = rbfs.get_divergence(posvectors)
        Units of B are given by [unit of cvectors]/[unit of posvectors and sigma vectors]**2

        For example, if c vectors are in nT-km^2 and posvectors are in km, Bvec is in nT
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        blankp = relposvec*0
        ip,jp,kp = blankp.copy(),blankp.copy(),blankp.copy()
        ip[:,:,0] = delta/2
        jp[:,:,1] = delta/2
        kp[:,:,2] = delta/2

        # Bi = (B_gauss3d_ith(relposvec+ip,self.cvectors,0,self.nuvectors)-B_gauss3d_ith(relposvec-ip,self.cvectors,0,self.nuvectors)).sum(axis=1)/delta
        # Bj = (B_gauss3d_ith(relposvec+jp,self.cvectors,1,self.nuvectors)-B_gauss3d_ith(relposvec-jp,self.cvectors,1,self.nuvectors)).sum(axis=1)/delta
        # Bk = (B_gauss3d_ith(relposvec+kp,self.cvectors,2,self.nuvectors)-B_gauss3d_ith(relposvec-kp,self.cvectors,2,self.nuvectors)).sum(axis=1)/delta

        dBidx = (B_gauss3d_ith(relposvec+ip,self.cvectors,0,self.nuvectors).sum(axis=1)-B_gauss3d_ith(relposvec-ip,self.cvectors,0,self.nuvectors).sum(axis=1))/delta
        dBjdy = (B_gauss3d_ith(relposvec+jp,self.cvectors,1,self.nuvectors).sum(axis=1)-B_gauss3d_ith(relposvec-jp,self.cvectors,1,self.nuvectors).sum(axis=1))/delta
        dBkdz = (B_gauss3d_ith(relposvec+kp,self.cvectors,2,self.nuvectors).sum(axis=1)-B_gauss3d_ith(relposvec-kp,self.cvectors,2,self.nuvectors).sum(axis=1))/delta

        # print("delta: ",delta)
        # print("Bi:")
        # print(Bi)
        # print("Bj:")
        # print(Bj)
        # print("Bk:")
        # print(Bk)
        # print("Bi+Bj+Bk:")
        # print(Bi+Bj+Bk)

        return dBidx+dBjdy+dBkdz


    def get_field_Gmatrix(self,posvectors):

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        return B_gauss3d_Gmatrix(relposvec,self.nuvectors)


    def get_fieldscalar_Gmatrix(self,posvectors,hatvectors):
        """
        Same as output for get_field_Gmatrix, except that now we assume each best-estimate vector is dotted with a hatvector
        Useful, for example, when you only measure one component of a divergence-free field
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        nuvec = self.nuvectors

        Gi = B_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
        Gj = B_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
        Gk = B_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

        G = Gi*hatvectors[:,0][:,np.newaxis] + \
            Gj*hatvectors[:,1][:,np.newaxis] + \
            Gk*hatvectors[:,2][:,np.newaxis]

        return G


    def get_curl_Gmatrix(self,posvectors):
        
        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        return J_gauss3d_Gmatrix(relposvec,self.nuvectors)


    def get_curlscalar_Gmatrix(self,posvectors,jhatvectors):
        
        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        nuvec = self.nuvectors

        Gi = J_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
        Gj = J_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
        Gk = J_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

        G = Gi*jhatvectors[:,0][:,np.newaxis] + \
            Gj*jhatvectors[:,1][:,np.newaxis] + \
            Gk*jhatvectors[:,2][:,np.newaxis]

        return G


# def J_gauss3d_ith(relposvec,cvec,i,nuvec):

#     j,k = (i+1)%3,(i+2)%3

#     jderiv = gauss3d_gradient_of_laplacian_ith(relposvec,j,nuvec)
#     kderiv = gauss3d_gradient_of_laplacian_ith(relposvec,k,nuvec)

#     current_ith = cvec[:,j]*kderiv-cvec[:,k]*jderiv

#     return current_ith


def J_gauss3d_Gmatrix_ith(relposvec,i,nuvec):

    j,k = (i+1)%3,(i+2)%3

    # ideriv = jderiv*0
    jderiv = gauss3d_gradient_of_laplacian_ith(relposvec,j,nuvec)
    kderiv = gauss3d_gradient_of_laplacian_ith(relposvec,k,nuvec)

    G_i = [np.zeros(kderiv.shape),
           kderiv,
           -jderiv]

    if i == 0:
        G_i = np.hstack(G_i)
    elif i == 1:
        G_i = np.hstack([G_i[2],G_i[0],G_i[1]])
    elif i == 2:
        G_i = np.hstack([G_i[1],G_i[2],G_i[0]])

    return G_i


def J_gauss3d_Gmatrix(relposvec,nuvec):
    
    Gi = J_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
    Gj = J_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
    Gk = J_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

    return np.vstack([Gi,Gj,Gk])


# def J_gauss3d_vec(relposvec,cvec,nuvec):
#     """
#     See RBFs.get_curl for comment on units
#     """

#     Ji = J_gauss3d_ith(relposvec,cvec,0,nuvec).sum(axis=1)
#     Jj = J_gauss3d_ith(relposvec,cvec,1,nuvec).sum(axis=1)
#     Jk = J_gauss3d_ith(relposvec,cvec,2,nuvec).sum(axis=1)

#     return np.vstack([Ji,Jj,Jk]).T


def B_gauss3d_ith(relposvec,cvec,i,nuvec):
    """
    Calculate the ith component (in Cartesian coords) of the div-free field at N locations due to
    M RBF nodes.
    """

    gaussvals = gaussian3d(relposvec,nuvec)

    j,k = (i+1)%3,(i+2)%3

    ci,cj,ck = cvec[:,i],cvec[:,j],cvec[:,k]
    xi,xj,xk = relposvec[:,:,i],relposvec[:,:,j],relposvec[:,:,k]
    nui,nuj,nuk = nuvec[:,i],nuvec[:,j],nuvec[:,k]

    # Faster?
    B_i = ( -ci * ( 2*nuj**2*(2*(nuj*xj)**2-1) + 2*nuk**2*(2*(nuk*xk)**2-1) ) \
            + cj * (4 * nui**2 * nuj**2 * xi * xj) \
            + ck * (4 * nui**2 * nuk**2 * xi * xk))*gaussvals

    return B_i


# def B_gauss3d_vec(relposvec,cvec,nuvec):
#     """
#     See RBFs.get_Bfield for comment on units
#     """
    
#     Bi = B_gauss3d_ith(relposvec,cvec,0,nuvec).sum(axis=1)
#     Bj = B_gauss3d_ith(relposvec,cvec,1,nuvec).sum(axis=1)
#     Bk = B_gauss3d_ith(relposvec,cvec,2,nuvec).sum(axis=1)

#     return np.vstack([Bi,Bj,Bk]).T


def B_gauss3d_Gmatrix_ith(relposvec,i,nuvec):
    
    gaussvals = gaussian3d(relposvec,nuvec)

    j,k = (i+1)%3,(i+2)%3

    xi,xj,xk = relposvec[:,:,i],relposvec[:,:,j],relposvec[:,:,k]
    nui,nuj,nuk = nuvec[:,i],nuvec[:,j],nuvec[:,k]

    G_i = [-( 2*nuj**2*(2*(nuj*xj)**2-1) + 2*nuk**2*(2*(nuk*xk)**2-1) )*gaussvals,
           (4 * nui**2 * nuj**2 * xi * xj)*gaussvals,
           (4 * nui**2 * nuk**2 * xi * xk)*gaussvals                             ]

    # Gotta make sure things end up in the right rows ...
    if i == 0:
        G_i = np.hstack(G_i)
    elif i == 1:
        G_i = np.hstack([G_i[2],G_i[0],G_i[1]])
    elif i == 2:
        G_i = np.hstack([G_i[1],G_i[2],G_i[0]])

    return G_i


def B_gauss3d_Gmatrix(relposvec,nuvec):
    
    Gi = B_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
    Gj = B_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
    Gk = B_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

    return np.vstack([Gi,Gj,Gk])

