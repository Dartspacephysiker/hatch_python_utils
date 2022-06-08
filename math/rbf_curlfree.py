import numpy as np
from hatch_python_utils.math.gauss3d import gaussian3d,gauss3d_gradient_of_laplacian_ith

MU0 = 1.25663706212e-06

class RBFs(object):

    def __init__(self,nodevectors,
                 cvectors,
                 sigmavectors):
        """
        nodevectors  : positions of curl-free RBF nodes in Cartesian coords

        cvectors     : weighting coefficients for each RBF node

        sigmavectors : Vector of sigma values

        For nodevectors, cvectors, and sigmavectors:
                        shape is [M,3], where
                        M is number of RBF nodes
                        '3' says these vectors are in 3D Cartesian coordinates

        SMH
        2022/06/07
        """

        self.nodevectors = nodevectors
        self.cvectors = cvectors
        self.nuvectors = 1/(np.sqrt(2)*sigmavectors)

        self.NRBFs = self.nodevectors.shape[0]

    def get_divergence(self,posvectors):
        """
        divVvec = rbfs.get_divE(posvectors)
        Units of divEvec are [unit of cvectors]/(FIND OUT)

        For example, if c vectors (which are ??) are in V-km and position is in km,
        divVvec has units V/km^2
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        GdivE = self.get_divE_Gmatrix(posvectors)

        return (GdivE@(self.cvectors.T.ravel())).reshape(posvectors.T.shape).T 

    def get_field(self,posvectors):
        """
        Evec = rbfs.get_field(posvectors)
        Units of E are given by [unit of cvectors]/[(position units)^2]

        For example, if c vectors are in V-km, Vvec is in V/km
        """

        # relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]
        # return V_gauss3d_vec(relposvec,self.cvectors,self.nuvectors)

        # MATRIX WAY
        G = self.get_field_Gmatrix(posvectors)
        return (G@(self.cvectors.T.ravel())).reshape(posvectors.T.shape).T

    def get_numerical_divergence(self,posvectors,delta=0.01):
        """
        div_E = rbfs.get_numerical_divergence(posvectors)
        Units of E are given by [unit of cvectors]/[unit of posvectors and sigma vectors]**2

        For example, if c vectors are in V-km and posvectors are in km, Evec is in V/km
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        blankp = relposvec*0
        ip,jp,kp = blankp.copy(),blankp.copy(),blankp.copy()
        ip[:,:,0] = delta/2
        jp[:,:,1] = delta/2
        kp[:,:,2] = delta/2

        dEidx = (V_gauss3d_ith(relposvec+ip,self.cvectors,0,self.nuvectors).sum(axis=1)-V_gauss3d_ith(relposvec-ip,self.cvectors,0,self.nuvectors).sum(axis=1))/delta
        dEjdy = (V_gauss3d_ith(relposvec+jp,self.cvectors,1,self.nuvectors).sum(axis=1)-V_gauss3d_ith(relposvec-jp,self.cvectors,1,self.nuvectors).sum(axis=1))/delta
        dEkdz = (V_gauss3d_ith(relposvec+kp,self.cvectors,2,self.nuvectors).sum(axis=1)-V_gauss3d_ith(relposvec-kp,self.cvectors,2,self.nuvectors).sum(axis=1))/delta

        return dEidx+dEjdy+dEkdz


    def get_numerical_curl(self,posvectors,delta=0.01):
        """
        div_E = rbfs.get_numerical_curl(posvectors)
        Units of E are given by [unit of cvectors]/[unit of posvectors and sigma vectors]**2

        For example, if c vectors are in V-km and posvectors are in km, Evec is in V/km
        """

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        blankp = relposvec*0
        ip,jp,kp = blankp.copy(),blankp.copy(),blankp.copy()
        ip[:,:,0] = delta/2
        jp[:,:,1] = delta/2
        kp[:,:,2] = delta/2

        i,j,k = 0,1,2

        # i term
        dEjdz = V_gauss3d_ith(relposvec+kp,self.cvectors,j,self.nuvectors) - V_gauss3d_ith(relposvec-kp,self.cvectors,j,self.nuvectors)
        dEkdy = V_gauss3d_ith(relposvec+jp,self.cvectors,k,self.nuvectors) - V_gauss3d_ith(relposvec-jp,self.cvectors,k,self.nuvectors)

        # j term
        dEkdx = V_gauss3d_ith(relposvec+ip,self.cvectors,k,self.nuvectors) - V_gauss3d_ith(relposvec-ip,self.cvectors,k,self.nuvectors)
        dEidz = V_gauss3d_ith(relposvec+kp,self.cvectors,i,self.nuvectors) - V_gauss3d_ith(relposvec-kp,self.cvectors,i,self.nuvectors)

        # k term
        dEidy = V_gauss3d_ith(relposvec+jp,self.cvectors,i,self.nuvectors) - V_gauss3d_ith(relposvec-jp,self.cvectors,i,self.nuvectors)
        dEjdx = V_gauss3d_ith(relposvec+ip,self.cvectors,j,self.nuvectors) - V_gauss3d_ith(relposvec-ip,self.cvectors,j,self.nuvectors)

        iterm = dEjdz-dEkdy
        jterm = dEkdx-dEidz
        kterm = dEidy-dEjdx

        return np.hstack([iterm,jterm,kterm])


    def get_field_Gmatrix(self,posvectors):

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        return V_gauss3d_Gmatrix(relposvec,self.nuvectors)


    def get_fieldscalar_Gmatrix(self,posvectors,hatvectors):

        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        nuvec = self.nuvectors

        Gi = V_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
        Gj = V_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
        Gk = V_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

        G = Gi*hatvectors[:,0][:,np.newaxis] + \
            Gj*hatvectors[:,1][:,np.newaxis] + \
            Gk*hatvectors[:,2][:,np.newaxis]

        return G

    def get_divergence_Gmatrix(self,posvectors):
        
        relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

        return divV_gauss3d_Gmatrix(relposvec,self.nuvectors)


    # def get_divEscalar_Gmatrix(self,posvectors,jhatvectors):
        
    #     relposvec = posvectors[:,np.newaxis,:]-self.nodevectors[np.newaxis,:,:]

    #     nuvec = self.nuvectors

    #     Gi = divV_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
    #     Gj = divV_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
    #     Gk = divV_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

    #     G = Gi*jhatvectors[:,0][:,np.newaxis] + \
    #         Gj*jhatvectors[:,1][:,np.newaxis] + \
    #         Gk*jhatvectors[:,2][:,np.newaxis]

    #     return G


# def divV_gauss3d_ith(relposvec,cvec,i,nuvec):

#     j,k = (i+1)%3,(i+2)%3

#     jderiv = gauss3d_gradient_of_laplacian_ith(relposvec,j,nuvec)
#     kderiv = gauss3d_gradient_of_laplacian_ith(relposvec,k,nuvec)

#     current_ith = cvec[:,j]*kderiv-cvec[:,k]*jderiv

#     return current_ith


def divV_gauss3d_Gmatrix_ith(relposvec,i,nuvec):

    j,k = (i+1)%3,(i+2)%3

    ideriv = gauss3d_gradient_of_laplacian_ith(relposvec,i,nuvec)
    jderiv = gauss3d_gradient_of_laplacian_ith(relposvec,j,nuvec)
    kderiv = gauss3d_gradient_of_laplacian_ith(relposvec,k,nuvec)

    G_i = [-ideriv,
           -jderiv,
           -kderiv]

    if i == 0:
        G_i = np.hstack(G_i)
    elif i == 1:
        G_i = np.hstack([G_i[2],G_i[0],G_i[1]])
    elif i == 2:
        G_i = np.hstack([G_i[1],G_i[2],G_i[0]])

    return G_i


def divV_gauss3d_Gmatrix(relposvec,nuvec):
    
    Gi = divV_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
    Gj = divV_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
    Gk = divV_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

    return np.vstack([Gi,Gj,Gk])


def V_gauss3d_ith(relposvec,cvec,i,nuvec):
    """
    Calculate the ith component (in Cartesian coords) of the curl-free field at N locations due to
    M RBF nodes.
    """

    gaussvals = gaussian3d(relposvec,nuvec)

    j,k = (i+1)%3,(i+2)%3

    ci,cj,ck = cvec[:,i],cvec[:,j],cvec[:,k]
    xi,xj,xk = relposvec[:,:,i],relposvec[:,:,j],relposvec[:,:,k]
    nui,nuj,nuk = nuvec[:,i],nuvec[:,j],nuvec[:,k]

    # Faster?
    E_i = ( ci * ( -4*(nui*nui)**2*xi*xi + 2*nui**2) + \
            cj * ( -4*(nui*nuj)**2*xi*xj           ) + \
            ck * ( -4*(nui*nuk)**2*xi*xk           ) )*gaussvals

    return E_i


def V_gauss3d_vec(relposvec,cvec,nuvec):
    
    Ei = np.sum(V_gauss3d_ith(relposvec,cvec,0,nuvec),axis=1)
    Ej = np.sum(V_gauss3d_ith(relposvec,cvec,1,nuvec),axis=1)
    Ek = np.sum(V_gauss3d_ith(relposvec,cvec,2,nuvec),axis=1)

    return np.vstack([Ei,Ej,Ek]).T


def V_gauss3d_Gmatrix_ith(relposvec,i,nuvec):
    
    gaussvals = gaussian3d(relposvec,nuvec)

    j,k = (i+1)%3,(i+2)%3

    xi,xj,xk = relposvec[:,:,i],relposvec[:,:,j],relposvec[:,:,k]
    nui,nuj,nuk = nuvec[:,i],nuvec[:,j],nuvec[:,k]

    G_i = [( -4*(nui*nui)**2*xi*xi + 2*nui**2)*gaussvals,
           ( -4*(nui*nuj)**2*xi*xj           )*gaussvals,
           ( -4*(nui*nuk)**2*xi*xk           )*gaussvals]

    # Gotta make sure things end up in the right rows ...
    if i == 0:
        G_i = np.hstack(G_i)
    elif i == 1:
        G_i = np.hstack([G_i[2],G_i[0],G_i[1]])
    elif i == 2:
        G_i = np.hstack([G_i[1],G_i[2],G_i[0]])

    return G_i


def V_gauss3d_Gmatrix(relposvec,nuvec):
    
    Gi = V_gauss3d_Gmatrix_ith(relposvec,0,nuvec)
    Gj = V_gauss3d_Gmatrix_ith(relposvec,1,nuvec)
    Gk = V_gauss3d_Gmatrix_ith(relposvec,2,nuvec)

    return np.vstack([Gi,Gj,Gk])

