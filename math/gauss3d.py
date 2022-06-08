import numpy as np

def gaussian3d(relposvec,nuvec):
    """
    relposvec: shape is [N,M,3], where 
                    N is number of calculation points, 
                    M is number of RBF nodes, and 
                    '3' says these vectors are in 3D Cartesian coordinates
    nuvec : shape is [M,3] — M vectors of np.sqrt(2)/sigma values

    """

    # frontfac = np.prod(nuvec/np.sqrt(np.pi),axis=1)
    frontfac = 1
    return frontfac*np.exp(-np.sum((relposvec*nuvec)**2,axis=2))


def laplacian_gaussian3d(relposvec,nuvec):

    gaussvals = gaussian3d(relposvec,nuvec)
    
    factor = np.zeros(relposvec.shape[:2])
    for i in range(3):
        nu = nuvec[:,i]
        posi = relposvec[:,:,i]
        factor += 4*nu**4*posi**2-2*nu**2
    
    return factor*gaussvals
    

def gauss3d_gradient_of_laplacian_ith(relposvec,i,nuvec):
    """
    ∂_i (∂_m ∂_m ϕ), that's what we calculate here
    """

    gaussvals = gaussian3d(relposvec,nuvec)
    laplace = laplacian_gaussian3d(relposvec,nuvec)
        
    nui,xi = nuvec[:,i],relposvec[:,:,i]

    # First Laplacian bit
    iterm = laplace*(-2*nui**2*xi)
    
    # Then the extra
    iterm += 8*nui**4*xi*gaussvals
    
    return iterm


