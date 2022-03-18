def normal_vector_polar(phi,rho,drhodphi):
    """
    This function returns the radial (or colatitudinal) and azimuthal components
    of the normal vector for an implicitly defined surface f(r,φ) = r-ρ(φ) == 0 in polar coordinates.

    The normal vector is given by the gradient of f: 
    ∇f/|∇f| = <∂f/∂r, 1/r ∂f/∂φ  > / |∇f| 
            = <  1  , -1/ρ ∂ρ/∂φ > /  sqrt(1+(∂ρ/∂φ)²/ρ²)

    phi       : Polar angle(s) in deg where we evaluate normal vector
    rho       : Vector-compatible function of phi of form   f(phi)
    drhodphi  : Vector-compatible function of phi of form   g(phi)


    NOTE: 
    The rho component can also represent the colatitude (or southward) component and 
    The phi component can represent the eastward component.

    RETURNS
    ======
    (rho component, phi component)

    """
    
    import numpy as np
    
    #Adjust phis that are near 0°, 90° 180°, and 270° to avoid singularity
    # checkem = np.array([0.,90.,180.,270.,360.])
    # for check in checkem:
    #     phi = np.where(np.isclose(phi%90,0),phi+0.0001,phi)

    # Get radii
    r = rho(phi)

    #Get deriv. of f wrt φ
    der = -1/r * drhodphi(phi)

    # Get norm of gradient
    vnorm = 1/np.sqrt(1+der**2)

    # r/colatitude/southward component
    rc = vnorm

    # phi/eastward component
    pc = der*vnorm

    return rc,pc


def full_normal_vector_polar(r,phi,f,dfdr,dfdphi):
    """
    This function returns the radial (or colatitudinal) and azimuthal components
    of the normal vector for an implicitly defined surface f(r,φ) == 0 in polar coordinates.

    The normal vector is given by the gradient of f: ∇f = (∂f/∂r, 1/r ∂f/∂φ).

    phi   : Polar angle(s) in deg where we evaluate normal vector
    f     : Vector-compatible function of r and phi of form f(r,phi)
    dfdr  : Vector-compatible function of r and phi of form g(r,phi)
    dfdphi: Vector-compatible function of r and phi of form h(r,phi)
    """
    
    import numpy as np
    
    #Adjust phis that are near 0°, 90° 180°, and 270° to avoid singularity
    # checkem = np.array([0.,90.,180.,270.,360.])
    # for check in checkem:
    #     phi = np.where(np.isclose(phi%90,0),phi+0.0001,phi)

    # r/colatitude component
    rc = dfdr(r,phi)

    # r/colatitude component
    phic = dfdphi(r,phi)    

