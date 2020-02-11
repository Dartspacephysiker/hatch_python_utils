import numpy as np

def closest_point_in_plane(n, d, p):
    """
    Parameters
    ==========
    n              : array, the vector [A,B,C] that defines the plane
    d              : float, the distance of the plane from the origin
    p              : array, 3XN vector of points in Cartesian <X,Y,Z> coordinates for which you want nearest plane point 

    print("EXAMPLE 1")
    print("=========")
    n = np.array([0,0,1])           # +z unit vector
    d = 10
    p = np.array([10,10,0])
    x = closest_point_in_plane(n,d,p)
    print("Calcked answer:",x)
    print("Right answer:",[p[0],p[1],d])
    print("")
    """

    assert n.shape[0] == 3,"Normal vector must have three components!"

    isSingleVector = len(p.shape) == 1
    if isSingleVector:
        assert p.shape[0] == 3
        p = p[np.newaxis].T

    elif len(p.shape) == 2:
        
        transposeme = p.shape[0] != 3

        if transposeme:
            assert p.shape[1] == 3,"Provide 3XN or Nx3 vector of points!"

            p = p.T

        else:
            if p.shape[0] == p.shape[1]:
                print("closest_point_in_plane: FARLIG! Provided 3x3 array!")

    v = (d - np.sum(p*n[:,np.newaxis],axis=0)) / np.sum(n*n)
    x = p + v[np.newaxis,:]*n[:,np.newaxis]

    if isSingleVector:
        return x.flatten()
    else:
        if transposeme:
            return x.T
        else:
            return x	

def get_plane_fit(xs,ys,zs,
                  reference_ax=None,
                  return_XYZ_meshgrid=True):
    """
    # Ripped from https://stackoverflow.com/questions/12299540/plane-fitting-to-4-or-more-xyz-points
    """

    # do fit
    tmp_A = []
    tmp_b = []
    for i in range(len(xs)):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_b.append(zs[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit = (A.T * A).I * A.T * b # DON'T MESS WITH FIT; you use it for calculating Z values below
    errors = b - A * fit
    residual = np.linalg.norm(errors)
    
    print("solution:")
    print("%f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
    # print("errors:",errors)
    # print("residual:",residual)
    
    normvec = np.array([fit[0].item(),fit[1].item(),1])
    normfactor = np.linalg.norm(normvec)
    normvec = normvec/normfactor

    outdict = dict(fit=fit,
                   errors=errors,
                   residual=residual,
                   x=fit[0].item(),
                   y=fit[1].item(),
                   z=1,
                   d=fit[2].item(),
                   normfactor=normfactor,
                   normvec=normvec)

    if (not return_XYZ_meshgrid) or (reference_ax is None):
        return outdict

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                      np.arange(ylim[0], ylim[1]))
    Z = np.zeros(X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]


    outdict['XYZ'] = [X,Y,Z]

    return outdict

    # yunk = ax.plot_wireframe(X,Y,Z, color='k')
    
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # plt.show()

def get_plane_points(a,b,d,xlim=[-10,10],ylim=[-10,10]):
    """
    Generate [x,y,z] points in a plane by providing a,b, and d in the equation

    a*x + b*y + d = z

    and xlim and ylim

    """

    X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                      np.arange(ylim[0], ylim[1]))
    Z = np.zeros(X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r,c] = a * X[r,c] + b * Y[r,c] + d

    return X,Y,Z

def make_rot_matrix(vec0,vecWant):
    """
    Let M be the vector normal to your current plane, and N be the vector normal to the plane you want to rotate into.
    If M == N you can stop now and leave the original points unchanged.

    Ripped from https://stackoverflow.com/questions/9423621/3d-rotations-of-a-plane
    """
    normvec = vec0/np.linalg.norm(vec0)

    costheta = np.dot(normvec,vecWant)/(np.linalg.norm(normvec)*np.linalg.norm(vecWant))

    unitcross = np.cross(normvec,vecWant)
    unitcross = unitcross/np.linalg.norm(unitcross)

    x,y,z = unitcross

    c = costheta
    s = np.sqrt(1-c*c)
    C = 1-c
    rmat = np.matrix([[ x*x*C+c  ,  x*y*C-z*s,  x*z*C+y*s ],
                     [ y*x*C+z*s,  y*y*C+c  ,  y*z*C-x*s ],
                     [ z*x*C-y*s,  z*y*C+x*s,  z*z*C+c   ]])
    return rmat
