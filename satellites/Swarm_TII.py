import numpy as np
# from hatch_python_utils.math.vectors import dotprod


def calculate_satellite_frame_vectors_in_NEC_coordinates(df,
                                                         vCol=['VsatN','VsatE','VsatC'],
                                                         assert_orthogonality=True):
    """
    Use the Swarm velocity components in NEC coordinates to calculate satellite-frame unit vectors in NEC coordinates

    ## Definition of satellite-track measurements according to "EFI TII Cross-Track Flow Data Release Notes" (Doc. no: SW-RN-UOC-GS-004, Rev: 6)
    "The coordinate system for the satellite-track measurements (along- and cross-track) is a right-handed orthonormal system
    * Defined in a frame of reference co-rotating with the Earth and having 
       * x in the direction of the satellite velocity, 
       * y perpendicular to x and horizontally to the right when facing the direction of motion, and 
       * z approximately downward completing the triad."

    "Measurements may be transformed into the north-east-centre (NEC) system using the supplied satellite velocity vector in NEC coordinates as a reference."

    2020-10-02
    SMH
    """

    def dotprod(a,b):
        return np.einsum('ij,ij->i',a,b)

    # Make sure we have all columns we need
    assert all(vCo in df.columns for vCo in vCol)

    vUnitCol = [vCo+'hat' for vCo in vCol]

    # Regne ut størrelsen til hastighetsvektorer
    df['VsatMag'] = np.sqrt((df[vCol]**2).sum(axis=1))

    # Calculate magnitude of "horizontal" (i.e., not-vertical) component of satellite velocity vector
    df['VsatHoriz'] = np.sqrt((df[['VsatN','VsatE']]**2).sum(axis=1))

    # The condition $|v_c|/v \ll \sqrt{v_n^2+v_e^2}/v$ must be fulfilled if the statement "z is approximately downward" is to be true
    # assert this condition
    assert (np.abs(df['VsatC'])/df['VsatHoriz']).max() < 0.01

    # Regne ut komponentene til hastighets-enhetsvektor
    for vCo,vUnitCo in zip(vCol,vUnitCol):
        df[vUnitCo] = df[vCo]/df['VsatMag']

    # Forsikre oss om at enhetsvektorene har faktisk størrelse 1 :)
    assert np.all(np.isclose(1,(df[vUnitCol]**2.).sum(axis=1)))

    #Da ifølge dokumentasjonen ...

    ########################################
    # i. Definere x-enhetsvektoren i NEC-koordinater
    
    xN = df['VsatNhat']
    xE = df['VsatEhat']
    xC = df['VsatChat']
    
    
    ########################################
    # ii. Definere y-enhetsvektoren i NEC-koordinater
    
    # Regne ut fortegn til yE.
    #Da y-enhetsvektoren er til høyre når man ser i bevegelsesretning, den peker østover (dvs yE > 0) når romfartøy går nordover, og vest (dvs yE < 0) når romfartøy går sørover
    yESign = np.int64(xN >= 0)*2-1

    yE = yESign / np.sqrt( (xE**2 / xN**2) + 1)
    yN = - xE / xN * yE
    yC = yE * 0
    
    # Renorm just to make sure magnitudes are 1
    yMag = np.sqrt(yN**2 + yE**2 + yC**2)
    yN,yE,yC = yN / yMag,yE / yMag,yC / yMag
    
    ########################################
    # iii. Definere z-enhetsvektoren i NEC-koordinater
    zN, zE, zC = np.cross(np.vstack([xN,xE,xC]).T,np.vstack([yN,yE,yC]).T).T

    # Renorm just to make sure magnitudes are 1
    zMag = np.sqrt(zN**2+zE**2+zC**2)
    zN,zE,zC = zN / zMag,zE / zMag,zC / zMag
    
    xhat = np.vstack([xN,xE,xC])
    yhat = np.vstack([yN,yE,yC])
    zhat = np.vstack([zN,zE,zC])

    ########################################
    # Assert orthogonality
    
    if assert_orthogonality:
        # Assert xhat-prikk-yhat == 0
        assert np.max(np.abs(dotprod(xhat.T,yhat.T))) < 0.00001
        # Assert xhat-prikk-zhat == 0
        assert np.max(np.abs(dotprod(xhat.T,zhat.T))) < 0.00001
        # Assert yhat-prikk-zhat == 0
        assert np.max(np.abs(dotprod(yhat.T,zhat.T))) < 0.00001

    return xhat,yhat,zhat
