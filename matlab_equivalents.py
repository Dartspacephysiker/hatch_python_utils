import numpy as np

def sphere(n):
    theta = np.flip(np.pi*np.arange(n+1)/n)
    phi = 2*np.pi*np.arange(n+1)/n
    
    # T,P = np.meshgrid(theta,phi)
    P,T = np.meshgrid(phi,theta)
    x = np.sin(T)*np.cos(P)
    y = np.sin(T)*np.sin(P)
    z = np.cos(T)

    return x,y,z

