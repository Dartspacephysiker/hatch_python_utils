import numpy as np

from scipy.special import logit

from matplotlib import pyplot as plt
plt.ion()

import hatch_python_utils.math.continuous_cecs
importlib.reload(hatch_python_utils.math.continuous_cecs)
from hatch_python_utils.math.continuous_cecs import BsplineCECS,get_rho,get_rhohat,get_phi,get_phihat

########################################
## Example 0: Make sure unit vectors are correct

xcecs = np.array([-0.5,0.5])
ycecs = np.array([0.,0.])
zcecs = np.array([0.,0.])

x = np.array([-1.,-0.5,0,0.5,1.,-1.,0,1.,-1.,-0.5,0.,0.5,1.])
y = np.array([-0.5,-0.5,-0.5,-0.5,-0.5,0.,0.,0.,0.5,0.5,0.5,0.5,0.5])


# x = np.array([-0.25,0.,0.25,-0.25,0.,0.25,-0.25,0.,0.25])
# y = np.array([-0.25,-0.25,-0.25,0,0,0,0.25,0.25,0.25])
z = np.zeros_like(y)

dfcecs = BsplineCECS(RBFX,RBFY,knotv,
                   polyorder=3,
                   current_type='divergence_free')


rho = get_rho(x,y,xcecs,ycecs)
phi = get_phi(x,y, xcecs, ycecs)
rhohat = get_rhohat(x, y, xcecs, ycecs)
phihat = get_phihat(x, y, xcecs, ycecs)


fig,ax = plt.subplots(1,1)

nMeas = len(x)
nCECS = len(xcecs)
for i in range(nMeas):
    ax.scatter([x[i]],[y[i]])
    ax.text(x[i],y[i],str(i),color='C'+str(i))

for i in range(nCECS):
    ax.scatter([xcecs[i]],[ycecs[i]],color='gray',marker='*')
ax.set_aspect('equal')

########################################
## Example 1: B-spline functions and their integrals

maxAltitude = 700

i_maxalt = np.where(z <= maxAltitude)[0].size
z = np.linspace(70,1200,33)
z = z[:i_maxalt]


## Select locations of RBFs
rbfdx,rbfdy = 40,40
rbfx = np.arange(-280,281,rbfdx)
rbfy = np.arange(-280,281,rbfdy)
rbfz = z

RBFX, RBFY = np.meshgrid(rbfx,rbfy,indexing='ij')

##Pick knot vector
# Select location of knots using logit function
k = 3                           # Polynomial order
Nknots = 15+2*k                 # Number of knots

tmpx = np.linspace(0.6,0.9999,Nknots)

# Select location of knots using logit function
hmin,hmax = 20,1500
heights = np.linspace(hmin,hmax,3001)

knotv = logit(tmpx)                 # Knots
knotv = (knotv - knotv.min()) / (knotv.max()-knotv.min()) * (hmax-hmin)+hmin  # Scale knots between hmin and hmax

dfcecs = BsplineCECS(RBFX,RBFY,knotv,
                   polyorder=3,
                   current_type='divergence_free')

cfcecs = BsplineCECS(RBFX,RBFY,knotv,
                   polyorder=3,
                   current_type='curl_free')


rbfpeakz = dfcecs.Bspl_func_peakheights


zeval = np.linspace(*cfcecs.internalknotinterval,1001)
showj = [0,2,4,6,9,12,15]

fig, ax = plt.subplots(1,1)
ax2 = ax.twinx()
for ij,j in enumerate(showj):
    ax.plot(zeval,cfcecs.bsplinefuncs[j](zeval),color='C'+str(ij),label='Func' if ij == 0 else None)
    integvals = np.array([cfcecs.bsplineintegfuncs[j](zer) for zer in zeval])
    ax2.plot(zeval,integvals,linestyle='--',linewidth=2,color='C'+str(ij),label='Integ' if ij == 0 else None)
