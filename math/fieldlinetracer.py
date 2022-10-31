import numpy as np

class FLTRACE_CART(object):

    def __init__(self,fx,fy,fz,
                 r0,
                 ds=0.1,smax=30,
                 xlim=(-np.inf,np.inf),
                 ylim=(-np.inf,np.inf),
                 zlim=(-np.inf,np.inf),
                 verbose=False,
                 DEBUG=False,
    ):
        """
        Cartesian field-line tracer

        Inputs
        ======
        fx   : Function that takes arguments (x,y,z) and returns the x component of the vector field
        fy   : Function that takes arguments (x,y,z) and returns the y component of the vector field
        fz   : Function that takes arguments (x,y,z) and returns the z component of the vector field
             
        r0   : Initial position(s). If multiple starting positions the shape of r0 must be (N,3), where
               N is the number of starting positions

        ds   : Step size
        smax : Max total distance to trace

        xlim : Boundaries for tracing in x direction (similar for ylim, zlim) 
               DEFAULT: No limit

        Example
        =======

        # from hatch_python_utils.math.fieldlinetracer import FLTRACE_CART
        import importlib
        from hatch_python_utils.math import fieldlinetracer
        importlib.reload(fieldlinetracer)

        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        plt.ion()

        def fx(x,y,z):
            return np.sin(y)

        def fy(x,y,z):
            return np.sin(x)

        def fz(x,y,z):
            return 0

        r0 = np.array([[-1,-1,0],[-1,1,0],[1,-1,0],[1,1,0]],dtype=np.float64)
        
        # ftracer = FLTRACE_CART(fx,fy,fz,r0,ds=0.01)
        ftracer = fieldlinetracer.FLTRACE_CART(fx,fy,fz,r0,ds=0.01)
        ftracer.trace(verbose=True,DEBUG=True)

        fig,ax = plt.subplots(1,1)
        for i in range(ftracer.NLines):
            ax.plot(ftracer.tracedicts[i]['r'][:,0],ftracer.tracedicts[i]['r'][:,1],marker='o')

        2022/10/12 SMH 
        """
        
        ####################
        # Put all arguments into this object
        
        assert callable(fx) and callable(fy) and callable(fz),"fx, fy, fz must be callable functions!"

        self.fx = fx
        self.fy = fy
        self.fz = fz

        # Vector, vector magnitude, and unit vector 
        self.fvec = lambda x,y,z: np.array([self.fx(x,y,z),self.fy(x,y,z),self.fz(x,y,z)])
        self.fmag = lambda x,y,z: np.sqrt(self.fx(x,y,z)**2+self.fy(x,y,z)**2+self.fz(x,y,z)**2)
        self.funit = lambda x,y,z: self.fvec(x,y,z) / self.fmag(x,y,z)

        if not isinstance(r0,np.ndarray):
            r0 = np.array(r0)

        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim

        # Get number of starting points, make sure r0 has shape (N,3) regardless
        if r0.ndim == 1:
            assert r0.size == 3,"If providing single vector as starting point r0, length must be three!"
            self.NLines = 1
            r0 = r0[np.newaxis,:]
        elif r0.ndim == 2:
            assert r0.shape[1] == 3,"If providing N starting points, shape must be (N,3)!"
            self.NLines = r0.shape[0]

        # Make sure all starting points are within limits
        for i in range(self.NLines):
            assert self.__iswithinlim(r0[i,:])

        self.r0 = np.float64(r0)

        self.ds = ds
        self.smax = smax
        self.maxnsteps = int((self.smax*2)/self.ds)

        self.tracedicts = [self.__make_trace_dict(self.r0[i,:],DEBUG=DEBUG) for i in range(self.NLines)]


    def trace(self,verbose=False,
              DEBUG=False,
    ):
        
        for i in range(self.NLines):
            self.__trace_single(i,verbose=verbose,
                                DEBUG=False)

        return self.tracedicts


    def __trace_single(self,i,
                       verbose=False,
                       DEBUG=False,
    ):

        tdict = self.tracedicts[i]

        if DEBUG:
            print(f"Performing tracing for field line #{i:03d}")

        n = 0
        r = tdict['r0']
        tdict['r'][n,:] = r
        while n < self.maxnsteps-1:

            if DEBUG:
                print("")
                print(f"*** Step {n:03d}/{self.maxnsteps:03d}***")
                print(f"r = <{r[0]:8.2f},{r[1]:8.2f},{r[2]:8.2f}>")

            if not self.__iswithinlim(r):
                break

            k1 = self.funit(*r)
            k2 = self.funit(*(r+self.ds * k1/2))
            k3 = self.funit(*(r+self.ds * k2/2))
            k4 = self.funit(*(r+self.ds * k3  ))
            
            rnext = r+self.ds/6 * (k1+2*k2+2*k3+k4)

            tdict['r'][n+1,:] = rnext

            if np.any(~np.isfinite(rnext)):
                print("Tracing from r0=<",tdict['r0'],"> stopped at r=<",tdict['r'][n,:],f"> after N={n} steps; encountered NaNs!")
                break

            if DEBUG:
                tdict['k1'][n,:] = k1
                tdict['k2'][n,:] = k2
                tdict['k3'][n,:] = k3
                tdict['k4'][n,:] = k4

            r = rnext
            n += 1

        self.tracedicts[i] = tdict


    def __make_trace_dict(self,r0,DEBUG=False):
        emptyvector = np.zeros((self.maxnsteps,3))
        emptyscalar = np.zeros(self.maxnsteps)
        tdict = {"r0":r0,
                 "r" :emptyvector,
                
        }

        if DEBUG:
            tdict['k1'] = emptyvector.copy()
            tdict['k2'] = emptyvector.copy()
            tdict['k3'] = emptyvector.copy()
            tdict['k4'] = emptyvector.copy()

        return tdict

    
    def __iswithinlim(self,r):
        x,y,z = r
        withinlimx = (self.xlim[0] <= x <= self.xlim[1])
        withinlimy = (self.ylim[0] <= y <= self.ylim[1])
        withinlimz = (self.zlim[0] <= z <= self.zlim[1])
        if not withinlimx and withinlimy and withinlimz:
            print(f"r = <{x},{y},{z}> out of bounds:")
            if not withinlimx:
                print(f"Not true that xlim[0] = {self.xlim[0]} <= x <= {self.xlim[1]} = xlim[1]")
            if not withinlimy:
                print(f"Not true that ylim[0] = {self.ylim[0]} <= y <= {self.ylim[1]} = ylim[1]")
            if not withinlimz:
                print(f"Not true that zlim[0] = {self.zlim[0]} <= z <= {self.zlim[1]} = zlim[1]")
            
            return False

        return True

if __name__ == "__main__":
    import importlib
    from hatch_python_utils.math import fieldlinetracer
    importlib.reload(fieldlinetracer)

    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.ion()

    def fx(x,y,z):
        return np.sin(y)

    def fy(x,y,z):
        return np.sin(x)

    def fz(x,y,z):
        return 0

    r0 = np.array([[-1,-1.2,0],[-1,-0.8,0],[-1,-1.,0],
                   [-1, 1.2,0],[-1, 0.8,0],[-1, 1.,0],
                   [ 1,-1.2,0],[ 1,-0.8,0],[ 1,-1.,0],
                   [1,1.2,0],[1,0.8,0],[1.,1.,0]],dtype=np.float64)
        
    colors = [*(['C0']*3),*(['C1']*3),*(['C2']*3),*(['C3']*3)]
    lstyles = [*(['-','--',':']*4)]
    lwidths = [*([1,2,3]*4)]

    # ftracer = FLTRACE_CART(fx,fy,fz,r0,ds=0.01)
    ftracer = fieldlinetracer.FLTRACE_CART(fx,fy,fz,r0,ds=0.05,smax=5)
    ftracer.trace(verbose=True,DEBUG=True)

    fig,ax = plt.subplots(1,1)
    for i in range(ftracer.NLines):
        tdict = ftracer.tracedicts[i]
        r0 = tdict['r0']
        ax.text(r0[0],r0[1],str(i))
        ax.scatter([r0[0]],[r0[1]],
                   marker='o',
                   color=colors[i],
                   label='Start' if i==0 else None)
        ax.scatter([tdict['r'][-1,0]],[tdict['r'][-1,1]],
                   marker='*',
                   color=colors[i],
                   label='Stop' if i==0 else None)
        ax.plot(tdict['r'][:,0],tdict['r'][:,1],
                # marker='.',
                linestyle=lstyles[i],
                color=colors[i],
                lw=lwidths[i],
                label=f'#{i:02d}: <{r0[0]:.1f},{r0[1]:.1f},{r0[2]:.1f}>')
    ax.legend()
    ax.set_aspect('equal')


