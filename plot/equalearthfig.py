import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs

from hatch_python_utils.earth.utils import get_noon_longitude
from hatch_python_utils.plot import colormaps as hCM
from pysymmetry.visualization.polarsubplot import Polarsubplot

class EqualEarthFig(object):

    def __init__(self,fig,gs,tStamp,totrows=18,totcols=32):

        # ncols=1
        # colwidth = 1
        cbcolfrac = 0.02
        
        polplotcols = int(np.ceil(totrows/2))
        
        cbcols = int(np.ceil(cbcolfrac*totcols))
        
        # earth panel
        earthpanelrowfrac = 0.7
        earthpanelcols    = totcols-cbcols-polplotcols
        earthpanelrows = int(np.ceil(earthpanelrowfrac*totrows))
        earthstartc,earthendc = 0,earthpanelcols
        earthstartr,earthendr = 0,earthpanelrows
        
        # polarplot panels
        pol1startc,pol1endc = earthpanelcols,earthpanelcols+polplotcols
        pol1startr,pol1endr = 0,polplotcols
        pol2startc,pol2endc = earthpanelcols,earthpanelcols+polplotcols
        pol2startr,pol2endr = polplotcols,polplotcols*2
        
        # colorbar
        cbstartc,cbendc = pol1endc,totcols
        cbstartr,cbendr = 0,totrows
        
        # t series
        tseriesrowfrac = (1-earthpanelrowfrac)/2.
        tseriesrows = int(np.ceil(tseriesrowfrac*totrows))
        
        ts1startc,ts1endc = earthstartc,earthendc
        ts2startc,ts2endc = earthstartc,earthendc
        ts1startr,ts1endr = earthendr,earthendr+tseriesrows
        ts2startr,ts2endr = earthendr+tseriesrows,earthendr+tseriesrows*2
        
        print(f"cbcols            : {cbcols}")
        print(f"polplotcols       : {polplotcols}")
        print(f"earthpanelrows    : {earthpanelrows}")
        print(f"earthpanelcols    : {earthpanelcols}")
        print(f"tseriesrows       : {tseriesrows}")
        
        self.projection,self.projectionName = ccrs.EqualEarth(central_longitude=get_noon_longitude(tStamp,verbose=False)[0]),'EqualEarth'

        self.axearth = fig.add_subplot(gs[earthstartr:earthendr,earthstartc:earthendc],projection=self.projection)
        self.axts1 = fig.add_subplot(gs[ts1startr:ts1endr,ts1startc:ts1endc])
        self.axts2 = fig.add_subplot(gs[ts2startr:ts2endr,ts2startc:ts2endc],sharex=self.axts1)
        self.axpol1 = fig.add_subplot(gs[pol1startr:pol1endr,pol1startc:pol1endc])
        self.axpol2 = fig.add_subplot(gs[pol2startr:pol2endr,pol2startc:pol2endc])
        self.cax = fig.add_subplot(gs[cbstartr:cbendr,cbstartc:cbendc])


    def make_polar_panel(self,tmpds,showquant,
                         cb_do_scientific=False,
                         cb_powlims=(4,4),
                         add_Weimer=False,
                         **contourkw):
    
        if add_Weimer:
            assert 2<0,"Not implemented!!!"

        reqdkw = ['norm','levels','cmap']
        for kw in reqdkw:
            assert kw in contourkw,f"Must provide '{kw}' as kw in call to make_polar_panel!"

        showindsN = (tmpds.mlat.values >= 45) & (tmpds.mlat.values < 90)
        showindsS = (tmpds.mlat.values <= -45) & (tmpds.mlat.values > -90)
    
        paxN = Polarsubplot(self.axpol1)
        paxS = Polarsubplot(self.axpol2)
    
        # contN = paxN.contourf(tmpds.mlat.values[showindsN],
        #                       tmpds.mlt.values[showindsN],
        #                       tmpds[showquant].values[showindsN],
        #                       levels=levels,
        #                       cmap=cmap)
        # contS = paxS.contourf(tmpds.mlat.values[showindsS],
        #                       tmpds.mlt.values[showindsS],
        #                       tmpds[showquant].values[showindsS],
        #                       levels=levels,
        #                       cmap=cmap)
        contN = paxN.contourf(tmpds.mlat.values[showindsN],
                              tmpds.mlt.values[showindsN],
                              tmpds[showquant].values[showindsN],
                              **contourkw)
        contS = paxS.contourf(tmpds.mlat.values[showindsS],
                              tmpds.mlt.values[showindsS],
                              tmpds[showquant].values[showindsS],
                              **contourkw)
    
        titN = self.axpol1.set_title('North')
        titS = self.axpol2.set_title('South')
        
        # if add_Weimer:
        #     dickie = get_Weimers_NS(allds,dateind)
        #     paxN.plotpins(dickie['N']['mlat'],dickie['N']['mlt'],dickie['N']['vN'],dickie['N']['vE'],SCALE=750)
        #     paxS.plotpins(dickie['S']['mlat'],dickie['S']['mlt'],dickie['S']['vN'],dickie['S']['vE'],SCALE=750)
    
    
        # plt.tight_layout()
        
        
        sfmt = mpl.ticker.ScalarFormatter(useMathText=True) 
        sfmt.set_powerlimits((0, 0))
        # cb = plt.colorbar(mpl.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,format='%.1e')

        # norm = mpl.colors.BoundaryNorm(contourkw['levels'],
        #                                ncolors=contourkw['cmap'].N,
        #                                clip=True)

        cb = mpl.pyplot.colorbar(
            mpl.cm.ScalarMappable(norm=contourkw['norm'],
                                  cmap=contourkw['cmap']),
            cax=self.cax,
            format=sfmt)
        _ = cb.set_label(tmpds[showquant].attrs['long_name']+' ['+tmpds[showquant].attrs['units']+']')
    
        if cb_do_scientific:
            cb.formatter.set_scientific(True)
            cb.formatter.set_powerlimits(cb_powlims)
            cb.update_ticks()
    
        return cb
    

    def plot_earth_panel(self,lon,lat,quant,**contourkw):

        # contN = self.axearth.contourf(tmpds['Lon'],
        #                               tmpds['Lat'],
        #                               tmpds[showquant],
        #                               transform=ccrs.PlateCarree(),
        #                               **contourkw)

        contN = self.axearth.contourf(lon,lat,quant,
                                      transform=ccrs.PlateCarree(),
                                      **contourkw)

        _ = self.axearth.coastlines()
        _ = self.axearth.set_global()

        return contN


    def plot_tseries_panel(self,statN,statS,
                           tstampvline=None,
                           NHcol='dodgerblue',
                           SHcol='lightcoral'):
    
        pN = self.axts1.plot(statN,color=NHcol)
        pS = self.axts1.plot(statS,color=SHcol)

        _ = [xlab.set_visible(False) for xlab in self.axts1.xaxis.get_ticklabels()]

        if tstampvline is not None:
            pvline = self.axts1.axvline(mpl.dates.date2num(tstampvline),color='gray')
            return (pN,pS,pvline)

        else:
            return (pN,pS)
        
    def plot_omni_panel(self,dfOMNI,
                        omniquant='SymH',
                        tstampvline=None,
                        t0=None,
                        t1=None,
                        ylabel=None):
    
        if t0 is None:
            assert isinstance(dfOMNI.index,pd.DatetimeIndex)
            t0 = pd.Timestamp(dfOMNI.index.values[0])
    
        if t1 is None:
            assert isinstance(dfOMNI.index,pd.DatetimeIndex)
            t1 = pd.Timestamp(dfOMNI.index.values[1])
    
        if ylabel is None:
            ylabel= omniquant
    
        pOMNI = self.axts2.plot(dfOMNI[omniquant])
        _ = self.axts2.set_ylabel(ylabel)
    
        self.axts2.set_xlim((mpl.dates.date2num(t0),mpl.dates.date2num(t1)))

        if tstampvline is not None:
            pvline = self.axts2.axvline(mpl.dates.date2num(tstampvline),color='gray')
            return (pOMNI,pvline)

        else:
            return pOMNI

    
