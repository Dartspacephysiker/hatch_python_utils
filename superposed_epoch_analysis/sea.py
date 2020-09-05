import numpy as np
import pandas as pd
from hatch_python_utils.earth.seasons import get_scaled_season_parameter,add_scaled_season_column

def screen_epochtime_db(dfepoch,
                        doScreenBySeason=False,
                        seasonOrMonths=None,
                        do_screen_by_tdiff=True,
                        mintdiffHours=60,
                        mintdiff_how_to_screen='keep_first',
                        data_latitudecol='alat110',
                        verbose=True):
    """
    screen_epochtime_db(**kws)
    
    ARGUMENTS
    =========
    dfepoch           (pandas.DataFrame): Dataframe with index == pd.DatetimeIndex for which to calculate relative epoch times
                      
    KEYWORDS          
    ========          
    doScreenBySeason              (bool): Shall I screen by season?
    seasonOrMonths         (list or str): How to do season screening? 
                                          Could be, e.g., either [3,4,9,10] (Mar, Apr, Sep, Oct), or "winter", or "Mar" for March equinox
    do_screen_by_tdiff            (bool): Do screen by temporal separation in hours
    mintdiffHours              (numeric): minimum separation, in hours, between adjacent sudden commencements
    mintdiff_how_to_screen         (str): One of "keep_first", "keep_last", or "keep_neither".

    NOTE: If seasonOrMonths = 

    SMH 2020-08-18
    Birkeland Centre for Space Science
    Universitetet i Bergen
    """
    
    ##############################
    # SAFETY CHECKS

    assert isinstance(dfepoch.index,pd.core.indexes.datetimes.DatetimeIndex),"dfepoch must have index of type pd.DatetimeIndex!"

    # Ensure valid input for screening by season, if any
    if doScreenBySeason:
    
        isnumericlist = isinstance(seasonOrMonths,list) and all([seas in np.arange(1,13) for seas in seasonOrMonths])
        isseason = isinstance(seasonOrMonths,str) and (seasonOrMonths[:3].lower() in ['mar','jun','sep','dec'])
        islocalseason = isinstance(seasonOrMonths,str) and (seasonOrMonths[:3].lower() in ['spr','sum','fal','win'])
        assertwarn = "Must have seasonOrMonths be one of [1..12,'mar','jun','sep','dec','spr','sum','fal','win']!"
        assert isnumericlist or isseason or islocalseason,assertwarn

    if do_screen_by_tdiff:
        assert np.isscalar(mintdiffHours),"mintdiffHours must be a scalar!"
        assert mintdiff_how_to_screen in ['keep_first','keep_last','keep_neither'],"mintdiff_how_to_screen must be one of ['keep_first','keep_last','keep_neither']"

    ##############################
    # SCREENING

    ##########
    # For starters, assume we keep all storms
    keepepoch = np.ones(dfepoch.shape[0],dtype=np.bool)        

    ##########
    # Screen so that separation between SCs is at least [mintdiffHours] hours
    if do_screen_by_tdiff:

        sc_tdiffs = np.diff(dfepoch.index).astype(np.int64)/1e9/3600  # hours
        sc_tdiffs = np.insert(sc_tdiffs,0,10000)

        if mintdiff_how_to_screen == 'keep_first':
            newscreen = keepepoch & (sc_tdiffs > mintdiffHours)
        elif mintdiff_how_to_screen in ['keep_last','keep_neither']:

            # sc_epochs = np.array([0,10,70,130,210])
            # sc_tdiffs = np.diff(sc_epochs)
            # sc_tdiffs = np.insert(sc_tdiffs,0,10000)

            # if mintdiff_how_to_screen == 'keep_first':
            #     keepheads = sc_tdiffs > mintdiffHours
            # else:
            newscreen = np.ones(len(sc_tdiffs),dtype=np.bool)

            for i,val in enumerate(sc_tdiffs > mintdiffHours): 
                # print(i,val)
                if not val:
                    if mintdiff_how_to_screen == 'keep_last':
                        # MIGHT BE THAT WE NEED TO DO SOMETHING SPECIAL WITH i==0 IN THIS CASE
                        newscreen[i-1] = False
                        newscreen[i] = True
                    elif mintdiff_how_to_screen == 'keep_neither':
                        # MIGHT BE THAT WE NEED TO DO SOMETHING SPECIAL WITH i==0 IN THIS CASE
                        newscreen[i-1] = False
                        newscreen[i] = False

        else:
            assert 2<0,"Huh?"

        nDropped = keepepoch.sum() - newscreen.sum()
        keepepoch = keepepoch & newscreen
        
        if verbose:
            print(f'mintdiff (hours)                 : {mintdiffHours}')
            print(f"(N dropped: {nDropped})")


    if doScreenBySeason:

        assert not islocalseason,"localseason should not make it here! To screen your epochDB by local season, use get_epoch_reltimes. There you can use, e.g., seasonOrMonths = 'fall' and make sure that your dataframe contains a latitude column (specified by KW data_latitudecol) so that get_epoch_reltimes() can set up the proper screening based on local season"

        if isnumericlist:

            print("Screening SCs by storm month, too. Alloweds:",SCmonths)
            dfepoch['month'] = dfepoch.index.month
            keepepoch = keepepoch & dfepoch['month'].isin(SCmonths)

        elif isseason:

            add_scaled_season_column(dfepoch,
                                     data_latitudecol=data_latitudecol,
                                     verbose=False)

            ctrpts = dict(mar=4.,
                          jun=1.,
                          sep=2.,
                          dec=3.)
            val0 = ctrpts[seasonOrMonths[:].lower()] - 0.5
            val1 = (ctrpts[seasonOrMonths[:].lower()] + 0.5) % 4
            print("Screening based on '{:s}' season ({:.2f} <= phi_s < {:.2f})".format(seasonOrMonths,val0,val1))

            if seasonOrMonths[:3].lower() == 'mar':
                newscreen = keepepoch & ((val0 <= dfepoch['scaledtidoffset']) | (dfepoch['scaledtidoffset'] < val1))
            else:
                newscreen = keepepoch & ((val0 <= dfepoch['scaledtidoffset']) & (dfepoch['scaledtidoffset'] < val1))

            nDropped = keepepoch.sum() - newscreen.sum()
            keepepoch = keepepoch & newscreen
        
            if verbose:
                print(f'season                 : {seasonOrMonths}')
                print(f"(N dropped: {nDropped})")

    print(f"N epochs to keep: {keepepoch.sum()}")
    dfepoch = dfepoch[keepepoch]
    
    return dfepoch


def get_epoch_reltimes(dfdata,dfepoch,
                       dtepochcol='dtEpoch',
                       data_latitudecol='alat110',
                       befaftEpochhours=np.array([60,60]),
                       doScreenBySeason=False,
                       seasonOrMonths=None,
                       mintdiffHours=60,
                       verbose=True):
    """
    get_epoch_reltimes(dfdata,dfepoch,**kws)
    
    ARGUMENTS
    =========
    dfdata            (pandas.DataFrame): Dataframe with type(index) == pd.DatetimeIndex for which to calculate relative epoch times
    dfepoch           (pandas.DataFrame): Dataframe with epochs (and type(index) == pd.DatetimeIndex) 
                      
    KEYWORDS          
    ========          
    befaftEpochhours     (np.array, size 2): Two-element array for minimum and maximum epoch time
    dtepochcol                        (str): Name of column in which to store epoch times in DataFrame dfdata
    doScreenBySeason                 (bool): Shall I screen by season?
    seasonOrMonths            (list or str): How to do season screening? 
                                             Could be, e.g., either [3,4,9,10] (Mar, Apr, Sep, Oct), or "winter", or "Mar" for March equinox
    data_latitudecol                  (str): Only needed if screening by season, and seasonOrMonths is one of ['winter','spring','summer','fall']


    SMH 2020-08-18
    Birkeland Centre for Space Science
    Universitetet i Bergen
    """

    assert isinstance(dfdata.index,pd.core.indexes.datetimes.DatetimeIndex),"dfdata must have index of type pd.DatetimeIndex!"

    # Make sure we can do what is asked first
    dolocalseasonscreening = False
    if doScreenBySeason:
    
        isnumericlist = isinstance(seasonOrMonths,list) and all([seas in np.arange(1,13) for seas in seasonOrMonths])
        isseason = isinstance(seasonOrMonths,str) and (seasonOrMonths[:3].lower() in ['mar','jun','sep','dec'])
        islocalseason = isinstance(seasonOrMonths,str) and (seasonOrMonths[:3].lower() in ['spr','sum','fal','win'])
        assertwarn = "Must have seasonOrMonths be one of [1..12,'mar','jun','sep','dec','spr','sum','fal','win']!"
        assert isnumericlist or isseason or islocalseason,assertwarn

        dolocalseasonscreening = islocalseason

    # Reset any existing calculation
    dfdata.loc[:,dtepochcol] = np.nan

    # Convert befaftEpochhours to two-element array of proper timedeltas
    befaftSC = np.array([pd.Timedelta(f'{befaftEpochhours[0]} hours').to_numpy(),
                         pd.Timedelta(f'{befaftEpochhours[1]} hours').to_numpy()])
    

    # Don't waste time with epochs occurring before or after timespan of dfdata
    dfepochnew = dfepoch[(dfepoch.index >= (dfdata.index[0]-befaftSC[0])) & (dfepoch.index <= (dfdata.index[-1]+befaftSC[1]))]

    # magind = (np.diff(dfepochnew.index).astype(np.int64)/1e9/3600).argmin()
    # # magind = np.where(np.diff(dfepochnew.index).astype(np.int64)/1e9/3600/24 < 0)[0]
    # dfepochnew.iloc[magind-1],dfepochnew.iloc[magind],dfepochnew.iloc[magind+1]
    
    if dolocalseasonscreening:
        
        assert data_latitudecol in dfdata.columns

        seasString = seasonOrMonths[:3].lower()
        # In this case, we need to treat each hemisphere separately
        dfepochnewN = dfepochnew.copy()
        dfepochnewS = dfepochnew.copy()
        
        Sseasmapper = dict(spr='sep',
                           sum='dec',
                           fal='mar',
                           win='jun')
        Nseasmapper = dict(spr='mar',
                           sum='jun',
                           fal='sep',
                           win='dec')
        
        print("Using {:s} for SH in {:s}".format(Sseasmapper[seasString],seasString))
        print("Using {:s} for NH in {:s}".format(Nseasmapper[seasString],seasString))
        dfepochnewS = screen_epochtime_db(dfepochnewS,
                                       doScreenBySeason=doScreenBySeason,
                                       seasonOrMonths=Sseasmapper[seasString],
                                       data_latitudecol=data_latitudecol,
                                       mintdiffHours=mintdiffHours,
                                       verbose=True)
        
        dfepochnewN = screen_epochtime_db(dfepochnewN,
                                       doScreenBySeason=doScreenBySeason,
                                       seasonOrMonths=Nseasmapper[seasString],
                                       data_latitudecol=data_latitudecol,
                                       mintdiffHours=mintdiffHours,
                                       verbose=True)

        ##############################
        # Get epoch times, SH
    
        for hemi,dfepochnew in zip(['NH','SH'],[dfepochnewN,dfepochnewS]):
            if hemi == 'NH':
                hemiInds = dfdata[data_latitudecol] > 0
            elif hemi == 'SH':
                hemiInds = dfdata[data_latitudecol] < 0

            print("===================")
            print("Doing {:s} ({:d} inds)".format(hemi,hemiInds.sum()))

            for index in dfepochnew.index:
                # print(index)
                indshere = ((index-dfdata.index) <= befaftSC[0]) & ((dfdata.index - index) <= befaftSC[1]) & hemiInds
                nIndsHere = indshere.sum()

                if nIndsHere > 0:
                    dfdata.loc[indshere,dtepochcol] = (dfdata.index[indshere]-index).to_numpy().astype(np.int64)/1e9/3600
                if verbose:
                    print(index,nIndsHere)

    else:

        dfepochnewUse = dfepochnew.copy()
        if doScreenBySeason:

            dfepochnewUse = screen_epochtime_db(dfepochnewUse,
                                          doScreenBySeason=doScreenBySeason,
                                          seasonOrMonths=seasonOrMonths[:3].lower(),
                                          mintdiffHours=mintdiffHours,
                                          verbose=True)

        ##############################
        # Get epoch times
    
        for index in dfepochnewUse.index:
            # print(index)
            indshere = ((index-dfdata.index) <= befaftSC[0]) & ((dfdata.index - index) <= befaftSC[1])
            nIndsHere = indshere.sum()
        
            if nIndsHere > 0:
                dfdata.loc[indshere,dtepochcol] = (dfdata.index[indshere]-index).to_numpy().astype(np.int64)/1e9/3600
            if verbose:
                print(index,nIndsHere)


def get_epochs_with_data(dfdata,dfepoch,
                         befaftEpochhours=np.array([60,60]),
                         verbose=True):
    """
    get_epochs_with_data(dfdata,dfepoch,**kws)
    
    ARGUMENTS
    =========
    dfdata            (pandas.DataFrame): Dataframe with type(index) == pd.DatetimeIndex for which to calculate relative epoch times
    dfepoch           (pandas.DataFrame): Dataframe with epochs (and type(index) == pd.DatetimeIndex) 
                      
    KEYWORDS          
    ========          
    befaftEpochhours     (np.array, size 2): Two-element array for minimum and maximum epoch time


    SMH 2020-08-31
    Birkeland Centre for Space Science
    Universitetet i Bergen
    """

    assert isinstance(dfdata.index,pd.core.indexes.datetimes.DatetimeIndex),"dfdata must have index of type pd.DatetimeIndex!"

    # Convert befaftEpochhours to two-element array of proper timedeltas
    befaftSC = np.array([pd.Timedelta(f'{befaftEpochhours[0]} hours').to_numpy(),
                         pd.Timedelta(f'{befaftEpochhours[1]} hours').to_numpy()])
    

    havedata = np.zeros(dfepoch.shape[0],dtype=np.bool)
    firsttime,lasttime = dfdata.index[0],dfdata.index[-1]
    for i,index in enumerate(dfepoch.index):

        if index < (firsttime-befaftSC[0]):
            # Save time, don't bother calculating
            continue
        elif index > (lasttime+befaftSC[1]):
            # Save time, don't bother calculating
            continue
        else:
            # print(index)
            indshere = ((index-dfdata.index) <= befaftSC[0]) & ((dfdata.index - index) <= befaftSC[1])
            nIndsHere = indshere.sum()
        
            if nIndsHere > 0:
                havedata[i] = True
            
    return havedata



