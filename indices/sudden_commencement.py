import pandas as pd
import numpy as np

def load_sc_database(minQualCode=0,
                     minNStationswithMinQualCode=0,
                     SCType=None,
                     SCType__beGenerousBefore2005=True,
                     verbose=True):
    """
    dfsc = load_sc_database(**kws)

    KEYWORDS          
    ========          
    minQualCode                 (scalar): Number between 1 and 3. (1 means "probably SC," and 3 means "very clearly SC." 
                                          See documentation at http://isgi.unistra.fr/data_download.php, and in metadata of SC lists downloaded from there.
    minNStationswithMinQualCode (scalar): Min number of stations (max 5) with quality code >= minQualCode.
    SCType                         (str): Impulse followed by storm ('SSC') or impulse without consequent geomagnetic storm activity ('SI')
    SCType__beGenerousBefore2005  (bool): Storm identification did not happen before 2006, so if this is true, don't drop these storms

    SMH 2020-08-21
    Birkeland Centre for Space Science
    Universitetet i Bergen
    """

    scdir = '/SPENCEdata/Research/database/sudden_commencement/'
    # output info
    outfile = 'SuddenCommencements__1999-2019.parquet'

    ##############################
    # SAFETY CHECKS

    assert minQualCode <= 3,"Max quality code is 3! Reduce 'minQualCode'"
    assert minNStationswithMinQualCode <= 5,"Max number of stations is 5! Reduce 'minNStationswithMinQualCode'"

    if (minNStationswithMinQualCode > 0):
        assert minQualCode > 0,"Invalid! You've set a threshold on minNStationswithMinQualCode, but minQualCode == 0!"

    if (minQualCode > 0):
        assert minNStationswithMinQualCode > 0,"Invalid! You've set a threshold on minimum quality code, but minNStationswithMinQualCode == 0!"
        
    if SCType is not None:
        assert SCType in ['SSC','SI'],"Valid SC types are either 'SSC' or 'SI'!"


    ##############################
    # LOAD DB
    dfsc = pd.read_parquet(scdir+outfile)

    ##############################
    # SCREENING

    ##########
    # For starters, assume we keep all storms
    keepepoch = np.ones(dfsc.shape[0],dtype=np.bool)        

    ##########
    # Making this separate "if" statement so that the code is clear
    if (minNStationswithMinQualCode > 0) and (minQualCode > 0):

        newscreen = keepepoch & (dfsc['n>='+str(minQualCode)] >= minNStationswithMinQualCode)
        nDropped = keepepoch.sum() - newscreen.sum()
        keepepoch = keepepoch & newscreen

        if verbose:
            print(f'minNStations w/ quality code >= {minQualCode}: {minNStationswithMinQualCode}')
            print(f"(N dropped: {nDropped})")

    if SCType is not None:
        
        if SCType__beGenerousBefore2005:
            bef2006 = dfsc.index.year <= 2005
            newscreen = keepepoch & (( (dfsc['SCType'] == SCType) & ~bef2006) | bef2006)
        else:
            newscreen = keepepoch & (dfsc['SCType'] == SCType)
        nDropped = keepepoch.sum() - newscreen.sum()
        keepepoch = keepepoch & newscreen

        if verbose:
            print(f'SCType                           : {SCType}',end='')
            if SCType__beGenerousBefore2005:
                print(' (being generous to SCs before 2006)')
            else:
                print("")

            print(f"(N dropped: {nDropped})")
        
    return dfsc
