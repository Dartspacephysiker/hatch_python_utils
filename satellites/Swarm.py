import ftplib
import os.path
import fnmatch
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from hatch_python_utils import omni_utils as hOMNI
from hatch_python_utils.earth import coordinates as hCoord
from hatch_python_utils.earth.bearing import mlt_bearing

import sys
needdir = '/SPENCEdata/Research/sandbox_and_journals/journals/Swarm/'
if needdir not in sys.path:
    sys.path.insert(0, needdir)
import swarmProcHelper as sPH
import conjunctions as swarmC

def load_Swarm_ephemeris(*args,**kws):
    """
    df = load_Swarm_ephemeris([sat,yr0,yr1])

    KWS
    ====
    dotiming (default False)  : Show timing stats for loading stuff
    """

    defaultArgs = dict(sat='A',yr0=2013,yr1=2019)
    defaultKWs = dict(dotiming=False)

    sat = args[0] if (len(args) >= 1) else defaultArgs['sat']
    yr0 = args[1] if (len(args) >= 2) else defaultArgs['yr0']
    yr1 = args[2] if (len(args) >= 3) else defaultArgs['yr1']

    dotiming = kws['dotiming'] if ('dotiming' in kws.keys()) else defaultKWs['dotiming']

    from hatch_python_utils.hatch_utils import get_basedir 

    basedir, isColtrane = get_basedir()
    outdir = basedir+'database/Swarm/ephemeris/'
    sampling_step = "PT1S"
    smoothWindowStr = '5s'
    strideval = int(smoothWindowStr[:-1])//int(sampling_step.replace('PT','')[:-1])
    
    collection = 'SW_OPER_EFI'+sat+'_LP_1B'
    parquetsuff = '_'+smoothWindowStr+'median_ephem_v1.parquet'           # nyeste
    
    # Get timestrs
    strftimefmt="%Y%m%d"
    timestrs = []
    iteryrs = range(yr0,yr1+1)
    if (2013 in iteryrs) or (2014 in iteryrs):
        start_time = datetime(2013, 12, 10)
        end_time = datetime(2014, 12, 31)
    
        timestrs.append("_{:s}-{:s}".format(
            start_time.strftime(strftimefmt),
            end_time.strftime(strftimefmt)))
    
    for yr in range(yr0,yr1+1):
        if yr in [2013,2014]:
            continue
    
        start_time = datetime(yr, 1, 1)
        end_time = datetime(yr, 12, 31)
    
        timestrs.append("_{:s}-{:s}".format(
            start_time.strftime(strftimefmt),
            end_time.strftime(strftimefmt)))
    
    fmtstr = "{:30s} : {:s}"
    print(fmtstr.format("outdir",outdir))
    print(fmtstr.format("sampling_step",sampling_step))
    print(fmtstr.format("smoothWindowStr",smoothWindowStr))
    print(fmtstr.format("Stride val",str(strideval)))
    print(fmtstr.format("years",",".join([str(yr) for yr in iteryrs])))
    print(fmtstr.format("satellite",sat))
    print(fmtstr.format("parquetinfo",sampling_step+timestrs[0]+parquetsuff))

    if dotiming:
        import time
        start_time = time.time()
    
    df = []

    for timestr in timestrs:
        parquetfile = collection+'_'+sampling_step+timestr+parquetsuff
    
        if not os.path.exists(outdir+parquetfile):
            print("No Existe: {:s}\nSkipping!".format(parquetfile))
            continue
    
        if dotiming:
            start_load_time = time.time()

        tmpdf = pd.read_parquet(outdir+parquetfile)

        if dotiming:
            stop_load_time = time.time()
            print("--- LOAD {:s}: {:.2f} seconds ---".format(parquetfile,stop_load_time-start_load_time))
    
        df.append(tmpdf)
    

    if dotiming:
        start_load_time = time.time()

    df = pd.concat(df)

    if dotiming:
        stop_load_time = time.time()
        print("--- CONCAT: {:.2f} seconds ---".format(stop_load_time-start_load_time))
        print("--- {:.2f} seconds total ---".format(time.time()-start_time))

    # Why on Earth did all QDLat values get divided by 1000?
    if dotiming:
        start_load_time = time.time()

    df.loc[:,'QDLat'] = df.loc[:,'QDLat']*1000

    df.loc[:,'mltbearing'] = mlt_bearing(df['QDLat'].values,
                                       df['MLT'].values,
                                       shift_to_90deg=False)
    df['mltbearingShift'] = df['mltbearing']
    more90ers = df['mltbearingShift'] > 90
    df.loc[more90ers,'mltbearingShift'] = df.loc[more90ers,'mltbearingShift'] -180
    df['mltbearingShift'] = np.abs(df['mltbearingShift'])

    if dotiming:
        stop_load_time = time.time()
        print("--- Fix QDLat, mltbearing: {:.2f} seconds ---".format(stop_load_time-start_load_time))

    if dotiming:
        print("--- {:.2f} seconds total ---".format(time.time()-start_time))

    return df


def load_2015_corr_db(date=None):
    """Get the Kalle 2015 correlation DB"""
    import pandas as pd
    # import numpy as np
    import os
    from functools import reduce
    datDir = '/SPENCEdata/Research/database/Swarm/201810-correlation/'
    fns = os.listdir(datDir)
    outFile = 'SwarmA-2015dBpar_Ne_correlation.pd'

    if date is not None:
        fil = date+'_corr.pd'
        if not os.path.isfile(datDir+'SwarmA/'+fil):
            print("Couldn't get " + fil + '!')
            return 0
        else:
            print("Loading " + fil + ' ...')
            dats = pd.read_pickle(datDir+'SwarmA/'+fil).dropna()
    else:
        if not os.path.isfile(datDir+outFile):
            print("Making " + outFile + ' ...')
            dats = reduce(lambda x, y: pd.concat((x, y)), (pd.read_pickle(
                datDir + 'SwarmA/' + fn).dropna() for fn in fns)).sort_index()
            print("Saving " + outFile + ' ...')
            dats.to_pickle(datDir+outFile)
        else:
            print("Reading " + outFile + ' ...')
            dats = pd.read_pickle(datDir+outFile)

    return dats


def get_Swarm_combo(dates,
                    sat='A',
                    get_dtypes=['MAG'],
                    dtype__add_Apex=[],
                    apex__geodetic2apexOpts=dict(min_time_resolution__sec=1,
                                                 max_N_months_twixt_apexRefTime_and_obs=3),
                    localSaveDir='/media/spencerh/data/Swarm/',
                    # check_for_existing_funcDict=dict(),
                    removeDownloadedZipAfterProcessing=False,
                    useKallesHardDrive=False,
                    only_download_zips=False,
                    only_list=False,
                    only_remove_zips=False,
                    opts_hurtigLast=dict()):

    if (localSaveDir == '/media/spencerh/data/Swarm/') and not useKallesHardDrive:
        localSaveDir = '/SPENCEdata/Research/database/Swarm/'

    if len(opts_hurtigLast.keys()) == 0:
        breakpoint()
        opts_hurtigLast = dict(FP__doCorrectTimestamps=False,
                               FP__doResample=False,
                               dont_touch_data=False,
                               dontInterp__justMag=False,
                               doDebug=False,
                               overwrite_existing=False,
                               use_existing=True,
                               removeCDF=True,
                               resampleString='500ms',
                               customSaveSuff='',
                               make_pickles=True)

    try:
        _ = dates[0]
        if isinstance(dates[0], str):
            dateStrs = dates
        else:
            dateStrs = [dater.strftime("%Y%m%d") for dater in dates]
    except:
        if isinstance(dates, str):
            dates = [dates]
            dateStrs = dates
        else:
            dates = [dates]
            dateStrs = [dates.strftime("%Y%m%d")]

    getFuncDict = dict(MAG=getMagFTP,
                       LP=getLPFTP,
                       FP=getFPFTP,
                       CT2HZ=getCT2HzFTP)

    # Make sure want dtypes are in getFuncDict
    if not isinstance(get_dtypes, list):
        if not isinstance(get_dtypes, str):
            print("Must provide a string or list of strings from among these: {:s}".format(
                ", ".join(list(getFuncDict.keys()))))
            return None

        get_dtypes = [get_dtypes]

    for dtyper in get_dtypes:

        if dtyper.upper() not in list(getFuncDict.keys()):

            print("{:s} not supported! Options are {:s}".format(
                dtyper.upper(),
                ", ".join(list(getFuncDict.keys()))))

            return None

    if only_list:
        print("Showing what I WOULD do")
    elif only_remove_zips:
        print("Only removing zips")
    elif only_download_zips:
        print("Only downloading zips")

    getFunc__onlylist = only_list or only_remove_zips
    outList = []
    deleteList = []
    for dtyper in get_dtypes:

        print("Getting {:s} ...".format(dtyper.upper()))

        dfList = []
        getFunc = getFuncDict[dtyper.upper()]

        gotFiles = getFunc(sat=sat,
                           dates=dateStrs,
                           only_list=getFunc__onlylist,
                           localSaveDir=localSaveDir,
                           append_dir=only_remove_zips or removeDownloadedZipAfterProcessing)

        if gotFiles is None:
            print("No files for {:s}!".format(",".join(dateStrs)))
            continue

        if only_list or only_download_zips:
            outList.append(gotFiles)
            continue
        elif (only_remove_zips or removeDownloadedZipAfterProcessing):
            # outList = outList + gotFiles
            deleteList = deleteList + gotFiles
            if only_remove_zips:
                continue
            # continue

        for ranDate in dates:

            df = sPH.hurtigLast(sat=sat,
                                doProcessMag=dtyper.upper() == 'MAG',
                                doProcessLP=dtyper.upper() == 'LP',
                                doProcessLowRes=False,
                                doProcessFP=dtyper.upper() == 'FP',
                                doProcessCT=dtyper.upper() == 'CT2HZ',
                                ranDate=ranDate,
                                useKallesHardDrive=useKallesHardDrive,
                                **opts_hurtigLast)

            if df is not None:
                dfList.append(df)

        if len(dfList) != 0:
            dfList = pd.concat(dfList)

            if (dtyper in dtype__add_Apex) and (dtyper != 'MAG'):

                print("Adding Apex coords to {:s}".format(dtyper))

                gdlat, gdalt_km = hCoord.geoclatR2geodlatheight(
                    dfList["Latitude"].values, dfList["Radius"].values/1000.)

                apexDict2 = hCoord.geodetic2apex(gdlat, dfList["Longitude"].values,
                                                 gdalt_km,
                                                 dfList.index.to_pydatetime(),
                                                 **apex__geodetic2apexOpts)

                dfList = dfList.assign(**apexDict2)

        else:
            dfList = pd.DataFrame()

        outList.append(dfList)

    if removeDownloadedZipAfterProcessing or only_remove_zips:

        # breakpoint()
            
        if len(deleteList) > 0:

            # if dtyper.upper() == 'MAG':
            #     lilsuff = 'MAGx_HR/Swarm_'+sat+'/'

            # elif dtyper.upper() == 'LP':
            #     lilsuff = '/EFI_LP/Swarm_'+sat+'/'
            # elif dtyper.upper() == 'FP':
            #     lilsuff = 'EFI_Faceplate_dens/Swarm_'+sat+'/'
            # elif dtyper.upper() == 'CT2HZ':
            #     lilsuff = '2Hz_TII_Cross-track/Swarm_'+sat+'/'

            # tmplocaldir = localSaveDir+lilsuff

            # for gotFile in gotFiles:
            #     if os.path.isfile(tmplocaldir+gotFile) and (gotFile[-4:].lower() == '.zip'):
            #         print("Removing {:s}".format(gotFile))
            #         os.remove(tmplocaldir+gotFile)

            for gotFile in deleteList:
                if isinstance(gotFile,list):
                    assert len(gotFile) == 1,"Not set up to handle multiple files this way! make gotFile single file!"
                    gotFil = gotFile[0]
                else:
                    gotFil = gotFile
                    # gotFile = gotFile[0]

                if os.path.isfile(gotFil) and (gotFil[-4:].lower() == '.zip'):
                    print("Removing {:s}".format(gotFil))
                    os.remove(gotFil)
                else:
                    print("Want to remove but don't find {:s}!".format(gotFil))

    return outList


def getMagFTP(sat='A',
              dates=None,
              # localSaveDir='/SPENCEdata/Research/Satellites/Swarm/rawFTP/MAG_HR/'):
              localSaveDir='/media/spencerh/data/Swarm/',
              only_list=False,
              check_for_existing_func=None,
              append_dir=False):
    # Swarm_A/SW_OPER_MAGA_

    localSaveDir += 'MAGx_HR/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"
    # outDir = '/SPENCEdata/Research/Satellites/Swarm/rawFTP/'
    # localSaveDir =

    subDir = '/Level1b/Entire_mission_data/MAGx_HR/Sat_'+sat+'/'

    # EXAMPLE: SW_OPER_MAGA_HR_1B_20131125T110251_20131125T235959_0505.CDF.ZIP
    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir,
                                only_list=only_list)

    if append_dir:
        gotFiles = [localSaveDir+gotFile for gotFile in gotFiles]

    return gotFiles


def getLPFTP(sat='A',
             dates=None,
             # localSaveDir='/SPENCEdata/Research/Satellites/Swarm/rawFTP/EFI_LP/'):
             localSaveDir='/media/spencerh/data/Swarm/',
             only_list=False,
             check_for_existing_func=None,
             append_dir=False):

    localSaveDir += '/EFI_LP/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = '/Level1b/Entire_mission_data/EFIx_LP/Sat_'+sat+'/'

    # EXAMPLE: SW_OPER_EFIA_LP_1B_20131202T101113_20131202T140109_0501.CDF.ZIP
    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir,
                                only_list=only_list)

    if append_dir and (gotFiles is not None):
        gotFiles = [localSaveDir+gotFile for gotFile in gotFiles]

    return gotFiles


def getFPFTP(sat='A',
             dates=None,
             localSaveDir='/media/spencerh/data/Swarm/',
             only_list=False,
             check_for_existing_func=None,
             append_dir=False):
    """
    Get a faceplate file
    """

    localSaveDir += 'EFI_Faceplate_dens/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = '/Advanced/Plasma_Data/16_Hz_Faceplate_plasma_density/Sat_'+sat+'/'

    # EXAMPLE: SW_EXTD_EFIA_LP_FP_20141002T092332_20141002T235958_0102.CDF
    # weirdSuff = '0102'
    # ftpFilePref = 'SW_EXTD_EFI'+sat+'_LP_FP_'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir,
                                only_list=only_list)

    if append_dir:
        gotFiles = [localSaveDir+gotFile for gotFile in gotFiles]

    return gotFiles


def getCT2HzFTP(sat='A',
                dates=None,
                localSaveDir='/media/spencerh/data/Swarm/',
                calversion='0201',
                only_list=False,
                check_for_existing_func=None,
                append_dir=False):
    """
    Get a cross-track 2-Hz file
    """

    VALID_CAL_VERSIONS = ['0201','0101']

    assert calversion in VALID_CAL_VERSIONS
    if calversion == '0201':
        subfolder = 'New_baseline' 
    elif calversion == '0101':
        subfolder = 'Old_baseline' 

    localSaveDir += '2Hz_TII_Cross-track/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = f'/Advanced/Plasma_Data/2Hz_TII_Cross-track_Dataset/{subfolder:s}/Sat_{sat:s}/'


    # EXAMPLE: SW_EXPT_EFIA_TIICT_20151101T155814_20151101T235004_0101.CDF.ZIP
    # ftpFilePref = 'SW_EXPT_EFI'+sat+'_TIICT_'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir,
                                only_list=only_list)

    if append_dir:
        if gotFiles is not None:
            gotFiles = [localSaveDir+gotFile for gotFile in gotFiles]

    return gotFiles


def _getFTP_dateGlob(dates, localSaveDir, subDir,
                     only_list=False):
    """
    Get a Swarm FTP file, genericizliaed
    """

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    ftp = ftplib.FTP(swarmFTPAddr)
    ftp.login()                 # Anonymous
    ftp.cwd(subDir)

    filz = ftp.nlst(subDir)
    ftpFiles = []

    if isinstance(dates, str):
        dates = [dates]
    elif isinstance(dates, list):

        for date in dates:
            if not isinstance(date, str):
                print("Must provide date strings 'YYYYMMDD'!")
                return None

    else:

        if only_list:
            ftp.close()
            return filz
        else:
            assert isinstance(
                dates, list), "Must provide list of date strings or a single date string (YYYYMMDD format)! (Or set kw only_list == True)"

    # Pick up all the files that match provided dates
    for date in dates:

        for f in filz:
            if fnmatch.fnmatch(f, '*'+date+'*'):
                # print(f)
                ftpFiles.append(f)
                # break

    if only_list:
        ftp.close()
        return ftpFiles

    # If no files found, exit
    if len(ftpFiles) == 0:
        print("Found no file! Exiting ...")
        ftp.close()
        return None

    # Junk the already-havers
    # ftpNotHavers = [ftpFile in ftpFiles if not os.path.isfile(localSaveDir + ftpFile)]
    ftpNotHavers = []
    for ftpFile in ftpFiles:
        fileExists = os.path.isfile(localSaveDir + ftpFile)

        if fileExists:
            if os.stat(localSaveDir+ftpFile).st_size == 0:
                print("{:s}:  File size is zero! Trying to get på nytt".
                      format(ftpFile))
                os.remove(localSaveDir+ftpFile)
        elif not fileExists:
            ftpNotHavers.append(ftpFile)


    if len(ftpNotHavers) == 0:
        print("Already have all {:d} files for {:d} date(s) provided! Exiting ...".format(
            len(ftpFiles), len(dates)))
        ftp.close()
        return ftpFiles

    print("Found {:d} files for the {:d} date(s) provided ({:d} are already downloaded)".format(
        len(ftpFiles), len(dates), len(ftpFiles)-len(ftpNotHavers)))

    # Get all matching files
    for ftpFile in ftpNotHavers:

        # Make sure we don't already have file
        if not os.path.isfile(localSaveDir + ftpFile):

            if not os.path.isdir(localSaveDir):
                os.mkdir(localSaveDir)

            with open(localSaveDir+ftpFile, "wb") as getFile:
                print("Trying to get " + ftpFile + ' ...')
                try:
                    ftp.retrbinary("RETR " + ftpFile, getFile.write)
                    print("Done!")
                except:
                    print("Couldn't get "+ftpFile+"!")

        else:
            print("Already have " + ftpFile + '!')

    ftp.close()

    return ftpFiles


def getMagFileDateRange(fName):
    """
    Skal begynne med 'SW_OPER_MAGA_HR_1B_'
    """

    yr0, mo0, day0 = int(fName[19:23]), int(fName[23:25]), int(fName[25:27])
    hr0, min0, sec0 = int(fName[28:30]), int(fName[30:32]), int(fName[32:34])
    yr1, mo1, day1 = int(fName[35:39]), int(fName[39:41]), int(fName[41:43])
    hr1, min1, sec1 = int(fName[44:46]), int(fName[46:48]), int(fName[48:50])
    return [datetime(yr0, mo0, day0, hr0, min0, sec0),
            datetime(yr1, mo1, day1, hr1, min1, sec1)]


def getLPFileDateRange(fName):
    """
    Skal begynne med 'SW_OPER_EFIB_LP_1B_'
    """

    yr0, mo0, day0 = int(fName[19:23]), int(fName[23:25]), int(fName[25:27])
    hr0, min0, sec0 = int(fName[28:30]), int(fName[30:32]), int(fName[32:34])
    yr1, mo1, day1 = int(fName[35:39]), int(fName[39:41]), int(fName[41:43])
    hr1, min1, sec1 = int(fName[44:46]), int(fName[46:48]), int(fName[48:50])
    return [datetime(yr0, mo0, day0, hr0, min0, sec0),
            datetime(yr1, mo1, day1, hr1, min1, sec1)]


def getFPFileDateRange(fName):
    """
    Skal begynne med 'SW_EXTD_EFIA_LP_FP_'
    """
    yr0, mo0, day0 = int(fName[19:23]), int(fName[23:25]), int(fName[25:27])
    hr0, min0, sec0 = int(fName[28:30]), int(fName[30:32]), int(fName[32:34])
    yr1, mo1, day1 = int(fName[35:39]), int(fName[39:41]), int(fName[41:43])
    hr1, min1, sec1 = int(fName[44:46]), int(fName[46:48]), int(fName[48:50])
    return [datetime(yr0, mo0, day0, hr0, min0, sec0),
            datetime(yr1, mo1, day1, hr1, min1, sec1)]


def getCT2HzFileDateRange(fName):
    """
    Skal begynne med 'SW_EXPT_EFIA_TIICT_'
    """
    yr0, mo0, day0 = int(fName[19:23]), int(fName[23:25]), int(fName[25:27])
    hr0, min0, sec0 = int(fName[28:30]), int(fName[30:32]), int(fName[32:34])
    yr1, mo1, day1 = int(fName[35:39]), int(fName[39:41]), int(fName[41:43])
    hr1, min1, sec1 = int(fName[44:46]), int(fName[46:48]), int(fName[48:50])
    return [datetime(yr0, mo0, day0, hr0, min0, sec0),
            datetime(yr1, mo1, day1, hr1, min1, sec1)]


def getMAGFilesDataFrame(**kws):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    df = _getFTP_FilesDataFrame(getType='MAG',
                                **kws)

    return df


def getLPFilesDataFrame(**kws):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    df = _getFTP_FilesDataFrame(getType='LP',
                                **kws)

    return df


def getFPFilesDataFrame(**kws):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    df = _getFTP_FilesDataFrame(getType='FP',
                                **kws)

    return df


def getCT2HzFilesDataFrame(**kws):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    df = _getFTP_FilesDataFrame(getType='CT2Hz',
                                **kws)

    return df


def _getFTP_FilesDataFrame(sat='A',
                           getType='CT2Hz',
                           t1=None,
                           t2=None,
                           add_OMNI_data=False,
                           OMNI_stats_list=['mean'],
                           OMNI_columns=None,
                           saveLocal=True,
                           localSaveDir='/SPENCEdata/Research/database/Swarm/buffered_lists/'):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    if OMNI_columns == None:
        OMNI_columns = ['Bz', 'By', 'Bx',
                        'AE', 'SymH', 'vsw', 'psw', 'borovsky',
                        'newell', 'F107', 'Kp', 'Dst']

    if getType == 'MAG':
        ftpFunc = getMagFTP
        descrip = 'B-field'
        fileDateParseFunc = getMagFileDateRange
    if getType == 'LP':
        ftpFunc = getLPFTP
        descrip = 'Langmuir Probe'
        fileDateParseFunc = getLPFileDateRange
    elif getType == 'FP':
        ftpFunc = getFPFTP
        descrip = 'faceplate'
        fileDateParseFunc = getFPFileDateRange
    elif getType == 'CT2Hz':
        ftpFunc = getCT2HzFTP
        descrip = '2Hz Cross-track'
        fileDateParseFunc = getCT2HzFileDateRange

    print(
        "Getting {:s} files dataframe for satellite {:s} ...".format(descrip,
                                                                     sat))

    omniStr = ''
    if add_OMNI_data:
        omniStr = '_w_OMNI'

    bufFile = 'Swarm'+sat+'_'+getType+'_files'+omniStr+'.pkl'

    try:
        these = ftpFunc(sat=sat,
                        only_list=True)
        starts = np.array([fileDateParseFunc(that)[0] for that in these])
        stops = np.array([fileDateParseFunc(that)[1] for that in these])

        lens = stops-starts
        lens_tSec = (stops-starts)/timedelta(seconds=1)

        df = pd.DataFrame(dict(start=starts, stop=stops,
                               dt=lens,
                               dt_sec=lens_tSec))

        if saveLocal:
            print("Saving {:s} list to {:s} ...".format(bufFile, localSaveDir))

            df.to_pickle(localSaveDir+bufFile)

    except:
        print("Couldn't get file from net!", end=' ')
        if os.path.isfile(localSaveDir + bufFile):
            print("Loading local file {:s} ...".format(bufFile))
            df = pd.read_pickle(localSaveDir+bufFile)
        else:
            print(
                "No local file! Sorry, no Swarm-{:s} {:s} files for you ...".format(sat, getType))

            return None

    havet1 = t1 is not None
    havet2 = t2 is not None

    if havet1:
        print("Filtering for dates >= {:s}".format(
            t1.strftime("%Y-%m-%d %H:%M:%S")))

    if havet2:
        print("Filtering for dates <= {:s}".format(
            t2.strftime("%Y-%m-%d %H:%M:%S")))

    if havet1 and havet2:
        df['match'] = (df['start'] >= t1) & (df['stop'] <= t2)
    elif havet1:
        df['match'] = df['start'] >= t1
    elif havet2:
        df['match'] = df['stop'] <= t2

    if add_OMNI_data:
        print("Adding OMNI data ...")
        these = hOMNI.omni_getter(df.iloc[0]['start'], df.iloc[-1]['stop'])

        these = these[OMNI_columns]

        addCols = list(these.columns)

        # addColsMean = [addCol + 'mean' for addCol in addCols]
        # addColsMedian = [addCol + 'median' for addCol in addCols]

        allAddCols = []
        for statType in OMNI_stats_list:
            allAddCols.append([addCol + statType for addCol in addCols])
            df[allAddCols[-1]
               ] = pd.DataFrame([[np.nan]*len(addCols)], index=df.index)

        for i, (star, stop) in enumerate(zip(df['start'].values, df['stop'].values)):
            # print(star,stop)
            thesens = (these.index >= star) & (these.index <= stop)
            zisMany = np.where(thesens)[0].size
            # print("{:3d}  {:s}--{:s} : Have {:4d}".format(i,
            #                                               pd.Timestamp(star).strftime("%Y-%m-%d %H:%M:%S"),
            #                                               pd.Timestamp(stop).strftime("%Y-%m-%d %H:%M:%S"),
            #                                               zisMany))
            if zisMany > 0:
                # df.loc[df.iloc[i].name, addCols] = these[thesens].mean().values
                for statType, statCols in zip(OMNI_stats_list, allAddCols):
                    df.loc[df.iloc[i].name, statCols] = these[thesens].agg(
                        statType).values

    return df
