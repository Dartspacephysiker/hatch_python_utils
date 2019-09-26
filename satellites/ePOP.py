import swarmProcHelper as sPH
import conjunctions as swarmC
import ftplib
import os.path
import fnmatch
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
from hatch_python_utils import omni_utils as hOMNI

import sys
needdir = '/SPENCEdata/Research/sandbox_and_journals/journals/Swarm/'
if needdir not in sys.path:
    sys.path.insert(0, needdir)


instrSuffDict = dict(CER=['zip', 'png', 'TEC.png'],
                     GAP=['14O', 'lv0', 'lv1', 'png'],
                     IRM=['mp4', 'lv0', 'lv1', 'png'],
                     MGF=['lv0', 'lv1', 'lv2', 'lv2Cal',
                          'lv3', 'GEI.lv3', 'SC.lv3',
                          'png'],
                     RRI=['summary.txt', 'h5.zip', 'lv0', 'png'],
                     SEI=['mp4', 'lv0', 'lv1', 'png'],
                     ephemeris=['.txt'])


def dateparse(yyyymmddstr):
    return yyyymmddstr[0:4], yyyymmddstr[4:6], yyyymmddstr[6:]


def getePOPFTP(dates=None,
               localSaveDir='/media/spencerh/data/e-POP/',
               instr='CER',
               get_ephemeris=True,
               instrsuffs=None,
               only_list=True):
    """
    Get e-POP files
    """

    instrDirs = ['CER', 'GAP', 'IRM', 'MGF', 'RRI', 'SEI', 'ephemeris']

    gotFiles = _getFTP_dateGlob_ePOP(dates, localSaveDir,
                                     instr=instr,
                                     get_ephemeris=get_ephemeris,
                                     instrsuffs=instrsuffs,
                                     only_list=only_list)

    return gotFiles


def print_instrSuffDict():
    global instrSuffDict

    for key in instrSuffDict.keys():
        print("{:10s} : {:s}".format(
            key, ", ".join(instrSuffDict[key])))


def _getFTP_dateGlob_ePOP(dates, localSaveDir,
                          instr='CER',
                          instrsuffs=None,
                          get_ephemeris=True,
                          only_list=True):
    """
    Get a Swarm-Echo FTP file, genericizliaed
    """

    global instrSuffDict

    if isinstance(instr, str):
        instr = [instr]

    instrDirs = ['CER', 'GAP', 'IRM', 'MGF', 'RRI', 'SEI', 'ephemeris']

    for ins in instr:
        if ins not in instrDirs:
            print("Must set 'instr' to one of {:s}!".format(
                ",".join(instrDirs)))
            return None

    if instrsuffs is not None:

        if not isinstance(instrsuffs, dict):
            print(
                "instrsuffs should be a dict with the types of file exts you want to pull up")
            print_instrSuffDict()

        for ins in instr:
            if ins in instrsuffs.keys():
                if isinstance(instrsuffs[ins], str):
                    instrsuffs[ins] = [instrsuffs[ins]]

                print("{:10s} : {:s}".format(ins, ", ".join(instrsuffs[ins])))

        # return None

    print(
        "Getting e-POP files for {:s} instruments ...".format(", ".join(instr)))

    instrFileListDict = {key: {} for key in instr}

    if 'ephemeris' in instr:
        instr.remove('ephemeris')
        get_ephemeris = True
        instrFileListDict['ephemeris'] = {}

    epopFTPAddr = 'swarm-diss.eo.esa.int'
    epopFTPDir = '/#CASSIOPE_e-POP'
    ftp = ftplib.FTP(epopFTPAddr)
    ftp.login()                 # Anonymous

    if isinstance(dates, str):
        dates = [dates]
    elif isinstance(dates, list):

        for date in dates:
            if not isinstance(date, str):
                print("Must provide date strings 'YYYYMMDD'!")
                return None

    # Pick up all the files that match provided dates
    goodDirs = []
    goodDates = []
    for date in dates:
        yr, mo, day = dateparse(date)

        print(yr, mo, day)

        tryDir = epopFTPDir+'/'+'/'.join([yr, mo, day])
        try:
            yunk = ftp.cwd(tryDir)
            goodDirs.append(tryDir)
            goodDates.append(date)
        except:
            print("Nothing for {:s}".format(date))

    for goodDate, goodDir in zip(goodDates, goodDirs):

        yunk = ftp.cwd(goodDir)

        tmpDirCont = ftp.nlst()

        if get_ephemeris:
            tmpephemlist = []
            for item in tmpDirCont:
                if 'ephem' in item:
                    tmpephemlist.append(item)

            instrFileListDict['ephemeris'][goodDate] = [
                goodDir+'/'] + tmpephemlist

        for ins in instr:

            if ins in tmpDirCont:

                yunk = ftp.cwd(goodDir+'/'+ins)

                if instrsuffs is not None:

                    if ins in instrsuffs.keys():
                        tmplist = ftp.nlst()
                        keepers = []
                        for instrsuff in instrsuffs[ins]:
                            for tmper in tmplist:
                                if instrsuff in tmper:
                                    keepers.append(tmper)

                        instrFileListDict[ins][goodDate] = keepers

                    else:

                        instrFileListDict[ins][goodDate] = ftp.nlst()

                else:
                    instrFileListDict[ins][goodDate] = ftp.nlst()

            if len(instrFileListDict[ins][goodDate]) > 0:
                instrFileListDict[ins][goodDate] = [
                    goodDir+'/'+ins+'/'] + instrFileListDict[ins][goodDate]

    if only_list:
        ftp.close()
        return instrFileListDict

    # If no files found, exit
    # if len(ftpFiles) == 0:
    #     print("Found no file! Exiting ...")
    #     ftp.close()
    #     return None

    # Now loop over everything
    ftpNotHavers = []
    localNotHavers = []
    totFileCount = 0
    for ins in instrFileListDict.keys():

        tmpDir = localSaveDir+ins+'/'

        for date in instrFileListDict[ins].keys():

            for ftpFile in instrFileListDict[ins][date]:

                if ftpFile[-1] == '/':
                    continue

                totFileCount += 1

                if not os.path.isfile(tmpDir + ftpFile):
                    print("Not have {:s}!".format(tmpDir+ftpFile))
                    localNotHavers.append(tmpDir+ftpFile)
                    ftpNotHavers.append(
                        instrFileListDict[ins][date][0]+ftpFile)

    if len(ftpNotHavers) == 0:
        print("Already have all {:d} files for {:d} date(s) provided! Exiting ...".format(
            totFileCount, len(dates)))
        ftp.close()
        return instrFileListDict

    print("Found {:d} files for the {:d} date(s) provided ({:d} are already downloaded)".format(
        totFileCount, len(dates), totFileCount-len(ftpNotHavers)))

    # return instrFileListDict
    # return localNotHavers, ftpNotHavers

    if totFileCount > 10:
        cont = False
        while not cont:
            svar = input(
                "More than 10 files requested! You sure you want to proceed (y/n/showme)?")

            if svar.lower() == 'n':
                print("OK, exiting ...")
                ftp.close()
                return instrFileListDict
            elif svar.lower() == 'showme':
                print("Værsågod:")
                for i, (localFile, ftpFile) in enumerate(zip(localNotHavers, ftpNotHavers)):
                    print(i, localFile, ftpFile)

            elif svar.lower() == 'y':
                cont = True

    # Get all matching files
    for localFile, ftpFile in zip(localNotHavers, ftpNotHavers):

        # Make sure we don't already have file
        if not os.path.isfile(localFile):

            # if not os.path.isdir(localSaveDir):
            #     os.mkdir(localSaveDir)

            with open(localFile, "wb") as getFile:
                print("Trying to get " + ftpFile + ' ...')
                try:
                    ftp.retrbinary("RETR " + ftpFile, getFile.write)
                    print("Done!")
                except:
                    print("Couldn't get "+ftpFile+"!")

        else:
            print("Already have " + ftpFile + '!')

    ftp.close()

    return instrFileListDict


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


def getMAGFilesDataFrame(**kws):
    """
    t1: start of desired interval (datetime-like object)
    t2: end   of desired interval (datetime-like object)
    """

    df = _getFTP_FilesDataFrame(getType='MAG',
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
