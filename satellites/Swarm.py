import ftplib
import os.path
import fnmatch
from datetime import datetime


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


def getMagFTP(sat='A',
              dates=None,
              # localSaveDir='/SPENCEdata/Research/Satellites/Swarm/rawFTP/MAG_HR/'):
              localSaveDir='/media/spencerh/data/Swarm/'):
    # Swarm_A/SW_OPER_MAGA_

    localSaveDir += 'MAGx_HR/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"
    # outDir = '/SPENCEdata/Research/Satellites/Swarm/rawFTP/'
    # localSaveDir =

    subDir = '/Level1b/Entire_mission_data/MAGx_HR/Sat_'+sat+'/'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir)

    return gotFiles

    # ftpFilePref = 'SW_OPER_MAG'+sat+'_HR_1B_'
    # ftpFile = ftpFilePref+date+'T000000_'+date+'T235959_0505.CDF.ZIP'

    # if not os.path.isfile(localSaveDir + ftpFile):
    #     ftp = ftplib.FTP(swarmFTPAddr)
    #     ftp.login()                 # Anonymous
    #     ftp.cwd(subDir)

    #     with open(localSaveDir+ftpFile, "wb") as getFile:
    #         print("Trying to get " + ftpFile + ' ...')
    #         try:
    #             ftp.retrbinary("RETR " + ftpFile, getFile.write)
    #             print("Done!")
    #         except:
    #             print("Couldn't get mag file!")

    #     ftp.close()

    # else:
    #     print("Already have " + ftpFile + '!')


def getLPFTP(sat='A',
             dates=None,
             # localSaveDir='/SPENCEdata/Research/Satellites/Swarm/rawFTP/EFI_LP/'):
             localSaveDir='/media/spencerh/data/Swarm/'):

    localSaveDir += '/EFI_LP/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = '/Level1b/Entire_mission_data/EFIx_LP/Sat_'+sat+'/'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir)

    return gotFiles

    # ftpFilePref = 'SW_OPER_EFI'+sat+'_LP_1B_'
    # ftpFile = ftpFilePref + date + 'T101113_'+date+'T140109_0501.CDF.ZIP'

    # ftpFile = None

    # ftp = ftplib.FTP(swarmFTPAddr)
    # ftp.login()                 # Anonymous
    # ftp.cwd(subDir)

    # filz = ftp.nlst(subDir)
    # for f in filz:
    #     if fnmatch.fnmatch(f, '*'+date+'*'):
    #         print(f)
    #         ftpFile = f
    #         break

    # if ftpFile is None:
    #     print("Found no file! Exiting ...")
    #     return

    # if not os.path.isfile(localSaveDir + ftpFile):

    #     with open(localSaveDir+ftpFile, "wb") as getFile:
    #         print("Trying to get " + ftpFile + ' ...')
    #         try:
    #             ftp.retrbinary("RETR " + ftpFile, getFile.write)
    #             print("Done!")
    #         except:
    #             print("Couldn't get EFI file!")

    #     ftp.close()

    # else:
    #     print("Already have " + ftpFile + '!')


def getFPFTP(sat='A',
             dates=None,
             localSaveDir='/media/spencerh/data/Swarm/'):
    """
    Get a faceplate file
    """

    localSaveDir += 'EFI_Faceplate_dens/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = '/Advanced/Plasma_Data/16_Hz_Faceplate_plasma_density/Sat_'+sat+'/'

    # EXAMPLE: SW_EXTD_EFIA_LP_FP_20141002T092332_20141002T235958_0102.CDF
    # weirdSuff = '0102'
    # ftpFilePref = 'SW_EXTD_EFI'+sat+'_LP_FP_'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir)

    return gotFiles

    # ftpFile = None

    # ftp = ftplib.FTP(swarmFTPAddr)
    # ftp.login()                 # Anonymous
    # ftp.cwd(subDir)

    # filz = ftp.nlst(subDir)
    # for f in filz:
    #     if fnmatch.fnmatch(f, '*'+date+'*'):
    #         print(f)
    #         ftpFile = f
    #         break

    # if ftpFile is None:
    #     print("Found no file! Exiting ...")
    #     return

    # if not os.path.isfile(localSaveDir + ftpFile):

    #     with open(localSaveDir+ftpFile, "wb") as getFile:
    #         print("Trying to get " + ftpFile + ' ...')
    #         try:
    #             ftp.retrbinary("RETR " + ftpFile, getFile.write)
    #             print("Done!")
    #         except:
    #             print("Couldn't get FP file!")

    #     ftp.close()

    # else:
    #     print("Already have " + ftpFile + '!')


def getCT2HzFTP(sat='A',
                dates=None,
                localSaveDir='/media/spencerh/data/Swarm/'):
    """
    Get a cross-track 2-Hz file
    """

    localSaveDir += '2Hz_TII_Cross-track/Swarm_'+sat+'/'

    swarmFTPAddr = "swarm-diss.eo.esa.int"

    subDir = '/Advanced/Plasma_Data/2Hz_TII_Cross-track_Dataset/Sat_{:s}/'.format(
        sat)

    # EXAMPLE: SW_EXPT_EFIA_TIICT_20151101T155814_20151101T235004_0101.CDF.ZIP
    # ftpFilePref = 'SW_EXPT_EFI'+sat+'_TIICT_'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir)

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

    if only_list:
        return filz

    if isinstance(dates, str):
        dates = [dates]
    else:
        assert isinstance(
            dates, list), "Must provide list of date strings or a single date string (YYYYMMDD format)!"

    # Pick up all the files that match provided dates
    for date in dates:

        for f in filz:
            if fnmatch.fnmatch(f, '*'+date+'*'):
                # print(f)
                ftpFiles.append(f)
                break

    # If no files found, exit
    if len(ftpFiles) == 0:
        print("Found no file! Exiting ...")
        ftp.close()
        return None

    # Junk the already-havers
    # ftpNotHavers = [ftpFile in ftpFiles if not os.path.isfile(localSaveDir + ftpFile)]
    ftpNotHavers = []
    for ftpFile in ftpFiles:
        if not os.path.isfile(localSaveDir + ftpFile):
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
