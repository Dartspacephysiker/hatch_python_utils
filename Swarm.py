import ftplib
import os.path
import fnmatch


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


def getMagFTP(sat='A', date=None):

    swarmFTPAddr = "swarm-diss.eo.esa.int"
    # outDir = '/SPENCEdata/Research/Satellites/Swarm/rawFTP/'
    outMagDir = '/SPENCEdata/Research/Satellites/Swarm/rawFTP/MAG_HR/'

    magHRDir = '/Level1b/Entire_mission_data/MAGx_HR/Sat_'+sat+'/'
    magHRPref = 'SW_OPER_MAG'+sat+'_HR_1B_'
    magHRFile = magHRPref+date+'T000000_'+date+'T235959_0505.CDF.ZIP'

    if not os.path.isfile(outMagDir + magHRFile):
        ftp = ftplib.FTP(swarmFTPAddr)
        ftp.login()                 # Anonymous
        ftp.cwd(magHRDir)

        with open(outMagDir+magHRFile, "wb") as getFile:
            print("Trying to get " + magHRFile + ' ...')
            try:
                ftp.retrbinary("RETR " + magHRFile, getFile.write)
                print("Done!")
            except:
                print("Couldn't get mag file!")
    else:
        print("Already have " + magHRFile + '!')

    ftp.close()


def getLPFTP(sat='A', date=None):

    swarmFTPAddr = "swarm-diss.eo.esa.int"
    outLPDir = '/SPENCEdata/Research/Satellites/Swarm/rawFTP/EFI_LP/'

    EFILPdir = '/Level1b/Entire_mission_data/EFIx_LP/Sat_'+sat+'/'
    EFILPPref = 'SW_OPER_EFI'+sat+'_LP_1B_'
    # EFILPFile = EFILPPref + date + 'T101113_'+date+'T140109_0501.CDF.ZIP'

    ftp = ftplib.FTP(swarmFTPAddr)
    ftp.login()                 # Anonymous
    ftp.cwd(EFILPdir)

    filz = ftp.nlst(EFILPdir)
    for f in filz:
        if fnmatch.fnmatch(f, '*'+date+'*'):
            print(f)
            EFILPFile = f
            break

    if not os.path.isfile(outLPDir + EFILPFile):

        with open(outLPDir+EFILPFile, "wb") as getFile:
            print("Trying to get " + EFILPFile + ' ...')
            try:
                ftp.retrbinary("RETR " + EFILPFile, getFile.write)
                print("Done!")
            except:
                print("Couldn't get EFI file!")
    else:
        print("Already have " + EFILPFile + '!')

    ftp.close()
