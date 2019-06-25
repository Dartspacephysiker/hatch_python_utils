
from spacepy import pycdf

import numpy as np
import pandas as pd
import datetime

import fnmatch
import ftplib
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import os
import os.path
import scipy.io as sio

import matplotlib as mpl
import tkinter
mplBkgrnd = 'TkAgg'
mpl.use(mplBkgrnd)


FASTFTPAddr = "cdaweb.gsfc.nasa.gov"
tight_rect = [0, 0.03, 1, 0.95]


def toTStamps(dts):
    return np.array([tid.timestamp() for tid in dts])


def getTeamsCDF(date,
                doDownload=True,
                quiet=False,
                attemptFTPDownloadForMissing=True):
    """
    Get path of local TEAMS CDF file
    date: 'YYYYMMDD' string (e.g., date.strftime("%Y%m%d")

    KWARGS
    ------
    doDownload : attempt to download from FTP if doesn't exist locally

    """

    localSaveDir = '/SPENCEdata/Research/database/FAST/TEAMS_example/'

    year = date[0:4]

    localDir = localSaveDir + year + '/'
    filz = os.listdir(localDir)

    matchFiles = []
    for f in filz:
        if fnmatch.fnmatch(f, '*'+date+'*'):
            matchFiles.append(localDir+f)

    if len(matchFiles) == 0:

        notHaveString = "{:s} {:s} doesn't exist locally!".format("TEAMS CDF for date",
                                                                  date)

        if attemptFTPDownloadForMissing:

            if not quiet:
                print(notHaveString + ' Attempting FTP download ...')

            matchFiles = getTeamsFTP(dates=date, localSaveDir=localSaveDir)

        else:
            print(notHaveString + ' Exiting ...')
            matchFiles = []

    # if not quiet:
    #     print(matchFiles)

    if len(matchFiles) == 1:
        matchFiles = matchFiles[0]

    return matchFiles


def getTeamsFTP(dates=None,
                localSaveDir='/SPENCEdata/Research/database/FAST/TEAMS_example/'):

    # fName = 'fa_k0_tms_19980923_v01.cdf'

    # localSaveDir += 'Swarm_'+sat+'/'

    subDir = '/pub/data/fast/teams/k0/'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir)

    return gotFiles


def _getFTP_dateGlob(dates, localSaveDir, subDir):
    """
    Get a FAST FTP file, genericizliaed
    """

    ftp = ftplib.FTP(FASTFTPAddr)
    ftp.login()                 # Anonymous
    junk = ftp.cwd(subDir)
    if isinstance(dates, str):
        dates = [dates]
    else:
        assert isinstance(
            dates, list), "Must provide list of date strings or a single date string (YYYYMMDD format)!"
    # How many years do we have?
    years = set([tring[0:4] for tring in dates])
    ftpFiles = []
    for year in years:
        tmpDir = subDir+year+'/'
        filz = ftp.nlst(tmpDir)

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
        # ftp.close()
        # return None

    # Junk the already-havers
    ftpNotHavers = []
    ftpLocals = []
    for ftpFile in ftpFiles:

        localFile = localSaveDir+ftpFile.lstrip(subDir)
        ftpLocals.append(localFile)

        chkDir = os.path.dirname(localSaveDir+ftpFile.lstrip(subDir))
        if not os.path.exists(chkDir):
            print("Making " + chkDir + ' ...')
            os.mkdir(chkDir)

        if not os.path.isfile(localFile):
            ftpNotHavers.append(ftpFile)
        else:
            junk = ftp.sendcmd("TYPE i")
            if ftp.size(ftpFile) != os.stat(localFile).st_size:
                print("Differing file sizes: ", ftp.size(ftpFile),
                      os.stat(localFile).st_size, " so get anyway ...")
                # print("Get anyway!")
                ftpNotHavers.append(ftpFile)

            junk = ftp.sendcmd("TYPE A")

    if len(ftpNotHavers) == 0:
        print("Already have all {:d} files for {:d} date(s) provided! Exiting ...".format(
            len(ftpFiles), len(dates)))
        ftp.close()
        return ftpLocals

    print("Found {:d} files for the {:d} date(s) provided ({:d} are already downloaded)".format(
        len(ftpFiles), len(dates), len(ftpFiles)-len(ftpNotHavers)))

    # Get all matching files
    for ftpFile in ftpNotHavers:

        localFile = localSaveDir+ftpFile.lstrip(subDir)

        # print(localFile)
        # print(ftpFile)

        # Make sure we don't already have file
        # if not os.path.isfile(localFile):

        with open(localFile, "wb") as getFile:
            print("Trying to get " + ftpFile + ' ... ', end='')
            try:
                ftp.retrbinary("RETR " + ftpFile, getFile.write)
                print("Done!")
            except:
                print("Couldn't get "+ftpFile+"!")

        # else:
        #     print("Already have " + ftpFile + '!')

    ftp.close()

    return ftpLocals


class TEAMS:

    # Epoch: CDF_EPOCH [2381]
    # H+: CDF_REAL4 [2381, 48]
    # H+_en: CDF_REAL4 [2381, 48]
    # H+_high: CDF_REAL4 [2381, 16]
    # H+_high_pa: CDF_REAL4 [2381, 16]
    # H+_low: CDF_REAL4 [2381, 16]
    # H+_low_pa: CDF_REAL4 [2381, 16]
    # He+: CDF_REAL4 [2381, 48]
    # He+_en: CDF_REAL4 [2381, 48]
    # He+_high: CDF_REAL4 [2381, 16]
    # He+_high_pa: CDF_REAL4 [2381, 16]
    # He+_low: CDF_REAL4 [2381, 16]
    # He+_low_pa: CDF_REAL4 [2381, 16]
    # O+: CDF_REAL4 [2381, 48]
    # O+_en: CDF_REAL4 [2381, 48]
    # O+_high: CDF_REAL4 [2381, 16]
    # O+_high_pa: CDF_REAL4 [2381, 16]
    # O+_low: CDF_REAL4 [2381, 16]
    # O+_low_pa: CDF_REAL4 [2381, 16]
    # Time_PB5: CDF_INT4 [2381, 3]
    # alt: CDF_REAL4 [2381]
    # cartesian: CDF_CHAR*1 [3] NRV
    # fa_spin_dec: CDF_REAL4 [2381]
    # fa_spin_ra: CDF_REAL4 [2381]
    # flat: CDF_REAL4 [2381]
    # flng: CDF_REAL4 [2381]
    # format_time: CDF_CHAR*2 [3] NRV
    # ilat: CDF_REAL4 [2381]
    # label_r: CDF_CHAR*5 [3] NRV
    # label_time: CDF_CHAR*27 [3] NRV
    # label_v: CDF_CHAR*6 [3] NRV
    # mlt: CDF_REAL4 [2381]
    # orbit: CDF_INT4 [2381]
    # post_gap_flag: CDF_INT4 [2381]
    # quality_flag: CDF_INT4 [2381]
    # r: CDF_REAL4 [2381, 3]
    # unit_time: CDF_CHAR*4 [3] NRV
    # unix_time: CDF_REAL8 [2381]
    # v: CDF_REAL4 [2381, 3]

    def __init__(self, dateOrCDFFile,
                 attemptFTPDownloadForMissing=True):

        isDateStr = False
        isOrbStr = False
        if isinstance(dateOrCDFFile, str):
            if len(dateOrCDFFile) == 8:
                isDateStr = True
            elif len(dateOrCDFFile) == 4:
                isOrbStr = True

        if isDateStr:
            cdfFile = getTeamsCDF(dateOrCDFFile,
                                  attemptFTPDownloadForMissing=attemptFTPDownloadForMissing)
        elif isOrbStr:
            orbCDFDir = '/SPENCEdata/Research/database/FAST/TEAMS_example/by_orbit__HOMEMADE/'
            orbCDFFile = 'fa_tms_'+str(dateOrCDFFile)+'.cdf'
            cdfFile = orbCDFDir+orbCDFFile
            print("Homemade file: "+orbCDFFile)
        else:
            cdfFile = dateOrCDFFile

        try:
            cdf = pycdf.CDF(cdfFile)
        except:
            print("Couldn't load CDF file for " + dateOrCDFFile + '!')

            return None

        # wantVarlist = ['tStamps', 'orbit', 'mlt', 'ilat', 'alt'
        #                'O', 'O_high', 'O_high_pa', 'O_low', 'O_low_pa', 'O_en', 'OF',
        #                'H', 'H_high', 'H_high_pa', 'H_low', 'H_low_pa', 'H_en', 'HF',
        #                'He', 'He_high', 'He_high_pa', 'He_low', 'He_low_pa', 'He_en', 'HeF',
        #                'quality_flag']

        keys = cdf.keys()

        # endeligVarList = []
        # for var in wantVarlist:
        #     if var in keys:
        #         endeligVarList.append(var)
        # self.varlist = endeligVarList

        wantOneDlist = ['tStamps', 'orbit',
                        'mlt', 'ilat', 'alt', 'quality_flag']

        self.oneDlist = []
        for var in wantOneDlist:
            if var in keys:
                self.oneDlist.append(var)

        wantDistVarList = ['O', 'O_high', 'O_high_pa', 'O_low', 'O_low_pa', 'O_en', 'OF',
                           'H', 'H_high', 'H_high_pa', 'H_low', 'H_low_pa', 'H_en', 'HF',
                           'He', 'He_high', 'He_high_pa', 'He_low', 'He_low_pa', 'He_en', 'HeF']
        self.distvarlist = []
        for var in wantDistVarList:
            if var in keys:
                self.distvarlist.append(var)

        self.varlist = self.oneDlist+self.distvarlist

        self.tStamps = pd.DatetimeIndex(
            list(map(pd.Timestamp, cdf["Epoch"][:])))

        if 'orbit' in keys:
            self.orbit = cdf["orbit"][:]
        elif isOrbStr:
            self.orbit = np.array([np.int64(dateOrCDFFile)]*self.tStamps.size)

        if 'mlt' in keys:
            self.mlt = cdf["mlt"][:]
        if 'ilat' in keys:
            self.ilat = cdf["ilat"][:]
        if 'alt' in keys:
            self.alt = cdf["alt"][:]

        self.O = cdf["O+"][:, :]          # Time, energy(?)
        self.O_high = cdf["O+_high"][:, :]          # Time, energy(?)
        self.O_high_pa = cdf["O+_high_pa"][:, :]          # Time, energy(?)
        self.O_low = cdf["O+_low"][:, :]          # Time, energy(?)
        self.O_low_pa = cdf["O+_low_pa"][:, :]          # Time, energy(?)
        self.O_en = cdf["O+_en"][:, :]          # Time, energy(?)

        self.H = cdf["H+"][:, :]          # Time, energy(?)
        self.H_high = cdf["H+_high"][:, :]          # Time, energy(?)
        self.H_high_pa = cdf["H+_high_pa"][:, :]          # Time, energy(?)
        self.H_low = cdf["H+_low"][:, :]          # Time, energy(?)
        self.H_low_pa = cdf["H+_low_pa"][:, :]          # Time, energy(?)
        self.H_en = cdf["H+_en"][:, :]          # Time, energy(?)

        self.He = cdf["He+"][:, :]          # Time, energy(?)
        self.He_high = cdf["He+_high"][:, :]          # Time, energy(?)
        self.He_high_pa = cdf["He+_high_pa"][:, :]          # Time, energy(?)
        self.He_low = cdf["He+_low"][:, :]          # Time, energy(?)
        self.He_low_pa = cdf["He+_low_pa"][:, :]          # Time, energy(?)
        self.He_en = cdf["He+_en"][:, :]          # Time, energy(?)

        # From ~/software/sdt/idl/convert_tms_units2.pro
        # To get to flux (not eflux) units, divide each channel by the relevant energy
        self.OF = self.O/self.O_en          # Time, energy(?)
        self.HF = self.H/self.H_en          # Time, energy(?)
        self.HeF = self.He/self.He_en          # Time, energy(?)

        self.init_consts()

        if 'quality_flag' in keys:
            self.quality_flag = cdf["quality_flag"][:]

        supercheck = (self.tStamps.size == self.O.shape[0]) \
            and (self.tStamps.size == self.O_high.shape[0]) \
            and (self.tStamps.size == self.O_high_pa.shape[0]) \
            and (self.tStamps.size == self.O_low.shape[0]) \
            and (self.tStamps.size == self.O_low_pa.shape[0]) \
            and (self.tStamps.size == self.O_en.shape[0]) \
            and (self.tStamps.size == self.H.shape[0]) \
            and (self.tStamps.size == self.H_high.shape[0]) \
            and (self.tStamps.size == self.H_high_pa.shape[0]) \
            and (self.tStamps.size == self.H_low.shape[0]) \
            and (self.tStamps.size == self.H_low_pa.shape[0]) \
            and (self.tStamps.size == self.H_en.shape[0]) \
            and (self.tStamps.size == self.He.shape[0]) \
            and (self.tStamps.size == self.He_high.shape[0]) \
            and (self.tStamps.size == self.He_high_pa.shape[0]) \
            and (self.tStamps.size == self.He_low.shape[0]) \
            and (self.tStamps.size == self.He_low_pa.shape[0]) \
            and (self.tStamps.size == self.He_en.shape[0])  # \
        # and (self.tStamps.size == self.orbit.size) \
        # and (self.tStamps.size == self.ilat.size) \
        # and (self.tStamps.size == self.mlt.size) \
        # and (self.tStamps.size == self.quality_flag.size)

        if not supercheck:
            print("these arrays are not all the same size!")

        self.O_denergy = (np.diff(self.O_en, axis=1)*-1)
        self.O_denergy = np.vstack([self.O_denergy[:, 0], self.O_denergy.T]).T

        self.H_denergy = (np.diff(self.H_en, axis=1)*-1)
        self.H_denergy = np.vstack([self.H_denergy[:, 0], self.H_denergy.T]).T

        self.He_denergy = (np.diff(self.He_en, axis=1)*-1)
        self.He_denergy = np.vstack(
            [self.He_denergy[:, 0], self.He_denergy.T]).T

    def init_consts(self):
        self.Omass = 0.165695  # Mass eV/(km/sec)^2
        self.OgeomFactor = 0.0015

        self.Hmass = 0.0104389  # Mass eV/(km/sec)^2
        self.HgeomFactor = 0.0015

        self.Hemass = 0.0414521  # Mass eV/(km/sec)^2

        self.OConst = np.sqrt(self.Omass * 1.6e-22/(2.*1.6e-12))
        self.HConst = np.sqrt(self.Hmass * 1.6e-22/(2.*1.6e-12))
        self.HeConst = np.sqrt(self.Hemass * 1.6e-22/(2.*1.6e-12))

        # Odensity = OConst*total(denergy*(energy^(-1.5))*sumdata)

    def set_keep_filter(self,
                        minMLT=None, maxMLT=None,
                        minILAT=None, maxILAT=None,
                        hemi='both',
                        minAlt=None, maxAlt=None,
                        t0=None, t1=None,
                        verbose=True):

        origTall = self.mlt.size
        runningTall = origTall

        junkString = "Junked {:8d} inds due to {:14s} ({:8d} matches)"
        freshKeeps = np.ones(self.mlt.size, dtype=np.bool)

        if hasattr(self, 'quality_flag'):
            newKeeps = (self.quality_flag == 0)
            if verbose:
                print(junkString.format(
                    runningTall - keeps[keeps].size, "quality_flag", origTall-newKeeps[newKeeps].size))
                runningTall = keeps[keeps].size
        else:
            newKeeps = np.ones(self.mlt.size, dtype=np.bool)

        keeps = newKeeps

        ####################
        # Times
        if (t0 is not None) or (t1 is not None):
            newKeeps = freshKeeps.copy()

            if t0 is not None:
                newKeeps = newKeeps & (self.tStamps >= t0)

                if t1 is not None:
                    newKeeps = newKeeps & (self.tStamps <= t1)

            keeps = keeps & newKeeps

            if verbose:
                print(junkString.format(
                    runningTall - keeps[keeps].size, "time restric", origTall-newKeeps[newKeeps].size))
                runningTall = keeps[keeps].size

        ####################
        # MLTs
        if (minMLT is not None) or (maxMLT is not None):
            newKeeps = freshKeeps.copy()

            if minMLT is not None:
                newKeeps = newKeeps & (self.mlt >= minMLT)

                if maxMLT is not None:
                    newKeeps = newKeeps & (self.mlt <= maxMLT)

            keeps = keeps & newKeeps

            if verbose:
                print(junkString.format(
                    runningTall - keeps[keeps].size, "MLT lims", origTall-newKeeps[newKeeps].size))
                runningTall = keeps[keeps].size

        ####################
        # alts
        if (minAlt is not None) or (maxAlt is not None):
            newKeeps = freshKeeps.copy()

            if minAlt is not None:
                newKeeps = newKeeps & (self.alt >= minAlt)

                if maxAlt is not None:
                    newKeeps = newKeeps & (self.alt <= maxAlt)

            keeps = keeps & newKeeps

            if verbose:
                print(junkString.format(
                    runningTall - keeps[keeps].size, "alt lims", origTall-newKeeps[newKeeps].size))
                runningTall = keeps[keeps].size

        ####################
        # ilats
        if (minILAT is not None) or (maxILAT is not None):
            newKeeps = freshKeeps.copy()

            if minILAT is not None:
                if hemi.lower()[0] == 'b':
                    newKeeps = newKeeps & (np.abs(self.ilat) >= minILAT)
                elif hemi.lower()[0] == 'n':
                    newKeeps = newKeeps & (self.ilat >= minILAT)
                elif hemi.lower()[0] == 's':
                    newKeeps = newKeeps & (self.ilat <= minILAT)

            if maxILAT is not None:
                if hemi.lower()[0] == 'b':
                    newKeeps = newKeeps & (np.abs(self.ilat) <= maxILAT)
                elif hemi.lower()[0] == 'n':
                    newKeeps = newKeeps & (self.ilat <= maxILAT)
                elif hemi.lower()[0] == 's':
                    newKeeps = newKeeps & (self.ilat >= maxILAT)

            keeps = keeps & newKeeps

            if verbose:
                print(junkString.format(
                    runningTall - keeps[keeps].size, "ilat lims", origTall-newKeeps[newKeeps].size))
                runningTall = keeps[keeps].size

        if verbose:
            print(
                "Junked {:8d}/{:8d} total inds".format(origTall - runningTall, origTall))

        self.keeps = keeps

    def addattrs(self, names, values,
                 verbose=False):

        if isinstance(names, str):
            setattr(self, names, values)
        else:

            if len(names) != len(values):
                print("Imbalance ...")
                return

            for name, value in zip(names, values):

                if verbose:
                    print("Adding " + name + ' ...')
                setattr(self, name, value)

    def interpthisquantandadd(self, quant, tStamps, quantname,
                              verbose=True,
                              **interpOpts):

        interpedQuant = np.interp(toTStamps(self.tStamps.to_pydatetime()),
                                  toTStamps(tStamps), quant, **interpOpts)

        # print(quantname)
        # print("quantsize  : {:d}".format(quant.size))
        # print("tstampssize: {:d}".format(tStamps.size))
        # print("interpqsize: {:d}".format(interpedQuant.size))
        # print("selftStsize: {:d}".format(self.tStamps.size))

        self.addattrs(quantname,
                      interpedQuant,
                      verbose=verbose)


class OutflowAlgorithm:

    def __init__(self, orbit):

        outflowAlgDir = '/SPENCEdata/Research/database/FAST/twoTypes_ion_identification/'
        fila = 'orbit_' + str(orbit) + \
            '__outflow_algorithm_and_beam_algorithm.sav'

        bro = sio.readsav(outflowAlgDir+fila, python_dict=True)

        self.tStamps = np.array([datetime.datetime.utcfromtimestamp(
            x) for x in bro['ionmomstruct'][0]['time'].byteswap().newbyteorder()])

        self.eranges = bro['ionmomstruct'][0]['erange'].byteswap(
        ).newbyteorder()
        self.types = bro['ionmomstruct']['type'][0].byteswap().newbyteorder()

        self.mlt = bro['ephemstruct'][0]['mlt'].byteswap().newbyteorder()
        self.ilat = bro['ephemstruct'][0]['ilat'].byteswap().newbyteorder()
        self.alt = bro['ephemstruct'][0]['alt'].byteswap().newbyteorder()
        self.lat = bro['ephemstruct'][0]['lat'].byteswap().newbyteorder()
        self.lng = bro['ephemstruct'][0]['lng'].byteswap().newbyteorder()

        supercheck = (self.tStamps.size == self.mlt.size) \
            and (self.tStamps.size == self.ilat.size) \
            and (self.tStamps.size == self.alt.size) \
            and (self.tStamps.size == self.lat.size) \
            and (self.tStamps.size == self.lng.size) \
            and (self.tStamps.size == self.types.size) \
            and (self.tStamps.size == self.eranges.shape[0])

        if not supercheck:
            print("these arrays are not all the same size!")


def plotTEAMSspecs(teams, outflowAlg, useWants, integwants, algWants,
                   doDiffFlux=False,
                   addHe=False,
                   t0=None,
                   t1=None,
                   doLogY=True,
                   title=None,
                   # y_limsO=[1,10000],
                   # y_limsH=[1,10000],
                   # cmap = 'inferno', #https://matplotlib.org/examples/color/colormaps_reference.html
                   cmap='magma',
                   showInterpedEranges=True,
                   quiet=True,
                   laptop=False):

    algTimes = mdates.date2num(outflowAlg.tStamps)

    if doDiffFlux:
        cbTitle = 'Log(#/cm$^2$-s-sr-eV)'
        z_lims = [0, 6]
        OQuant = teams.OFSub.T
        HQuant = teams.HFSub.T
        HeQuant = teams.HeFSub.T
    else:
        cbTitle = 'Log(eV/cm$^2$-s-sr-eV)'
        z_lims = [2, 7]
        OQuant = teams.OSub.T
        HQuant = teams.HSub.T
        HeQuant = teams.HeSub.T

    # START (needifusingimshow)

    # if doLogY:
    #     y_limsO = np.log10(y_limsO)
    #     y_limsH = np.log10(y_limsH)

    if t0 is not None:
        startTime = t0
    else:
        startTime = teams.tStamps[0]

    if t0 is not None:
        stopTime = t0
    else:
        stopTime = teams.tStamps[-1]

    tlims = mdates.date2num([startTime, stopTime])

    erangeInds = (algTimes >= tlims[0]) & (algTimes <= tlims[1])

    # STOP (needifusingimshow)

    if laptop:
        figSize = (12, 6)
    else:
        figSize = (16, 8)

    specStamps = mdates.date2num(teams.tStamps[useWants])
    specStampsV, specEnV = np.meshgrid(specStamps, teams.O_enSub[0, :])

    if addHe:
        fig, (ax, bx, cx) = plt.subplots(3, 1, sharex=True, figsize=figSize)
    else:
        fig, (ax, bx) = plt.subplots(2, 1, sharex=True, figsize=figSize)

    if title is not None:
        junk = fig.suptitle(title)

    OplusEnSpec = ax.pcolormesh(specStampsV, np.log10(specEnV), np.log10(OQuant), cmap=cmap,
                                vmin=z_lims[0],
                                vmax=z_lims[1])

    HplusEnSpec = bx.pcolormesh(specStampsV, np.log10(specEnV), np.log10(HQuant), cmap=cmap,
                                vmin=z_lims[0],
                                vmax=z_lims[1])

    junk = ax.set_ylabel("$O^+$ Log Energy (eV)")
    ax.xaxis_date()
    cba = fig.colorbar(OplusEnSpec, ax=ax)
    junk = cba.set_label(cbTitle)

    junk = bx.set_ylabel("$H^+$ Log Energy (eV)")
    bx.xaxis_date()
    cbb = fig.colorbar(HplusEnSpec, ax=bx)
    junk = cbb.set_label(cbTitle)

    if addHe:
        HeplusEnSpec = cx.pcolormesh(specStampsV, np.log10(specEnV), np.log10(HeQuant), cmap=cmap,
                                     vmin=z_lims[0],
                                     vmax=z_lims[1])

        junk = cx.set_ylabel("$He^+$ Log Energy (eV)")
        cx.xaxis_date()
        cbc = fig.colorbar(HeplusEnSpec, ax=cx)
        junk = cbc.set_label(cbTitle)

    if showInterpedEranges:
        specStampsERanges = mdates.date2num(teams.tStamps[integwants])
        junk = ax.plot(specStampsERanges, np.log10(
            teams.eranges[integwants, 0]), color='blue')
        junk = bx.plot(specStampsERanges, np.log10(
            teams.eranges[integwants, 0]), color='blue')

        junk = ax.plot(specStampsERanges, np.log10(
            teams.eranges[integwants, 1]), color='orange')
        junk = bx.plot(specStampsERanges, np.log10(
            teams.eranges[integwants, 1]), color='orange')

        if addHe:
            junk = cx.plot(specStampsERanges, np.log10(
                teams.eranges[integwants, 0]), color='blue')
            junk = cx.plot(specStampsERanges, np.log10(
                teams.eranges[integwants, 1]), color='orange')

    else:
        junk = ax.plot(algTimes[erangeInds], np.log10(
            outflowAlg.eranges[erangeInds, 0]), color='blue')
        junk = bx.plot(algTimes[erangeInds], np.log10(
            outflowAlg.eranges[erangeInds, 0]), color='blue')

        junk = ax.plot(algTimes[erangeInds], np.log10(
            outflowAlg.eranges[erangeInds, 1]), color='red')
        junk = bx.plot(algTimes[erangeInds], np.log10(
            outflowAlg.eranges[erangeInds, 1]), color='red')

        if addHe:
            junk = cx.plot(algTimes[erangeInds], np.log10(
                outflowAlg.eranges[erangeInds, 0]), color='blue')
            junk = cx.plot(algTimes[erangeInds], np.log10(
                outflowAlg.eranges[erangeInds, 1]), color='red')

    plt.tight_layout(rect=tight_rect)

    if addHe:
        return fig, (ax, bx, cx)
    else:
        return fig, (ax, bx)


def plotTEAMShists(teams, integwants,
                   logFluxBins=None,
                   logRatioBins=None,
                   doDiffFlux=False,
                   title=None):

    fig, (ax, bx) = plt.subplots(2, 1)

    if title is not None:
        junk = fig.suptitle(title)

    alpha = 0.5
    Oopts = dict(color='blue', alpha=alpha, label="$O^+$")
    Hopts = dict(color='red', alpha=alpha, label="$H^+$")
    Heopts = dict(color='orange', alpha=alpha, label="$He^+$")

    ratioAlpha = 0.8
    OHopts = dict(color=Oopts['color'], label="$O^+ / H^+$", alpha=ratioAlpha)
    HeHopts = dict(color=Heopts['color'],
                   label="$He^+ / H^+$", alpha=ratioAlpha)

    Ohistdat = np.log10(teams.OJ[np.isfinite(np.log10(teams.OJ))])
    Hhistdat = np.log10(teams.HJ[np.isfinite(np.log10(teams.HJ))])
    Hehistdat = np.log10(teams.HeJ[np.isfinite(np.log10(teams.HeJ))])

    junk = ax.hist(Ohistdat, bins=logFluxBins, **Oopts)
    junk = ax.hist(Hhistdat, bins=logFluxBins, **Hopts)
    junk = ax.hist(Hehistdat, bins=logFluxBins, **Heopts)

    aleg = ax.legend()

    if doDiffFlux:
        xTitle = 'Log(#/cm$^2$-s)'
    else:
        xTitle = 'Log(eV/cm$^2$-s)'

    junk = ax.set_xlabel(xTitle)

    junk = bx.hist(np.log10(teams.OHratio[np.isfinite(
        teams.OHratio)]), bins=logRatioBins, **OHopts)
    junk = bx.hist(np.log10(teams.HeHratio[np.isfinite(
        teams.HeHratio)]), bins=logRatioBins, **HeHopts)

    junk = bx.axvline(np.log10(teams.OHratio_median),
                      label="O/H Median", linestyle='--', color='black')
    junk = bx.axvline(np.log10(teams.HeHratio_median),
                      label="He/H Median", linestyle='-.', color='black')

    junk = bx.set_xlabel("Log(Species Ratio)")

    junk = bx.legend()

    plt.tight_layout(rect=tight_rect)

    return fig, (ax, bx)
