
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


def getESAL2CDF(date,
                instr='ees',
                doDownload=True,
                quiet=False,
                attemptFTPDownloadForMissing=True):
    """
    Get path of local ESA L2 CDF file
    date: 'YYYYMMDD' string (e.g., date.strftime("%Y%m%d")

    KWARGS
    ------
    doDownload : attempt to download from FTP if doesn't exist locally

    """

    localSaveDir = '/media/spencerh/data/FAST/ftp/ESA/'+instr+'/'

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

            matchFiles = getESAFTP(
                dates=date, localSaveDir=localSaveDir, instr=instr)

        else:
            print(notHaveString + ' Exiting ...')
            matchFiles = []

    # if not quiet:
    #     print(matchFiles)

    if len(matchFiles) == 1:
        matchFiles = matchFiles[0]

    return matchFiles


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


def getESAFTP(dates=None,
              instr='ees',
              localSaveDir='/media/spencerh/data/FAST/ftp/ESA/'):

    assert instr in ['ees', 'eeb', 'ies', 'ieb'], "Must use one of zese!"

    # fName = 'fa_esa_l2_ees_19961101005800_00774_v01.cdf'

    # localSaveDir += 'Swarm_'+sat+'/'

    # https://spdf.gsfc.nasa.gov/pub/data/fast/esa/l2/ees/1996/11/
    subDir = '/pub/data/fast/esa/l2/'+instr+'/'

    gotFiles = _getFTP_dateGlob(dates, localSaveDir, subDir,
                                descendToMonth=True)

    return gotFiles


def _getFTP_dateGlob(dates, localSaveDir, subDir,
                     descendToMonth=False):
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

        if descendToMonth:

            months = set([tring[4:6] for tring in dates if tring[0:4] == year])

            for month in months:

                tmpMonthDir = subDir+year+'/'+month+'/'
                filz = ftp.nlst(tmpDir)

                for date in dates:

                    for f in filz:
                        if fnmatch.fnmatch(f, '*'+date+'*'):
                            # print(f)
                            ftpFiles.append(f)
                            break

        else:

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


# This one is for TEAMS cdfs
class TEAMSdata:

    class TEAMSstruct:
        def __init__(self, bigdata, index):

            if index > bigdata.nStruct:
                print("Can't use index={:d} (only {:d} structs possible)! Returning ...".format(
                    index, bigdata.nStruct))
                return None

            self.valid = bigdata.valid[index]
            self.data_name = bigdata.data_name[index]

            # Time stuff
            self.time = bigdata.time[index]
            self.end_time = bigdata.end_time[index]
            self.integ_t = bigdata.integ_t[index]
            self.delta_time = bigdata.delta_time[index]

            # Units and mass stuff
            self.units_name = bigdata.units_name[index].decode("utf-8")
            self.units_procedure = bigdata.units_procedure[index].decode(
                "utf-8")
            self.project_name = bigdata.project_name[index].decode("utf-8")
            self.mass = bigdata.mass[index]

            # Data stuff
            self.nenergy = bigdata.nenergy[index]
            self.nbins = bigdata.nbins[index]
            self.orbit = bigdata.orbit[index]
            self.geomfactor = bigdata.geomfactor[index]

            self.data = bigdata.data[:, :, index]
            self.energy = bigdata.energy[:, :, index]
            self.theta = bigdata.theta[:, :, index]
            self.phi = bigdata.phi[:, :, index]
            self.denergy = bigdata.denergy[:, :, index]
            self.dtheta = bigdata.dtheta[:, index]
            self.dphi = bigdata.dphi[:, index]
            self.domega = bigdata.domega[:, index]

            self.geom = bigdata.geom[:, :, index]
            self.eff = bigdata.eff[:, :, index]
            self.spin_fract = bigdata.spin_fract[:, :, index]

            self.pt_limits = bigdata.pt_limits[:, index]

            self.alt = bigdata.alt[index]
            self.mlt = bigdata.mlt[index]
            self.ilat = bigdata.ilat[index]
            self.foot_lat = bigdata.foot_lat[index]
            self.foot_lng = bigdata.foot_lng[index]
            self.fa_spin_ra = bigdata.fa_spin_ra[index]
            self.fa_spin_dec = bigdata.fa_spin_dec[index]

            # Housekeeping
            self.header_bytes = bigdata.header_bytes[index]
            self.eff_version = bigdata.eff_version[index]

        def __repr__(self):
            attrs = [a for a in dir(self) if not a.startswith('__')]

            def stringOutput(a):
                if hasattr(getattr(self, a), 'shape'):
                    if getattr(self, a).shape:
                        return "{:15s}: {:10s} {:10s}".format(a,
                                                              str(getattr(
                                                                  self, a).shape),
                                                              str(getattr(self, a).dtype))
                    else:
                        return "{:15s}: {:10s} {:10s}".format(a,
                                                              "{:.2f}".format(
                                                                  getattr(self, a)),
                                                              str(getattr(self, a).dtype))

                elif isinstance(getattr(self, a), str):
                    return "{:15s}: {:10s} {:10s}".format(a,
                                                          getattr(self, a),
                                                          str(type(getattr(self, a))))
                elif isinstance(getattr(self, a), datetime.datetime):
                    return "{:15s}: {:10s} {:10s}".format(a,
                                                          str(getattr(self, a)),
                                                          str(type(getattr(self, a))))
                else:
                    return "{:15s}: {:10s} {:10s}".format(a,
                                                          'scalar',
                                                          str(type(getattr(self, a))))

            # fullStrs = ["{:15s}: {:10s} {:10s}".format(a,
            #                                            str(getattr(self,a).shape),
            #                                            str(getattr(self,a).dtype)) for a in attrs]
            fullStrs = [stringOutput(a) for a in attrs]

            # fullStrs = ["{:s}: {:s}".format(getattr(self,a),str(getattr(self,a).shape)) for a in attrs]
            return "\n".join(fullStrs)

    def __init__(self, filename):

        self.fileLoaded = False

        # filename = '/SPENCEdata/software/sdt/batch_jobs/TEAMS/orb1421_teams.sav'
        if not os.path.isfile(filename):
            print("Doesn't exist: "+filename)
            return

        this = sio.readsav(filename)

        self.fileLoaded = True

        navn = list(this.keys())[0]
        self.nStruct = len(this[navn])

        bigger = np.zeros(
            (*this[navn][0]['data'].T.shape, self.nStruct), dtype=np.float32)
        bigger2 = np.zeros(
            (*this[navn][0]['dtheta'].shape, self.nStruct), dtype=np.float32)

        self.valid = this[navn]['valid'].byteswap().newbyteorder()
        self.data_name = this[navn]['data_name'].byteswap().newbyteorder()

        # Time stuff
        self.time = np.array([datetime.datetime.utcfromtimestamp(
            x) for x in this[navn].time.byteswap().newbyteorder()])
        self.end_time = [datetime.datetime.utcfromtimestamp(
            x) for x in this[navn].end_time.byteswap().newbyteorder()]
        # self.end_time = this[navn]['end_time'].byteswap().newbyteorder()
        self.integ_t = this[navn]['integ_t'].byteswap().newbyteorder()
        self.delta_time = this[navn]['delta_time'].byteswap(
        ).newbyteorder()

        # Units and mass stuff
        self.units_name = this[navn]['units_name'].byteswap(
        ).newbyteorder()
        self.units_procedure = this[navn]['units_procedure'].byteswap(
        ).newbyteorder()
        # self.project_name = this[navn]['project_name'].byteswap().newbyteorder()
        self.project_name = this[navn]['project_name']
        # eV/c^2 with c in (km/s)^2 I think
        self.mass = this[navn]['mass'].byteswap().newbyteorder()

        # Data stuff
        self.nenergy = this[navn]['nenergy'].byteswap().newbyteorder()
        self.nbins = this[navn]['nbins'].byteswap().newbyteorder()
        self.orbit = this[navn]['orbit'].byteswap().newbyteorder()
        self.geomfactor = this[navn]['geomfactor'].byteswap(
        ).newbyteorder()

        # Ephemeris stuff
        self.alt = this[navn]['alt'].byteswap().newbyteorder()
        self.mlt = this[navn]['mlt'].byteswap().newbyteorder()
        self.ilat = this[navn]['ilat'].byteswap().newbyteorder()
        self.foot_lat = this[navn]['foot_lat'].byteswap().newbyteorder()
        self.foot_lng = this[navn]['foot_lng'].byteswap().newbyteorder()
        self.fa_spin_ra = this[navn]['fa_spin_ra'].byteswap(
        ).newbyteorder()
        self.fa_spin_dec = this[navn]['fa_spin_dec'].byteswap(
        ).newbyteorder()

        self.fa_pos = np.vstack(
            this[navn]['fa_pos'].byteswap().newbyteorder()).T
        self.fa_vel = np.vstack(
            this[navn]['fa_vel'].byteswap().newbyteorder()).T
        self.b_model = np.vstack(
            this[navn]['b_model'].byteswap().newbyteorder()).T
        self.b_foot = np.vstack(
            this[navn]['b_foot'].byteswap().newbyteorder()).T

        # Housekeeping
        self.header_bytes = this[navn]['header_bytes'].byteswap(
        ).newbyteorder()
        self.eff_version = this[navn]['eff_version'].byteswap(
        ).newbyteorder()

        # Have to just initialize these
        self.data = bigger.copy()
        self.energy = bigger.copy()
        self.denergy = bigger.copy()
        self.theta = bigger.copy()
        self.phi = bigger.copy()

        # biggers
        self.geom = bigger.copy()
        self.eff = bigger.copy()
        self.spin_fract = bigger.copy()

        # biggers2
        self.domega = bigger2.copy()
        self.dtheta = bigger2.copy()
        self.dphi = bigger2.copy()

        self.pt_limits = np.zeros(
            (*this[navn][0]['pt_limits'].shape, self.nStruct))

        for i in range(self.nStruct):

            self.data[:, :, i] = this[navn][i]['data'].byteswap(
            ).newbyteorder().T
            self.energy[:, :, i] = this[navn][i]['energy'].byteswap(
            ).newbyteorder().T
            self.denergy[:, :, i] = this[navn][i]['denergy'].byteswap(
            ).newbyteorder().T
            self.theta[:, :, i] = this[navn][i]['theta'].byteswap(
            ).newbyteorder().T
            self.phi[:, :, i] = this[navn][i]['phi'].byteswap().newbyteorder().T
            self.geom[:, :, i] = this[navn][i]['geom'].byteswap(
            ).newbyteorder().T
            self.eff[:, :, i] = this[navn][i]['eff'].byteswap().newbyteorder().T
            self.spin_fract[:, :, i] = this[navn][i]['spin_fract'].byteswap(
            ).newbyteorder().T

            self.dtheta[:, i] = this[navn][i]['dtheta'].byteswap(
            ).newbyteorder()
            self.dphi[:, i] = this[navn][i]['dphi'].byteswap().newbyteorder()
            self.domega[:, i] = this[navn][i]['domega'].byteswap(
            ).newbyteorder()

            self.pt_limits[:, i] = this[navn][0]['pt_limits'].byteswap(
            ).newbyteorder()

            # print(this[navn][i]['energy'].byteswap().newbyteorder()[0,:])
            # print(self.energy[:,:,i])

    def __call__(self, index):

        return self.TEAMSstruct(self, index)

    def apply_inds(self, keepInds, verbose=False):

        attrs = [a for a in dir(self) if not a.startswith('__')]

        for a in attrs:
            tmpAttr = getattr(self, a)
            if verbose:
                print(a, type(tmpAttr), isinstance(tmpAttr, np.ndarray))

            if isinstance(tmpAttr, np.ndarray):
                shape = np.array(tmpAttr.shape)
                # print(shape.size)
                if shape.size == 1:
                    if tmpAttr.size == self.nStruct:
                        if verbose:
                            print("Resizing " + a + "...")
                        setattr(self, a, tmpAttr[keepInds])
                elif shape.size == 2:
                    if verbose:
                        print("Resizing " + a + ": ",
                              tmpAttr[:, keepInds].shape)
                    setattr(self, a, tmpAttr[:, keepInds])
                elif shape.size == 3:
                    if verbose:
                        print("Resizing " + a + ": ",
                              tmpAttr[:, :, keepInds].shape)
                    setattr(self, a, tmpAttr[:, :, keepInds])

        self.nStruct = np.sum(keepInds)


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

        junkedStr = "Junked {:8d} inds due to {:14s} ({:8d} matches)"
        freshKeeps = np.ones(self.mlt.size, dtype=np.bool)

        if hasattr(self, 'quality_flag'):
            newKeeps = (self.quality_flag == 0)
            if verbose:
                print(junkedStr.format(
                    runningTall - newKeeps[newKeeps].size, "quality_flag", origTall-newKeeps[newKeeps].size))
                runningTall = newKeeps[newKeeps].size
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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


class EES:

    # angle_labl_64: CDF_CHAR*12 [64] NRV
    # bins: CDF_UINT1 [2, 64, 96] NRV
    # bins_ind: CDF_INT2 [1272]
    # bkg: CDF_FLOAT [1272]
    # bkg_arr: CDF_FLOAT [10, 64, 96] NRV
    # charge: CDF_INT2 [] NRV
    # compno_64: CDF_INT2 [64] NRV
    # compno_96: CDF_INT2 [96] NRV
    # data: CDF_UINT1 [1272, 64, 96]
    # data_level: CDF_CHAR*7 [] NRV
    # data_name: CDF_CHAR*11 [] NRV
    # data_quality: CDF_UINT1 [1272]
    # dead: CDF_FLOAT [] NRV
    # denergy: CDF_FLOAT [3, 64, 96] NRV
    # denergy_full: CDF_FLOAT [1272, 64, 96]
    # domega: CDF_FLOAT [1272, 64, 96]
    # dtheta: CDF_FLOAT [3, 64, 96] NRV
    # eff: CDF_FLOAT [3, 64, 96] NRV
    # eflux: CDF_FLOAT [1272, 64, 96]
    # eflux_byA_atE: CDF_FLOAT [0]
    # eflux_byE_atA: CDF_FLOAT [0]
    # eflux_byenergy_labl: CDF_CHAR*33 [96] NRV
    # eflux_bypitch_labl: CDF_CHAR*29 [64] NRV
    # eflux_movie: CDF_FLOAT [0]
    # energy: CDF_FLOAT [3, 64, 96] NRV
    # energy_full: CDF_FLOAT [1272, 64, 96]
    # energy_labl_96: CDF_CHAR*13 [96] NRV
    # energy_median: CDF_FLOAT [0]
    # epoch: CDF_EPOCH [1272]
    # geom_factor: CDF_FLOAT [1272]
    # gf: CDF_FLOAT [10, 64, 96] NRV
    # gf_ind: CDF_INT2 [1272]
    # header_bytes: CDF_UINT1 [1272, 44]
    # mass: CDF_FLOAT [] NRV
    # mode_ind: CDF_UINT1 [1272]
    # nbins: CDF_UINT1 [1272]
    # nenergy: CDF_UINT1 [1272]
    # num_dists: CDF_INT4 [] NRV
    # orbit_end: CDF_INT4 [] NRV
    # orbit_number: CDF_INT4 [2]
    # orbit_number_epoch: CDF_EPOCH [2]
    # orbit_number_time: CDF_DOUBLE [2]
    # orbit_start: CDF_INT4 [] NRV
    # pitch_angle: CDF_FLOAT [1272, 64, 96]
    # pitch_angle_median: CDF_FLOAT [0]
    # project_name: CDF_CHAR*4 [] NRV
    # sc_pot: CDF_FLOAT [1272]
    # theta: CDF_FLOAT [3, 64, 96] NRV
    # theta_max: CDF_FLOAT [1272]
    # theta_min: CDF_FLOAT [1272]
    # theta_shift: CDF_FLOAT [1272]
    # time_delta: CDF_DOUBLE [1272]
    # time_end: CDF_DOUBLE [1272]
    # time_integ: CDF_DOUBLE [1272]
    # time_start: CDF_DOUBLE [1272]
    # time_unix: CDF_DOUBLE [1272]
    # units_name: CDF_CHAR*10 [] NRV
    # units_procedure: CDF_CHAR*20 [] NRV
    # valid: CDF_INT2 [1272]

    def __init__(self, dateOrCDFFile,
                 attemptFTPDownloadForMissing=True):
        """
        Infos in
        spedas/idl/general/missions/fast/fa_esa/cdf_load/fa_load_esa_l1.pro
        """
        isDateStr = False
        isOrbStr = False
        if isinstance(dateOrCDFFile, str):
            if len(dateOrCDFFile) == 8:
                isDateStr = True
            elif len(dateOrCDFFile) == 4:
                isOrbStr = True

        if isDateStr:
            cdfFile = getESAL2CDF(dateOrCDFFile,
                                  attemptFTPDownloadForMissing=attemptFTPDownloadForMissing, instr='ees')
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

        keys = cdf.keys()

        # endeligVarList = []
        # for var in wantVarlist:
        #     if var in keys:
        #         endeligVarList.append(var)
        # self.varlist = endeligVarList

        wantOneDlist = ['bkg', 'compno_64', 'compno_96', 'time_start', 'time_end', 'sc_pot', 'epoch',
                        'orbit_number', 'data_quality', 'mode_ind', 'valid', 'geom_factor',
                        'bins_ind', 'gf_ind', 'nbins', 'nenergy',
                        'theta_shift', 'theta_min', 'theta_max']

        self.oneDlist = []
        for var in wantOneDlist:
            if var in keys:
                self.oneDlist.append(var)

        wantDistVarList = ['bkg_arr', 'bins', 'data', 'energy', 'theta', 'denergy', 'dtheta', 'domega', 'eff', 'gf',
                           'pitch_angle']
        self.distvarlist = []
        for var in wantDistVarList:
            if var in keys:
                self.distvarlist.append(var)

        self.varlist = self.oneDlist+self.distvarlist

        self.tStamps = pd.DatetimeIndex(
            list(map(pd.Timestamp, cdf["epoch"][:])))

        if 'orbit_number' in keys:
            self.orbit = cdf["orbit_number"][:]
        elif isOrbStr:
            self.orbit = np.array([np.int64(dateOrCDFFile)]*self.tStamps.size)

        for var in self.oneDlist:
            if var in ['time_start', 'time_end', 'time_unix']:
                setattr(self, var,
                        np.array([datetime.datetime.utcfromtimestamp(x) for x in cdf[var][:]]))
            else:
                setattr(self, var, cdf[var][:])

        for var in self.distvarlist:
            setattr(self, var, cdf[var][:, :, :])

        # breakpoint()
        # if 'mlt' in keys:
        #     self.mlt = cdf["mlt"][:]
        # if 'ilat' in keys:
        #     self.ilat = cdf["ilat"][:]
        # if 'alt' in keys:
        #     self.alt = cdf["alt"][:]

        # self.O = cdf["O+"][:, :]          # Time, energy(?)
        # self.O_high = cdf["O+_high"][:, :]          # Time, energy(?)
        # self.O_high_pa = cdf["O+_high_pa"][:, :]          # Time, energy(?)
        # self.O_low = cdf["O+_low"][:, :]          # Time, energy(?)
        # self.O_low_pa = cdf["O+_low_pa"][:, :]          # Time, energy(?)
        # self.O_en = cdf["O+_en"][:, :]          # Time, energy(?)

        # self.H = cdf["H+"][:, :]          # Time, energy(?)
        # self.H_high = cdf["H+_high"][:, :]          # Time, energy(?)
        # self.H_high_pa = cdf["H+_high_pa"][:, :]          # Time, energy(?)
        # self.H_low = cdf["H+_low"][:, :]          # Time, energy(?)
        # self.H_low_pa = cdf["H+_low_pa"][:, :]          # Time, energy(?)
        # self.H_en = cdf["H+_en"][:, :]          # Time, energy(?)

        # self.He = cdf["He+"][:, :]          # Time, energy(?)
        # self.He_high = cdf["He+_high"][:, :]          # Time, energy(?)
        # self.He_high_pa = cdf["He+_high_pa"][:, :]          # Time, energy(?)
        # self.He_low = cdf["He+_low"][:, :]          # Time, energy(?)
        # self.He_low_pa = cdf["He+_low_pa"][:, :]          # Time, energy(?)
        # self.He_en = cdf["He+_en"][:, :]          # Time, energy(?)

        # # From ~/software/sdt/idl/convert_tms_units2.pro
        # # To get to flux (not eflux) units, divide each channel by the relevant energy
        # self.OF = self.O/self.O_en          # Time, energy(?)
        # self.HF = self.H/self.H_en          # Time, energy(?)
        # self.HeF = self.He/self.He_en          # Time, energy(?)

        # self.init_consts()

        if 'data_quality' in keys:
            self.data_quality = cdf["data_quality"][:]

        # supercheck = (self.tStamps.size == self.O.shape[0]) \
        #     and (self.tStamps.size == self.O_high.shape[0]) \
        #     and (self.tStamps.size == self.O_high_pa.shape[0]) \
        #     and (self.tStamps.size == self.O_low.shape[0]) \
        #     and (self.tStamps.size == self.O_low_pa.shape[0]) \
        #     and (self.tStamps.size == self.O_en.shape[0]) \
        #     and (self.tStamps.size == self.H.shape[0]) \
        #     and (self.tStamps.size == self.H_high.shape[0]) \
        #     and (self.tStamps.size == self.H_high_pa.shape[0]) \
        #     and (self.tStamps.size == self.H_low.shape[0]) \
        #     and (self.tStamps.size == self.H_low_pa.shape[0]) \
        #     and (self.tStamps.size == self.H_en.shape[0]) \
        #     and (self.tStamps.size == self.He.shape[0]) \
        #     and (self.tStamps.size == self.He_high.shape[0]) \
        #     and (self.tStamps.size == self.He_high_pa.shape[0]) \
        #     and (self.tStamps.size == self.He_low.shape[0]) \
        #     and (self.tStamps.size == self.He_low_pa.shape[0]) \
        #     and (self.tStamps.size == self.He_en.shape[0])  # \
        # and (self.tStamps.size == self.orbit.size) \
        # and (self.tStamps.size == self.ilat.size) \
        # and (self.tStamps.size == self.mlt.size) \
        # and (self.tStamps.size == self.quality_flag.size)

        # if not supercheck:
        #     print("these arrays are not all the same size!")

        # self.O_denergy = (np.diff(self.O_en, axis=1)*-1)
        # self.O_denergy = np.vstack([self.O_denergy[:, 0], self.O_denergy.T]).T

        # self.H_denergy = (np.diff(self.H_en, axis=1)*-1)
        # self.H_denergy = np.vstack([self.H_denergy[:, 0], self.H_denergy.T]).T

        # self.He_denergy = (np.diff(self.He_en, axis=1)*-1)
        # self.He_denergy = np.vstack(
        #     [self.He_denergy[:, 0], self.He_denergy.T]).T

    # def init_consts(self):
        # self.Omass = 0.165695  # Mass eV/(km/sec)^2
        # self.OgeomFactor = 0.0015

        # self.Hmass = 0.0104389  # Mass eV/(km/sec)^2
        # self.HgeomFactor = 0.0015

        # self.Hemass = 0.0414521  # Mass eV/(km/sec)^2

        # self.OConst = np.sqrt(self.Omass * 1.6e-22/(2.*1.6e-12))
        # self.HConst = np.sqrt(self.Hmass * 1.6e-22/(2.*1.6e-12))
        # self.HeConst = np.sqrt(self.Hemass * 1.6e-22/(2.*1.6e-12))

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

        junkedStr = "Junked {:8d} inds due to {:14s} ({:8d} matches)"
        freshKeeps = np.ones(self.mlt.size, dtype=np.bool)

        if hasattr(self, 'quality_flag'):
            newKeeps = (self.quality_flag == 0)
            if verbose:
                print(junkedStr.format(
                    runningTall - newKeeps[newKeeps].size, "quality_flag", origTall-newKeeps[newKeeps].size))
                runningTall = newKeeps[newKeeps].size
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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
                print(junkedStr.format(
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
