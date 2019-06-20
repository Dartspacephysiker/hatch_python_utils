import ftplib
import os.path
import fnmatch
from datetime import datetime

FASTFTPAddr = "cdaweb.gsfc.nasa.gov"


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
    for ftpFile in ftpFiles:

        localFile = localSaveDir+ftpFile.lstrip(subDir)
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
        # ftp.close()
        # return ftpFiles

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

    return ftpFiles
