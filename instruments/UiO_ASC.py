from bs4 import BeautifulSoup

import datetime
import imageio
import numpy as np
import os
import requests
import scipy.io as sio
import urllib.request
import time
import re

UiOASCBaseAddr = 'http://tid.uio.no/plasma/aurora'


def get_UiO_ASC_calfile_list_from_net(site='lyr5',
                                      date='20190104',
                                      boelgelengde='5577',
                                      calAddr=None):

    try:
        nDates = len(date)
        batchMode = True
    except:
        batchMode = False

    if not batchMode:

        if calAddr is None:
            if isinstance(date, datetime.datetime) or isinstance(date, datetime.date):
                date = date.strftime("%Y%m%d")

            calAddr = '/'.join((UiOASCBaseAddr, site,
                                boelgelengde,
                                date[0:4], date))+'/'

        availFiles = [filer for filer in get_url_paths(
            calAddr) if filer.endswith('.dat')]

    else:

        availFiles = []
        print("Getting list for {:d} dates ...".format(nDates))

        for tmpdate in date:

            # if calAddr is None:
            if isinstance(tmpdate, datetime.datetime) \
               or isinstance(tmpdate, datetime.date) \
               or isinstance(tmpdate, numpy.ndarray):
                tmpdate = tmpdate.strftime("%Y%m%d")

            calAddr = '/'.join((UiOASCBaseAddr, site,
                                boelgelengde,
                                tmpdate[0:4], tmpdate))+'/'

            try:
                theseFiles = get_url_paths(calAddr)
            except:
                print(
                    "Couldn't get files for {:s}! Continuing ...".format(tmpdate))
                continue

            availFiles += [
                filer for filer in theseFiles if filer.endswith('.dat')]

    return availFiles


def get_UiO_ASC_image_list_from_net(site='lyr5',
                                    date='20190104',
                                    boelgelengde='5577',
                                    add_datetimes=False):

    try:
        nDates = len(date)
        batchMode = True
    except:
        batchMode = False

    if not batchMode:

        return "BUNK"

    else:

        availFiles = []
        print("Getting image list for {:d} dates ...".format(nDates))

        for tmpdate in date:

            # if calAddr is None:
            if isinstance(tmpdate, datetime.datetime) \
               or isinstance(tmpdate, datetime.date) \
               or isinstance(tmpdate, numpy.ndarray):
                tmpdate = tmpdate.strftime("%Y%m%d")

            print("{:s}".format(tmpdate), end=' ')

            remoteDir = '/'.join((UiOASCBaseAddr, site,
                                  boelgelengde,
                                  tmpdate[0:4], tmpdate))+'/'
            try:
                blig = get_url_paths(remoteDir)
            except:
                print(
                    "Couldn't get files! Continuing ...".format(tmpdate))
                continue

            dirs = [dirr for dirr in blig if dirr[-5:-3] == 'ut']

            if len(dirs) == 0:
                print(
                    "No directories! Continuing ...".format(tmpdate))

            for dirr in dirs:

                time.sleep(0.1)

                try:
                    theseFiles = get_url_paths(dirr)
                except:
                    print(
                        "Couldn't get files for {:s}! Continuing ...".format(dirr))
                    continue

                if add_datetimes:

                    myRE = re.compile(
                        site+"_"+tmpdate+"_([0-9]{2})([0-9]{2})([0-9]{2})_"+boelgelengde+'_cal.png')

                    dt_file_list = []
                    for filer in theseFiles:

                        this = myRE.search(filer.split("/")[-1])
                        if this is None:
                            continue

                        hr, mt, sec = this.groups()

                        dt_file_list.append((datetime.datetime(int(tmpdate[0:4]),
                                                               int(tmpdate[4:6]),
                                                               int(tmpdate[6:8]),
                                                               int(hr), int(mt), int(sec)),
                                             filer))

                    availFiles += dt_file_list

                else:

                    availFiles += [
                        filer for filer in theseFiles if filer.endswith('.png')]

            print("")

    return availFiles


class UiO_ASC_cal(object):

    def __init__(self,
                 *args,
                 site='lyr5',
                 date='20190104',
                 boelgelengde='5577',
                 calDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                 fetch_from_net=True):

        if len(args) == 1:
            assert (isinstance(args[0], str) or
                    isinstance(args[0], datetime.date) or
                    isinstance(
                        args[0], datetime.datetime)), "Argument should be of type 'YYYYMMDD' or of type date/datetime!"

            if isinstance(args[0], str):
                date = args[0]
            elif (isinstance(args[0], datetime.date) or isinstance(args[0], datetime.datetime)):
                date = args[0].strftime("%Y%m%d")

        calFile = site+'_'+date+'_'+boelgelengde+'_cal.dat'

        if not os.path.isfile(calDir+calFile):
            print(calFile+" not found!", end='')

            if not fetch_from_net:
                print(" Returning ...")
                return
            else:
                print("Checking UiO website ...")

                calAddr = '/'.join((UiOASCBaseAddr, site,
                                    boelgelengde,
                                    date[0:4], date))+'/'

                availFiles = get_calfile_list_from_net(calAddr=calAddr)

                # availFiles = [filer for filer in get_url_paths(
                #     calAddr) if filer.endswith('.dat')]

                # if list_from_net:
                #     return availFiles

                if (len(availFiles) == 0) or (not (calAddr+calFile) in availFiles):
                    print("Couldn't find cal file ...")
                    return
                else:
                    print("Found cal file! Downloading to " + calDir + ' ...')

                    # Download the file from `url` and save it locally under `file_name`:
                    urllib.request.urlretrieve(calAddr+calFile, calDir+calFile)

        calStruc = sio.readsav(calDir+calFile)

        self.site = calStruc['site'].decode("utf-8")
        self.glats = calStruc['glats']
        self.glons = calStruc['glons']
        self.mlats = calStruc['mlats']
        self.mlons = calStruc['mlons']
        self.alts = calStruc['alts']
        self.gazms = calStruc['gazms']
        self.mazms = calStruc['mazms']
        self.elevs = calStruc['elevs']
        self.x0 = calStruc['x0']
        self.y0 = calStruc['y0']
        self.r0 = calStruc['r0']
        self.angle = calStruc['angle']
        self.mirror_xaxis = calStruc['mirror_xaxis']
        self.conversion_factor = calStruc['conversion_factor']
        self.station = calStruc['station'].decode("utf-8")
        self.location = calStruc['location'].decode("utf-8")
        self.comment = calStruc['comment'].decode("utf-8")
        self.version = calStruc['version']
        self.time_modified = calStruc['time_modified'].decode("utf-8")
        self.date = date
        self.utc = datetime.date(
            int(date[0:4]), int(date[4:6]), int(date[6:8]))

    def __repr__(self):
        attrs = [a for a in dir(self) if not a.startswith('__')]

        def stringOutput(a):
            if hasattr(getattr(self, a), 'shape'):
                if getattr(self, a).shape:
                    return "{:20s}: {:10s} {:10s}".format(a,
                                                          str(getattr(
                                                              self, a).shape),
                                                          str(getattr(self, a).dtype))
                else:
                    return "{:20s}: {:10s} {:10s}".format(a,
                                                          "{:.2f}".format(
                                                              getattr(self, a)),
                                                          str(getattr(self, a).dtype))

            elif isinstance(getattr(self, a), str):
                return "{:20s}: {:10s} {:10s}".format(a,
                                                      getattr(self, a),
                                                      str(type(getattr(self, a))))
            elif isinstance(getattr(self, a), datetime.datetime):
                return "{:20s}: {:10s} {:10s}".format(a,
                                                      str(getattr(self, a)),
                                                      str(type(getattr(self, a))))
            else:
                return "{:20s}: {:10s} {:10s}".format(a,
                                                      'scalar',
                                                      str(type(getattr(self, a))))

        # fullStrs = ["{:20s}: {:10s} {:10s}".format(a,
        #                                            str(getattr(self,a).shape),
        #                                            str(getattr(self,a).dtype)) for a in attrs]
        fullStrs = [stringOutput(a) for a in attrs]

        # fullStrs = ["{:s}: {:s}".format(getattr(self,a),str(getattr(self,a).shape)) for a in attrs]
        return "\n".join(fullStrs)


class UiO_ASC_image(object):

    def __init__(self, *args,
                 site='lyr5',
                 date='20190104',
                 UTCHHMMSS='093600',
                 boelgelengde='5577',
                 verbose=False,
                 saveToLocal=True,
                 saveDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                 lookForNearestTime=True):

        getHHMMSSFromArg0OrKeyword = True
        getDateFromArg0OrKeyword = True
        if len(args) == 2:
            # args[0] = date

            assert (isinstance(
                args[1], str)), "Second argument should be UTC string 'HHMMSS' (or provide datetime as first arg)!"
            # isinstance(args[1], datetime.date) or
            # isinstance(
            #     args[1], datetime.datetime)), "Second argument should be of type 'YYYYMMDD' or of type date/datetime!"

            self.UTCHHMMSS = args[1]
            getHHMMSSFromArg0OrKeyword = False

        if len(args) >= 1:

            assert (isinstance(args[0], str) or
                    isinstance(args[0], datetime.date) or
                    isinstance(
                        args[0], datetime.datetime)), "First argument should be of type 'YYYYMMDD' or of type date/datetime!"

            if isinstance(args[0], str):
                self.date = args[0]
            elif (isinstance(args[0], datetime.date) or isinstance(args[0], datetime.datetime)):
                self.date = args[0].strftime("%Y%m%d")

            getDateFromArg0OrKeyword = False

            if getHHMMSSFromArg0OrKeyword:
                if isinstance(args[0], datetime.datetime):
                    self.UTCHHMMSS = args[0].strftime("%H%M%S")
                else:
                    self.UTCHHMMSS = '060000'

                getHHMMSSFromArg0OrKeyword = False

                # self.data = get_UiO_ASC_image_from_net(site=site,
                #                                        date=date,
                #                                        UTCHHMMSS=UTCHHMMSS,
                #                                        boelgelengde=boelgelengde,
                #                                        verbose=verbose,
                #                                        saveToLocal=saveToLocal,
                #                                        saveDir=saveDir,
                #                                        lookForNearestTime=lookForNearestTime)

        self.saveDir = saveDir
        self.site = site

        if getDateFromArg0OrKeyword:
            self.date = date
        if getHHMMSSFromArg0OrKeyword:
            self.UTCHHMMSS = UTCHHMMSS

        self.boelgelengde = boelgelengde

        self.availFiles = "No websearch performed"
        self.nfiles = 0

        self.remoteDir = '/'.join((UiOASCBaseAddr, self.site,
                                   boelgelengde,
                                   self.date[0:4], self.date, 'ut'+self.UTCHHMMSS[0:2]))+'/'

        fName = '_'.join((self.site, self.date, self.UTCHHMMSS,
                          self.boelgelengde, 'cal.png'))
        self.request_fName = fName

        self.remoteFile = self.remoteDir + fName

        # Check local dir
        if os.path.isfile(saveDir+fName):
            if verbose:
                print("Found file locally ... Loading!")

            self.data = imageio.imread(saveDir+fName)
            self.utc = datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                         int(self.UTCHHMMSS[0:2]), int(self.UTCHHMMSS[2:4]), int(self.UTCHHMMSS[4:6]))

            return

        self.get_UiO_image_times()

        dtWant = datetime.datetime(int(self.date[0:4]), int(self.date[4:6]), int(self.date[6:8]),
                                   int(self.UTCHHMMSS[0:2]), int(self.UTCHHMMSS[2:4]), int(self.UTCHHMMSS[4:6]))

        if len(self.availFiles) == 0:
            self.data = None
            return

        if self.remoteFile in self.availFiles:
            if verbose:
                print("Found "+self.remoteFile)

            closest = np.argmin(np.abs(np.array(self.dtlist)-dtWant))
            self.utc = dtWant

        elif lookForNearestTime:

            if verbose:
                print(self.remoteFile+' not found!')
                print('Finding nearest ... ', end='')

            closest = np.argmin(np.abs(np.array(self.dtlist)-dtWant))

            if verbose:
                print("Closest is {:s}".format(
                    self.dtlist[closest].strftime("%Y%m%d/%H:%M:%S")))

            self.remoteFile = self.availFiles[closest]
            self.utc = self.dtlist[closest]

        else:
            if verbose:
                print("Didn't find "+self.remoteFile+'! Returning ...')
            return None

        if verbose:
            print("remoteDir : " + self.remoteDir)
            print("remoteFile: " + self.remoteFile)

        if saveToLocal:
            # Download the file from `url` and save it locally under `file_name`:
            self.localFile = saveDir + \
                self.remoteFile.replace(self.remoteDir, '')
            if not os.path.isfile(self.localFile):

                if verbose:
                    print("Saving to " + self.localFile)

                urllib.request.urlretrieve(self.remoteFile, self.localFile)

            else:
                if verbose:
                    print("Already have locally!")

        self.data = imageio.imread(self.remoteFile)

    def get_UiO_image_times(self):
        self.availFiles = [filer for filer in get_url_paths(
            self.remoteDir) if filer.endswith('.png')]

        killStrL = self.remoteDir+self.site+'_'+self.date+'_'
        killStrR = '_'+self.boelgelengde+'_cal.png'
        pngTimes = np.array([filer.replace(killStrL, '').replace(
            killStrR, '') for filer in self.availFiles])

        keepFiles = [(len(time) == 6) for time in pngTimes]
        if np.where(keepFiles)[0].size == 0:
            print("Couldn't find any valid files on the 'net! Returning ...")
            self.availFiles = []
            self.nfiles = 0
            return
        else:
            pngTimes = pngTimes[keepFiles]
            self.availFiles = [self.availFiles[i]
                               for i in range(len(keepFiles)) if keepFiles[i]]
            self.nfiles = len(self.availFiles)

        this = [(int(time[0:2]), int(time[2:4]), int(time[4:6]))
                for time in pngTimes]
        self.dtlist = np.array([datetime.datetime(int(self.date[0:4]), int(self.date[4:6]), int(self.date[6:8]),
                                                  these[0], these[1], these[2]) for these in this])


def get_UiO_ASC_image_from_net(site='lyr5',
                               date='20190104',
                               UTCHHMMSS='093600',
                               boelgelengde='5577',
                               verbose=True,
                               saveToLocal=True,
                               saveDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                               lookForNearestTime=True):
    """
    site               : 'lyr5' or 'nya6'
    date               : Format YYYYMMDD
    UTCHHMMSS          : Format HHMMSS
    boelgelengde       : '5577' or '6300'
    savetoLocal        : Save image from UiO website to "saveDir"
    lookForNearestTime : If time specified by UTCHHMMSS not found, search UiO 'date' directory for nearest time
    """
    fName = '_'.join((site, date, UTCHHMMSS, boelgelengde, 'cal.png'))

    # Check local dir
    if os.path.isfile(saveDir+fName):
        print("Found file locally ... Loading!")
        return imageio.imread(saveDir+fName)

    remoteDir = '/'.join((UiOASCBaseAddr, site,
                          boelgelengde,
                          date[0:4], date, 'ut'+UTCHHMMSS[0:2]))+'/'

    remoteFile = remoteDir + fName
    # remoteFile = '/'.join((UiOASCBaseAddr,site,
    #                    boelgelengde,
    #                    date[0:4],date,'ut'+UTCHHMMSS[0:2],
    #                    fName))

    availFiles = [filer for filer in get_url_paths(
        remoteDir) if filer.endswith('.png')]

    if remoteFile in availFiles:
        print("Found "+remoteFile)
    elif lookForNearestTime:
        print(remoteFile+' not found!')
        print('Finding nearest ... ', end='')

        killStrL = remoteDir+site+'_'+date+'_'
        killStrR = '_'+boelgelengde+'_cal.png'
        pngTimes = np.array([filer.replace(killStrL, '').replace(
            killStrR, '') for filer in availFiles])

        keepFiles = [(len(time) == 6) for time in pngTimes]
        if np.where(keepFiles)[0].size == 0:
            print("Couldn't find any valid files on the 'net! Returning ...")
            return None
        else:
            pngTimes = pngTimes[keepFiles]
            availFiles = [availFiles[i]
                          for i in range(len(keepFiles)) if keepFiles[i]]

        this = [(int(time[0:2]), int(time[2:4]), int(time[4:6]))
                for time in pngTimes]
        dtlist = np.array([datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                             these[0], these[1], these[2]) for these in this])

        dtWant = datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                   int(UTCHHMMSS[0:2]), int(UTCHHMMSS[2:4]), int(UTCHHMMSS[4:6]))

        closest = np.argmin(np.abs(np.array(dtlist)-dtWant))

        print("Closest is {:s}".format(
            dtlist[closest].strftime("%Y%m%d/%H:%M:%S")))

        remoteFile = availFiles[closest]

    else:
        print("Didn't find "+remoteFile+'! Returning ...')
        return None

    if verbose:
        print("remoteDir : " + remoteDir)
        print("remoteFile: " + remoteFile)

    if saveToLocal:
        # Download the file from `url` and save it locally under `file_name`:
        localFile = saveDir+remoteFile.replace(remoteDir, '')
        if not os.path.isfile(localFile):

            print("Saving to " + localFile)
            urllib.request.urlretrieve(remoteFile, localFile)

        else:
            print("Already have locally!")

    return imageio.imread(remoteFile)


def get_url_paths(url, ext='', params={}):
    response = requests.get(url, params=params)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()
    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a')
              if node.get('href').endswith(ext)]
    return parent
