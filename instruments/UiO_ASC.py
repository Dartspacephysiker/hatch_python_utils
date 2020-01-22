from bs4 import BeautifulSoup

import copy
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

allsitelist = ['lyr5', 'nya4', 'nya6', 'skn4', 'and3', 'lyr1', 'all']


def get_UiO_data_for_time(wantTime,
                          wantsite='lyr5',
                          checksites=None,
                          wantbl='6300',
                          maxtdelt=datetime.timedelta(minutes=10),
                          strictwant=False):

    boelgelengde = ['5577', '6300']
    if strictwant:
        boelgelengde = wantbl

    if checksites is None:
        sites = ['lyr5', 'nya4', 'nya6', 'skn4', 'and3', 'lyr1']
    else:
        sites = checksites
        try:
            _ = sites[0]
        except:
            sites = [sites]

    this = get_UiO_ASC_image_list_from_net(site=sites,
                                           date=wantTime,
                                           boelgelengde=boelgelengde,
                                           add_datetimes=True,
                                           return_dict=True)
    # imgdts = [that[0] for that in this]

    # roundConj = pd.Timestamp(wantTime).round('10min')
    # np.abs(this[0][0]-wantTime)

    if this is None:
        print("Couldn't get UiO nothing!")
        return None

    if len(this.keys()) == 0:
        print("Couldn't get UiO nothing!")
        return None

    if strictwant:
        if wantsite not in this.keys():
            print("Couldn't get UiO-{:s} data!".format(wantsite))
            return None

    UiOCloseDict = copy.deepcopy(this)

    # breakpoint()

    for site in this.keys():
        for bl in this[site].keys():
            for date in this[site][bl].keys():
                keeps = [(tim, np.abs(tim-wantTime), fila) for tim, fila in this[site]
                         [bl][date] if (np.abs(tim-wantTime) <= maxtdelt)]

                if len(keeps) > 0:
                    keepIkeep = -9
                    bestetdiff = datetime.timedelta(minutes=1000)
                    for ikeep, (tim, tdiff, fila) in enumerate(keeps):
                        if tdiff < bestetdiff:
                            bestetdiff = tdiff
                            keepIkeep = ikeep
                    if bestetdiff <= maxtdelt:
                        UiOCloseDict[site][bl] = keeps[keepIkeep]
                    else:
                        UiOCloseDict[site][bl] = None
                else:
                    UiOCloseDict[site][bl] = None

                this[site][bl][date] = keeps

        #         break
        #     break
        # break

    havesites = list(this.keys())

    havebls = []
    for site in havesites:
        blshere = this[site].keys()
        for bl in blshere:
            if bl not in havebls:
                havebls.append(bl)

    cangetimage = (len(havesites) > 0) and (len(havebls) > 0)
    wantsAreInClosedict = UiOCloseDict[wantsite][wantbl] is not None

    haveAnyCloseDictStuff = wantsAreInClosedict
    if not haveAnyCloseDictStuff:
        haveAnyCloseDictStuff = False
        doQuit = False
        for site in this.keys():
            for bl in this[site].keys():
                if UiOCloseDict[site][bl] is not None:
                    haveAnyCloseDictStuff = True
                    useUiObl = bl
                    useUiOsite = site
                    doQuit = True
                    break
                if doQuit:
                    break
            if doQuit:
                break

    if not (cangetimage and haveAnyCloseDictStuff):
        print("Couldn't get UiO image!")
        return None

    if cangetimage:
        if (wantsite in havesites) and wantsAreInClosedict:
            useUiOsite = wantsite
            print("Have {:s} UiO images!".format(wantsite))
        else:
            useUiOsite = havesites[0]
            print(
                "Don't have {:s} UiO images! Defaulting to {:s} ...".format(useUiOsite))

        if (wantbl in this[useUiOsite].keys()) and wantsAreInClosedict:
            useUiObl = wantbl
            print("Have {:s}-Å UiO images at {:s}!".format(wantbl, useUiOsite))
        else:
            # useUiObl = list(this[useUiOsite].keys())[0]
            print(
                "Don't have {:s}-Å UiO images at {:s}! Defaulting to {:s} ...".format(wantbl, useUiObl))

        uioCal = UiO_ASC_cal(site=useUiOsite,
                             date=wantTime,
                             boelgelengde=useUiObl,
                             fetch_from_net=True,
                             verbose=True)

        closedt, closedeltat, closefile = UiOCloseDict[useUiOsite][useUiObl]

        img = UiO_ASC_image(site=useUiOsite,
                            boelgelengde=useUiObl,
                            date=closedt.strftime("%Y%m%d"),
                            UTCHHMMSS=closedt.strftime("%H%M%S"),
                            verbose=True,
                            saveDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/')

        return dict(closedt=closedt, closedeltat=closedeltat, closefile=closefile,
                    cal=uioCal, img=img)


# def get_UiO_ASC_calfile_list_from_net(site='lyr5',
#                                       date='20190104',
#                                       boelgelengde='5577',
#                                       calAddr=None):

#     try:
#         if not isinstance(date, str):
#             nDates = len(date)
#             batchMode = True
#         else:
#             nDates = 1
#             date = [date]
#             batchMode = True
#     except:
#         batchMode = False

#     if not batchMode:

#         if calAddr is None:
#             if isinstance(date, datetime.datetime) or isinstance(date, datetime.date):
#                 date = date.strftime("%Y%m%d")

#             calAddr = '/'.join((UiOASCBaseAddr, site,
#                                 boelgelengde,
#                                 date[0:4], date))+'/'

#         availFiles = [filer for filer in get_url_paths(
#             calAddr) if filer.endswith('.dat')]

#     else:

#         availFiles = []
#         print("Getting list for {:d} dates ...".format(nDates))

#         for tmpdate in date:

#             # if calAddr is None:
#             if isinstance(tmpdate, datetime.datetime) \
#                or isinstance(tmpdate, datetime.date) \
#                or isinstance(tmpdate, np.ndarray):
#                 tmpdate = tmpdate.strftime("%Y%m%d")

#             calAddr = '/'.join((UiOASCBaseAddr, site,
#                                 boelgelengde,
#                                 tmpdate[0:4], tmpdate))+'/'

#             try:
#                 theseFiles = get_url_paths(calAddr)
#             except:
#                 print(
#                     "Couldn't get files for {:s}! Continuing ...".format(tmpdate))
#                 continue

#             availFiles += [
#                 filer for filer in theseFiles if filer.endswith('.dat')]

#     return availFiles


def get_UiO_ASC_image_list_from_net(site='lyr5',
                                    date='20190104',
                                    boelgelengde='5577',
                                    add_datetimes=False,
                                    return_dict=False):

    global allsitelist

    print("REMEMBER that you need to visit http://tid.uio.no/plasma/aurora and sign in FIRST; otherwise you'll never be able access the data...")

    try:
        if not isinstance(site, str):
            nSites = len(site)
        else:
            site = [site]
            nSites = 1

    except:
        site = [site]
        nSites = 1

    for sit in site:
        if sit not in allsitelist:
            print("No such thing as {:s}! please choose one of {:s} as a site!".format(
                sit,
                ", ".join(allsitelist)))
            return None

    try:
        if not isinstance(boelgelengde, str):
            nBLer = len(boelgelengde)
        else:
            boelgelengde = [boelgelengde]
            nBLer = 1

    except:
        boelgelengde = [boelgelengde]
        nBLer = 1

    blList = ['5577', '6300']
    for bl in boelgelengde:
        if bl not in blList:
            print("please choose one of {:s} as a boelgelengde!".format(
                ", ".join(blList)))
            return None

    batchMode = False
    try:
        if not isinstance(date, str):
            nDates = len(date)
            batchMode = True
        else:
            nDates = 1
            date = [date]
            batchMode = True
    except:
        print("Single date! {:s}".format(date.strftime("%Y-%m-%d")))
        date = [date]
        nDates = 1
        batchMode = True

    if not batchMode:

        return "BUNK"

    else:

        availFiles = []
        if return_dict:
            availFiles = {sit: {bl: {} for bl in boelgelengde} for sit in site}

        print("Getting image list for {:d} dates ...".format(nDates))

        for bl in boelgelengde:

            print("BL: {:s}".format(bl))
            print("#########")
            for tmpdate in date:

                # if calAddr is None:
                if isinstance(tmpdate, datetime.datetime) \
                   or isinstance(tmpdate, datetime.date) \
                   or isinstance(tmpdate, numpy.ndarray):
                    tmpdate = tmpdate.strftime("%Y%m%d")

                print("{:s}".format(tmpdate), end=' ')

                for sit in site:

                    print("{:s}".format(sit), end=' ')

                    remoteDir = '/'.join((UiOASCBaseAddr, sit,
                                          bl,
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
                                sit+"_"+tmpdate+"_([0-9]{2})([0-9]{2})([0-9]{2})_"+bl+'_cal.png')

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

                            if return_dict:
                                if tmpdate not in availFiles[sit][bl].keys():
                                    availFiles[sit][bl][tmpdate] = dt_file_list
                                else:
                                    availFiles[sit][bl][tmpdate] += dt_file_list
                            else:
                                availFiles += dt_file_list

                        else:

                            if return_dict:
                                if tmpdate not in availFiles[sit][bl].keys():
                                    availFiles[sit][bl][tmpdate] = [
                                        filer for filer in theseFiles if filer.endswith('.png')]
                                else:
                                    availFiles[sit][bl][tmpdate] += [
                                        filer for filer in theseFiles if filer.endswith('.png')]
                            else:
                                availFiles += [
                                    filer for filer in theseFiles if filer.endswith('.png')]

                print("")

    if return_dict:
        siteList = []
        for sit in site:
            keepsite = False
            for bl in boelgelengde:
                if len(availFiles[sit][bl].keys()) > 0:
                    keepsite = True
            if keepsite:
                siteList.append(sit)
        availFiles = {sit: availFiles[sit] for sit in siteList}

        # tmpkeepkeylist = []
        # for key in tmpdict.keys():

        # availFiles[bl][sit] = {k:v for k,v in tmpdict.items if vif len(tmpdict.keys()) > 0
        # removeMe =

    return availFiles


def download_UiO_ASC_image_list_from_net(availFiles,
                                         saveDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                                         verbose=False):

    try:
        doubleInder = len(availFiles[0]) > 1
    except:
        doubleInder = False

    if doubleInder:
        availFiles = [tup[1] for tup in availFiles]

    print("Getting {:d} files".format(len(availFiles)))

    # Download the file from `url` and save it locally under `file_name`:
    for remoteFile in availFiles:
        filNavn = remoteFile.split('/')[-1]

        localFile = saveDir + filNavn
        if not os.path.isfile(localFile):
            if verbose:
                print("Saving to " + localFile)

            urllib.request.urlretrieve(remoteFile, localFile)

    print("Done!")


class UiO_ASC_cal(object):

    def __init__(self,
                 *args,
                 site='lyr5',
                 date='20190104',
                 boelgelengde='5577',
                 calDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                 fetch_from_net=True,
                 verbose=False):

        if len(args) == 1:
            assert (isinstance(args[0], str) or
                    isinstance(args[0], datetime.date) or
                    isinstance(
                        args[0], datetime.datetime)), "Argument should be of type 'YYYYMMDD' or of type date/datetime!"

            date = args[0]

        # if isinstance(date, str):
        #     date = date
        if (isinstance(date, datetime.date) or isinstance(date, datetime.datetime)):
            date = date.strftime("%Y%m%d")

        calFile = site+'_'+date+'_'+boelgelengde+'_cal.dat'

        if not os.path.isfile(calDir+calFile):
            print(calFile+" not found!", end=' ')

            if not fetch_from_net:
                print(" Returning ...")
                return
            else:
                print("Checking UiO website ...")

                calAddr = '/'.join((UiOASCBaseAddr, site,
                                    boelgelengde,
                                    date[0:4], date))+'/'

                # availFiles = get_UiO_ASC_calfile_list_from_net(calAddr=calAddr)

                availFiles = [filer for filer in get_url_paths(
                    calAddr) if filer.endswith('.dat')]

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

        self.calFile = calFile
        self.calDir = calDir

        if 'site' in calStruc:
            self.site = calStruc['site'].decode("utf-8")
        if 'glats' in calStruc:
            self.glats = calStruc['glats']
        if 'glons' in calStruc:
            self.glons = calStruc['glons']
        if 'mlats' in calStruc:
            self.mlats = calStruc['mlats']
        if 'mlons' in calStruc:
            self.mlons = calStruc['mlons']
        if 'alts' in calStruc:
            self.alts = calStruc['alts']
        if 'gazms' in calStruc:
            self.gazms = calStruc['gazms']
        if 'mazms' in calStruc:
            self.mazms = calStruc['mazms']
        if 'elevs' in calStruc:
            self.elevs = calStruc['elevs']
        if 'x0' in calStruc:
            self.x0 = calStruc['x0']
        if 'y0' in calStruc:
            self.y0 = calStruc['y0']
        if 'r0' in calStruc:
            self.r0 = calStruc['r0']
        if 'angle' in calStruc:
            self.angle = calStruc['angle']
        if 'mirror_xaxis' in calStruc:
            self.mirror_xaxis = calStruc['mirror_xaxis']
        if 'conversion_factor' in calStruc:
            self.conversion_factor = calStruc['conversion_factor']
        if 'station' in calStruc:
            self.station = calStruc['station'].decode("utf-8")
        if 'location' in calStruc:
            self.location = calStruc['location'].decode("utf-8")
        if 'comment' in calStruc:
            self.comment = calStruc['comment'].decode("utf-8")
        if 'version' in calStruc:
            self.version = calStruc['version']
        if 'time_modified' in calStruc:
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

        self.get_UiO_image_times()

        # Check local dir
        if os.path.isfile(saveDir+fName):
            if verbose:
                print("Found file locally ... Loading!")

            self.data = imageio.imread(saveDir+fName)
            self.utc = datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                         int(self.UTCHHMMSS[0:2]), int(self.UTCHHMMSS[2:4]), int(self.UTCHHMMSS[4:6]))

            return

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

    # breakpoint()

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


def get_url_paths(url, ext='', params={}, debug=False):
    if debug:
        print("DEBUG: url = {:s}".format(url))
    response = requests.get(url, params=params, timeout=10)
    if response.ok:
        response_text = response.text
    else:
        return response.raise_for_status()

    soup = BeautifulSoup(response_text, 'html.parser')
    parent = [url + node.get('href') for node in soup.find_all('a')
              if node.get('href').endswith(ext)]
    return parent
