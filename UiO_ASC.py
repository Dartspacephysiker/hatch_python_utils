from bs4 import BeautifulSoup

import datetime
import imageio
import numpy as np
import os
import requests
import scipy.io as sio
import urllib.request

UiOASCBaseAddr = 'http://tid.uio.no/plasma/aurora'


class UiO_ASC_cal(object):

    def __init__(self,
                 site='lyr5',
                 date='20190104',
                 boelgelengde='5577',
                 calDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                 fetch_from_net=True):

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

                availFiles = [filer for filer in get_url_paths(
                    calAddr) if filer.endswith('.dat')]

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


def get_UiO_ASC_image_from_net(site='lyr5',
                               date='20190104',
                               UTCHHMMSS='093600',
                               boelgelengde='5577',
                               verbose=True,
                               saveToLocal=True,
                               saveDir='/SPENCEdata/Research/database/Rockets/CAPER2/nordlyskamera/',
                               lookForNearestTime=True):
    """
    site        : 'lyr5' or 'nya6'
    date        : Format YYYYMMDD
    UTCHHMMSS   : Format HHMMSS
    boelgelengde: '5577' or '6300'
    """
    fName = '_'.join((site, date, UTCHHMMSS, boelgelengde, 'cal.png'))

    # Check local dir
    if os.path.isfile(saveDir+fName):
        print("Found file locally ... Loading!")
        return imageio.imread(saveDir+fName)

    imDir = '/'.join((UiOASCBaseAddr, site,
                      boelgelengde,
                      date[0:4], date, 'ut'+UTCHHMMSS[0:2]))+'/'

    imAddr = imDir + fName
    # imAddr = '/'.join((UiOASCBaseAddr,site,
    #                    boelgelengde,
    #                    date[0:4],date,'ut'+UTCHHMMSS[0:2],
    #                    fName))

    availFiles = [filer for filer in get_url_paths(
        imDir) if filer.endswith('.png')]

    if imAddr in availFiles:
        print("Found "+imAddr)
    elif lookForNearestTime:
        print(imAddr+' not found!')
        print('Finding nearest ... ', end='')

        killStrL = imDir+site+'_'+date+'_'
        killStrR = '_'+boelgelengde+'_cal.png'
        pngTimes = [filer.replace(killStrL, '').replace(
            killStrR, '') for filer in availFiles]

        this = [(int(time[0:2]), int(time[2:4]), int(time[4:6]))
                for time in pngTimes if (len(time) == 6)]
        dtlist = np.array([datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                             these[0], these[1], these[2]) for these in this])

        dtWant = datetime.datetime(int(date[0:4]), int(date[4:6]), int(date[6:8]),
                                   int(UTCHHMMSS[0:2]), int(UTCHHMMSS[2:4]), int(UTCHHMMSS[4:6]))

        closest = np.argmin(np.abs(np.array(dtlist)-dtWant))

        print("Closest is {:s}".format(
            dtlist[closest].strftime("%Y%m%d/%H:%M:%S")))

        imAddr = availFiles[closest]

    else:
        print("Didn't find "+imAddr+'! Returning ...')
        return None

    if verbose:
        print("imDir : " + imDir)
        print("imAddr: " + imAddr)

    if saveToLocal:
        # Download the file from `url` and save it locally under `file_name`:
        localFile = saveDir+imAddr.replace(imDir, '')
        if not os.path.isfile(localFile):

            print("Saving to " + localFile)
            urllib.request.urlretrieve(imAddr, localFile)

        else:
            print("Already have locally!")

    return imageio.imread(imAddr)


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
