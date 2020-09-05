#!/usr/bin/env python
# coding: utf-8

########################################
# SHEET 1     (20191025): Read conjunction lists from https://sscweb.gsfc.nasa.gov/cgi-bin/Query.cgi, make dataframe
########################################

from datetime import datetime,timedelta
import importlib
import numpy as np
import pandas as pd

from hatch_python_utils.hatch_utils import get_basedir
baseDir,isColtrane = get_basedir()

from hatch_python_utils.earth.utils import sphDist

import apexpy
from pytt.earth import geodesy

import hatch_python_utils
from hatch_python_utils.earth import coordinates as hCoord
junk = importlib.reload(hatch_python_utils)
junk = importlib.reload(hCoord)


def read_nasa_sscweb_query_output(txtfiles,debug=False):
    """
    Read output from a conjunction list generated by NASA Satellite Situation Center's Query tool* and output a pandas DataFrame

    NOTE: This will probably not work without modification based on the structure of your query output.
    In particular I assume this sort of header:

    "       Time       Sat.   |     GEO     | Radius|  Trace GEO    ArcLen|Dst(lead)-(name"
    "yy/mm/dd hh:mm:ss        |  Lat   Long |  (km) |  Lat   Long    (km) | Dist.  Name   "


    * https://sscweb.gsfc.nasa.gov/cgi-bin/Query.cgi

    """

    apex_refh_km = 110

    numerics = ['glat','glon','r_km','tracelat','tracelon','arclenkm','dist']
    datetimes = ['date']

    if isinstance(txtfiles,str):
        txtfiles = [txtfiles]

    # doSave = True
    # outdir = '/SPENCEdata/Research/database/Swarm/konjunksjoner/'
    # outcsv = 'Swarm'+leadSatLetter+'B_konj_'+datesuff+'.csv'

    df = pd.DataFrame(columns=['date','satellite','glat','glon','r_km','tracelat','tracelon','arclenkm','dist','leadname'])
    keepcount = 0

    totcnt = 0
    for infile in txtfiles:
        cnt = 0
        datacnt = 0
        with open(infile) as fp:

            for line in fp:

                # line = fp.readline()

                cnt += 1
                totcnt += 1

                if debug:
                    print("Line {}: {}".format(cnt, line.strip()),end='')
                else:
                    if (cnt % 100) == 0:
                        print("Line {}: {}".format(cnt, line.strip()))

                try:
                    # yr = int(line[0:2])
                    yr,mo,day = [int(ta) for ta in line[:9].split('/')]

                except:
                    # print("noparse! Skipt!")
                    if debug:
                        print(" [SKIP]")
                    continue

                if debug:
                    print(" [DATA]")

                parseline = [this for this in line.rstrip().split(' ') if (this != '')]
                linelen = len(parseline)
                if linelen == 11:
                    date, timestr, satellite, glat, glon, r_km, tracelat, tracelon, arclenkm, dist,leadname = parseline
                    leadSatName = leadname
                    prevsat = satellite
                elif linelen == 10:
                    date, timestr, glat, glon, r_km, tracelat, tracelon, arclenkm, dist,leadname = parseline
                    satellite = prevsat
                elif linelen == 9:
                    date, timestr, satellite, glat, glon, r_km, tracelat, tracelon, arclenkm = parseline
                    dist = '0'
                    leadname = leadSatName
                    prevsat = satellite
                elif linelen == 8:
                    date, timestr, glat, glon, r_km, tracelat, tracelon, arclenkm = parseline
                    dist = '0'
                    leadname = leadSatName
                    satellite = prevsat
                else:
                    assert 2<0,"AARRRG"

                if yr <= 30:     # For years greater than 30, assume we're talking about the 1900s
                    centuryMinus1 = '20'
                else:
                    centuryMinus1 = '19'

                dt = datetime.strptime(centuryMinus1+date+' '+timestr,"%Y/%m/%d %H:%M:%S")
                glat = float(glat)
                glon = float(glon)
                r_km = int(r_km)
                tracelat = float(tracelat)
                tracelon = float(tracelon)
                arclenkm = int(arclenkm)
                dist = int(dist.replace('\n',''))

                df.loc[keepcount] = [dt,satellite,glat,glon,r_km,tracelat,tracelon,arclenkm,dist,leadname.replace('\n','')]

                keepcount += 1

                # break
                # if cnt > 93:
                #     break

    for num in numerics:
        df[num] = pd.to_numeric(df[num])
    for dat in datetimes:
        df[dat] = pd.to_datetime(df[dat])

    # df = df.astype({'r_km': 'int64','arclenkm':'int64'}).copy()

    gdlat = geodesy.geocentric2geodeticlat(df['glat'].values)
    gdlat2, heights = hCoord.geoclatR2geodlatheight(df['glat'].values,df['r_km'].values)
    df['alt'] = heights
    df['gdlat'] = gdlat

    mlats = []
    mlons = []
    mlts = []
    for i in range(df.shape[0]):
        # print(i)
        a = apexpy.Apex(date=df.iloc[i].date,refh=apex_refh_km)
        mlat, mlon = a.geo2apex(df.iloc[i]['gdlat'],
                                df.iloc[i]['glon'],
                                df.iloc[i]['alt'])

        mlt = a.mlon2mlt(mlon,df.iloc[i].date)

        mlats.append(mlat)
        mlons.append(mlon)
        mlts.append(mlt)

    df['mlat'] = mlats
    df['mlon'] = mlons
    df['mlt'] = mlts

    # df.mlt.hist()
    df.set_index('date',inplace=True)
    df.sort_index(inplace=True)

    return df

    # print("Maybe what you'd like to do is make a multi-index, with timestamp as index0 and not-lead satellite as index1")

    # dfA = df[df['satellite'] == leadSatName].copy()
    # dfB = df[df['satellite'] != leadSatName].copy()

    # dfA.drop(labels='satellite',axis=1,inplace=True)
    # dfB.drop(labels='satellite',axis=1,inplace=True)

    # combo = dfA.join(dfB,lsuffix=leadSatLetter,rsuffix='B')

    # # get sphdist
    # combo['sphDist_mag'] = sphDist(combo['mlat'+leadSatLetter], combo['mlon'+leadSatLetter], combo['mlatB'], combo['mlonB'],
    #                                mltMode=False)
    # combo['sphDist_geoc'] = sphDist(combo['glat'+leadSatLetter], combo['glon'+leadSatLetter], combo['glatB'], combo['glonB'],
    #                                mltMode=False)
    # combo['sphDist_geod'] = sphDist(combo['gdlat'+leadSatLetter], combo['glon'+leadSatLetter], combo['gdlatB'], combo['glonB'],
    #                                mltMode=False)


    # return combo

    # if doSave:
    #     print("Saving to {:s}".format(outcsv))
    #     combo.to_csv(outdir+outcsv)