import numpy as np
import pandas as pd
from datetime import datetime 
import gzip

from hatch_python_utils.hatch_utils import get_basedir
from getpass import getpass
import os
import urllib
import requests

def load_cluster_parquet(clusterparquets):
    """
    Load cluster parquet files, gjerne utdata fra download_process_save_cluster_lobe_file()
    """

    if isinstance(clusterparquets,str):
        clusterparquets = [clusterparquets]

    dfs = []
    for clusterparquet in clusterparquets:
        dfs.append(pd.read_parquet(clusterparquet))
    dfs = pd.concat(dfs)

    return dfs


def download_process_save_cluster_lobe_file(
        sat='C1',
        years=[],
        filetype='A',
        verbose=False,
        overwrite_existing=False,
        only_local=True,
        remove_txtgz_files=True):
    
    if (not hasattr(years,'__len__')) or isinstance(years,str):
        years = [years]
    
    assert all([isinstance(year,int) for year in years]), \
        "Must provide year as four-digit int between 2001 and 2020!"
        
    assert sat in ['C1','C2','C3','C4']
    assert all([year in range(2001,2021) for year in years])
    assert filetype in ['A','C']

    basedir, isColtrane = get_basedir()	
    datadir = basedir+'database/Cluster/lobe/'
    polarcapdbdirREMOTE = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/lobe_ver_2020-03-23/'

    if not only_local:
        user = getpass('Username: ')
        pw = getpass('PW: ')

    # polarcapdbfileREMOTE = 'Lobe_C1_2001_A_med1min.txt.gz'
    polarcapdbfilesREMOTE = [f"Lobe_{sat:s}_{year:d}_{filetype}_med1min.txt.gz" for year in years]
    
    outparquets = []
    for polarcapdbfileREMOTE in polarcapdbfilesREMOTE:
        url = polarcapdbdirREMOTE+polarcapdbfileREMOTE
        filename = os.path.basename(urllib.parse.urlparse(url)[2])
        parquetout = filename.replace('.txt','.parquet').replace('.gz','')

        if os.path.exists(datadir+parquetout) and (not overwrite_existing):
            if not only_local:
                print(f"Already have {parquetout:s}, not overwriting!")

            outparquets.append(datadir+parquetout)
            continue

        # Skip if only local files wanted
        if only_local:
            continue

        r = requests.get(url, auth=(user,pw))
        
        if not (r.status_code == 200):
            print(f"Zing! Couldn't get {polarcapdbfileREMOTE}!")
            continue
        else:
           with open(datadir+filename, 'wb') as out:
              for bits in r.iter_content():
                  _ = out.write(bits)
        
        df = read_cluster_txt_file(datadir+filename,
                              verbose=verbose)
        
        print(f"Saving to {parquetout:s}")
        df.to_parquet(datadir+parquetout)
        outparquets.append(datadir+parquetout)

    if remove_txtgz_files:
        toberemoved = []

        for polarcapdbfileREMOTE in polarcapdbfilesREMOTE:
            url = polarcapdbdirREMOTE+polarcapdbfileREMOTE
            filename = os.path.basename(urllib.parse.urlparse(url)[2])

            if os.path.exists(datadir+filename):
                toberemoved.append(datadir+filename)

        if len(toberemoved) > 0:
            print("Removing .txt.gz files")
            for toremove in toberemoved:
                os.remove(toremove)

    return outparquets


def parse_cluster_txt_header(headerline,verbose=False):

    from itertools import groupby
    
    # Get column names
    namedict = dict()
    positiondict = dict()
    positioncount = 0
    for k, g in groupby(enumerate(headerline), lambda x: not x[1].isspace()):
        if k:
    
            pos, first_item = next(g)
            name = (first_item + ''.join([x for _, x in g])).replace('#','')
            if verbose:
                print(pos, name)
    
            namedict[positioncount] = name
            positiondict[positioncount] = pos
    
            positioncount += 1

    return namedict,positiondict


def read_cluster_txt_file(clusterfile,
                          verbose=False):
    linecount = 0
    dataz = []
    headerpos = -1
    haveHeader = False
    headerstr = '#_ISO_Time'

    dtfmt = '%Y-%m-%dT%H:%M:%S.%f'
    integerparams = [1,12,32,66]

    if '.gz' in clusterfile[-4:]:
        openfunc = lambda f: gzip.open(f,mode='rt')
        if verbose:
            print("Opening with gzip ...")
    else:
        openfunc = open

    with openfunc(clusterfile) as f:
    
        for line in f:
    
            if haveHeader:
                if line[0] == '#':
                    if verbose:
                        print(f"End of file at line {linecount:d}")
                    break
    
                dataz.append(line.split())
    
            else:
                if headerstr in line:
                    if verbose:
                        print(f"Found header at line {linecount:d}")
                    headerpos = linecount
                    
                    headerline = line[:]
                    haveHeader = True
                    
                    linecount += 1
                    
                    units = next(f).split()
                    
                    namedict,positiondict = parse_cluster_txt_header(headerline,
                                                                     verbose=verbose)
    
            linecount += 1
    
    datadict = {'Time':[datetime.strptime(dataz[i][0],dtfmt) for i in range(len(dataz))]}    
    
    for k in range(1,len(namedict.keys())):
        if k in integerparams:
            if verbose:
                print(f"INTEGER: {k:2d} {namedict[k]:s}")
            datadict[namedict[k]] = np.array([int(dataz[i][k]) for i in range(len(dataz))]).astype(np.int32)
        else:
            if verbose:
                print(f"FLOAT  : {k:2d} {namedict[k]:s}")
    
            datadict[namedict[k]] = np.array([float(dataz[i][k]) for i in range(len(dataz))])
            
    df = pd.DataFrame(datadict)
    df.set_index('Time',inplace=True)

    return df


