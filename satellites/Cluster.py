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
    Load cluster parquet files, 
        gjerne utdata fra download_process_save_cluster_{lobe,polarcap,all}_file()
    """

    if isinstance(clusterparquets,str):
        clusterparquets = [clusterparquets]

    dfs = []
    for clusterparquet in clusterparquets:
        dfs.append(pd.read_parquet(clusterparquet))
    dfs = pd.concat(dfs)

    return dfs


def download_process_save_cluster_file_byregion(**kws):
    
    assert 'region' in kws,"Use filetype == [one of 'lobe','polarcap','all']"

    if kws['region'].lower() == 'lobe':
        print("Lobe-region files ...")
        outparquets = download_process_save_cluster_lobe_file(**kws)
    elif kws['region'].lower() == 'polarcap':
        print("Polar cap-region files ...")
        outparquets = download_process_save_cluster_polarcap_file(**kws)
    elif kws['region'].lower() == 'all':
        print("All-region files ...")
        outparquets = download_process_save_cluster_all_file(**kws)
    elif kws['region'].lower() == 'hires':
        print("Hi-res files ...")
        outparquets = download_process_save_cluster_hires_file(**kws)

    return outparquets 

def download_process_save_cluster_hires_file(**kws):
    
    basedir, isColtrane = get_basedir()	
    datadir = basedir+'database/Cluster/hires/'
    remotedir = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/hires_no_filter_ver_2020-03-23/'

    filepref,filesuff = 'All_','_hires.txt.xz'
    clusterfiles = get_cluster_filenames(filepref,filesuff,**kws)

    outparquets = _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,**kws)

    return outparquets


def download_process_save_cluster_lobe_file(**kws):
    
    basedir, isColtrane = get_basedir()	
    datadir = basedir+'database/Cluster/lobe/'
    remotedir = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/lobe_ver_2020-03-23/'

    filepref,filesuff = 'Lobe_','_med1min.txt.gz'
    clusterfiles = get_cluster_filenames(filepref,filesuff,**kws)

    outparquets = _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,**kws)

    return outparquets


def download_process_save_cluster_polarcap_file(**kws):
    
    basedir, isColtrane = get_basedir()	
    datadir = basedir+'database/Cluster/polarcap/'
    remotedir = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/polar_cap_ver_2020-03-23/'

    filepref,filesuff = 'PolarCap_','_med1min.txt.gz'
    clusterfiles = get_cluster_filenames(filepref,filesuff,**kws)

    outparquets = _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,**kws)

    return outparquets


def download_process_save_cluster_all_file(**kws):
    
    basedir, isColtrane = get_basedir()	
    datadir = basedir+'database/Cluster/all/'
    remotedir = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/all_no_filter_ver_2020-03-23/'

    filepref,filesuff = 'All_','_med1min.txt.gz'
    clusterfiles = get_cluster_filenames(filepref,filesuff,**kws)

    outparquets = _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,**kws)

    return outparquets


def download_process_save_cluster_sum_file(**kws):
    
    location = kws['location'] if 'location' in kws else 'All'

    basedir, isColtrane = get_basedir()	

    baseremotedir = 'http://fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/current_version/'
    if 'all' in location.lower():
        remotedir = baseremotedir+'all_no_filter_ver_2020-03-23/'
        filepref,filesuff = 'All_','_med1min.txt.gz'
        
    elif 'polarcap' in location.lower():
        remotedir = baseremotedir+'polar_cap_ver_2020-03-23/'
        filepref,filesuff = 'PolarCap_','_med1min.txt.gz'
    elif 'lobe' in location.lower():
        remotedir = baseremotedir+'lobe_ver_2020-03-23/'
        filepref,filesuff = 'Lobe_','_med1min.txt.gz'
    else:
        assert 2 < 0,"Must provide 'all', 'polarcap', or 'lobe' as location keyword!"
        
    datadir = basedir+'database/Cluster/'+location.lower()+'/'

    clusterfiles = ['Sum_'+filepref+'C1234_2001-2020'+filesuff]

    outparquets = _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,**kws)

    return outparquets


def get_cluster_filenames(filepref,filesuff,**kws):

    years = kws['years']
    sats = kws['sat']
    filetypes = kws['filetype']

    if (not hasattr(years,'__len__')) or isinstance(years,str):
        years = [years]
    
    assert all([isinstance(year,int) for year in years]), \
        "Must provide year as four-digit int between 2001 and 2020!"
        
    if isinstance(sats,str):
        sats = [sats]

    if isinstance(filetypes,str):
        filetypes = [filetypes]

    assert all([sat in ['C1','C2','C3','C4'] for sat in sats])
    assert all([year in range(2001,2021) for year in years])
    assert all([filetype in ['A','S'] for filetype in filetypes])

    # OLD
    # clusterfiles = [f"{filepref:s}{sat:s}_{year:d}_{filetype}{filesuff:s}" for year in years]
    
    # NY
    clusterfiles = []
    for sat in sats:
        for year in years:
            clusterfiles += [f"{filepref:s}{sat:s}_{year:d}_{filetype}{filesuff:s}" for filetype in filetypes]

    return clusterfiles


def _download_process_save_cluster_file(
        datadir,remotedir,clusterfiles,
        **kws):
    
    only_local = kws['only_local'] if 'only_local' in kws else False
    remove_txtgz_files = kws['remove_txtgz_files'] if 'remove_txtgz_files' in kws else True
    overwrite_existing = kws['overwrite_existing'] if 'overwrite_existing' in kws else False
    verbose = kws['verbose'] if 'verbose' in kws else False

    if not only_local:
        user = getpass('Username: ')
        pw = getpass('PW: ')

    outparquets = []
    for clusterfile in clusterfiles:
        url = remotedir+clusterfile
        filename = os.path.basename(urllib.parse.urlparse(url)[2])
        parquetout = filename.replace('.txt','.parquet').replace('.gz','').replace('.xz','')

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
            print(f"Zing! Couldn't get {clusterfile}!")
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

        for clusterfile in clusterfiles:
            url = remotedir+clusterfile
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
    elif '.xz' in clusterfile[-4:]:
        import lzma
        openfunc = lambda f: lzma.open(f, mode='rt', encoding='utf-8')
        if verbose:
            print("Opening with lzma ...")
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
                    assert len(namedict) == 69

            linecount += 1
    
    ####################
    # OLD
    # 
    # datadict = {'Time':[datetime.strptime(dataz[i][0],dtfmt) for i in range(len(dataz))]}    

    ####################
    # NY (fordi me hadde problemer med spørsmåltegn i noen tidsstempler)

    tids = []
    jerks = []
    for i in range(len(dataz)):
        try:
            tids.append(datetime.strptime(dataz[i][0],dtfmt))
        except:
            try:
                tids.append(datetime.strptime(dataz[i][0][:-1],dtfmt))
            except:
                jerks.append(i)
                breakpoint()

    datadict = {'Time':tids}
    
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


def get_bjoern_pc_inds(df):
    # From the documentation:
    # (fys-server-web.uio.no/plasma/cluster/polarCapLobeDataBase/documents/PolarCap_LobeDataBase_version6.pdf)
    # • -5.0 <= Xgsm <= 5.0 [Re]
    # • [-10.0 , -5.0] <= Ygsm <= [10.0 , 5.0] [Re]
    # • [-15.0 , -10.0] <= Zgsm <= [15.0 , 10.0] [Re] and |Zgsm| >= 1.0 [Re]
    # 'X(GSM,GSE)', 'Y(GSM)', 'Z(GSM)'


    # Bjørn's polar cap box
    x1right,y1right = -5,10
    x2right,y2right = 5,5
    mright = (y2right-y1right)/(x2right-x1right)
    
    x1left,y1left = -5,-10
    x2left,y2left = 5,-5
    mleft = (y2left-y1left)/(x2left-x1left)
    
    x1up,z1up = -5,15
    x2up,z2up = 5,10
    mup = (z2up-z1up)/(x2up-x1up)
    
    get_ygsmright = lambda xgsm: mright * (xgsm - x1right) + y1right
    get_ygsmleft = lambda xgsm: mleft * (xgsm - x1left) + y1left
    get_zgsmupper = lambda xgsm: mup * (xgsm - x1up) + z1up
    # print(get_ygsmleft(np.arange(-5,5.1,1)),get_ygsmright(np.arange(-5,5.1,1)))

    inds = (-5 <= df['X(GSM,GSE)']) & (df['X(GSM,GSE)'] <= 5)
    inds = inds & (1 <= np.abs(df['Z(GSM)'])) \
        & (np.abs(df['Z(GSM)']) <= get_zgsmupper(df['X(GSM,GSE)']))

    inds = inds & (get_ygsmleft(df['X(GSM,GSE)']) <= df['Y(GSM)']) \
        & (df['Y(GSM)'] <= get_ygsmright(df['X(GSM,GSE)']))

    indsN = inds & (df['Z(GSM)'] > 0)
    indsS = inds & (df['Z(GSM)'] < 0)
    return indsN,indsS


def get_haaland_2017_fig3a_inds(clusterdf):

    # from journal__20200207__pysymmetry_tsyganenko__polar_cap_fieldlines.ipynb

    gsecols = ['X(GSM,GSE)','Y(GSE)','Z(GSE)']
    # Haaland et al. Figure 3a box
    xL, x0, xR = -13.33333333, 10, 6.66666667
    zU, z0, zL = 6.66666667, 5, 3.333333333 # ACTUAL BOX
    yL, y0, yR = -3.333333333, 0, 3.33333333

    # Get inds in this box
    xinds = (clusterdf[gsecols[0]] >= xL) & (clusterdf[gsecols[0]] <= xR)
    yinds = (clusterdf[gsecols[1]] >= yL) & (clusterdf[gsecols[1]] <= yR)
    zindsN = (clusterdf[gsecols[2]] >= zL) & (clusterdf[gsecols[2]] <= zU)
    zindsS = (clusterdf[gsecols[2]] >= (-1)*zU) & (clusterdf[gsecols[2]] <= (-1)*zL)
    indsN = xinds & yinds & zindsN
    indsS = xinds & yinds & zindsS

    return indsN,indsS


# Apply Stein criteria from Stein's python script
# (see ~/Desktop/Spence_paper_drafts/2020/Cluster_Swarm/Steins_Ne_vs_xz.py)
def get_steinscreen_inds(df,
                         do_NeCond=False,
                         do_VpotCond=True,
                         do_BxCond=True,
                         do_betaCond=True):
    # SQL = "SELECT Cl_X as x,Cl_Y as y, Cl_Z as z,"\
    # " NeEFW as Ne, sqrt(Bx_Clu*Bx_Clu+By_Clu*By_Clu+Bz_Clu*Bz_Clu)"\
    # " as B FROM Lybekk \
    # WHERE Vpot < -10 \
    # AND abs(Bx_Clu) > 5\
    # AND abs(Beta) < 0.05\
    # AND NeEFW < 0.5\
    # AND Cl_X > -12 AND CL_x < -6"

    steinconds = np.isfinite(df['Vpot(EFW)'])

    if do_VpotCond:
        print("Vpot(EFW)   <= -10")
        cond_vpot = df['Vpot(EFW)'] <= -10
        steinconds &= cond_vpot

    if do_BxCond:
        print("Bx(FGM)     >= 5")
        cond_bx = np.abs(df['Bx(FGM)']) >= 5
        steinconds &= cond_bx

    if do_betaCond:
        print("plasma_beta <= 0.05")
        cond_beta = np.abs(df['plasma_beta']) <= 0.05
        steinconds &= cond_beta

    if do_NeCond:
        print("Ne(EFW)     <= 0.5")
        cond_NeEFW = df['Ne(EFW)'] <= 0.5
        steinconds &= cond_NeEFW

    return steinconds
