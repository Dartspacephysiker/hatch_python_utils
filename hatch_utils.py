# 2017/02/22
import numpy as np
import scipy.io as sio
# import ipdb

def get_basedir(return_bool=True,verbose=True):

    import os
    thisIsColtrane = os.path.isdir("/Data/ift/ift_romfys1/Q1/folk/spencer/")
    if thisIsColtrane:
        if verbose:
            print("This is Coltrane!")
        baseDir = '/Data/ift/ift_romfys1/Q1/folk/spencer/Coltrane/home/Research/'
    else:
        if verbose:
            print("This ain't Coltrane!")
        baseDir = '/SPENCEdata/Research/'

    if return_bool:
        return baseDir, thisIsColtrane
    else:
        return baseDir

def print_dict(d):
    if isinstance(d, dict):
        print('\n')
        for key in d:
            print('******' + key.upper() + '******')
            print('{!s:15} : {!s:25}, {!s:<10}'.format(
                key, type(d[key]), len(d[key])))
            if type(d[key]) is np.ndarray:
                print('ndim : {!s:>10}'.format(d[key].ndim))
                print('shape: {!s:>10}'.format(d[key].shape))
                print('dtype: {!s:>10}'.format(d[key].dtype))
                print('flags: ')
                print(d[key].flags)
            print('')
    else:
        print("Not a dictionary")
        return


def read_mat_file(matFile='', datDir=''):
    if not (isinstance(matFile, str) and isinstance(matFile, str)):
        print("That's totally illegal.")
        return
    temp_list = [matFile]
    temp_list = [datDir + s for s in temp_list]
    matdat = sio.loadmat(temp_list[0])

    print("Loaded {}".format(matFile))

    print('matdat has : ', matdat.keys())

    print_dict(matdat)
    return matdat


def _get_args_dict(fn, args, kwargs):
    args_names = fn.__code__.co_varnames[:fn.__code__.co_argcount]
    return dict(zip(args_names, args))
