# 2018/03/13
import os
import idlsave
from pathlib import Path


def load_omniDB():
    dir = '/SPENCEdata/Research/database/OMNI/'
    file = 'culled_OMNI_magdata.dat'
    # omnifile = Path(dir+file)

    if os.path.exists(dir+file):
        print("Opening " + file + ' ...')
        omni = idlsave.read(dir+file)
        return omni
    else:
        # doesn't exist
        print("Couldn't get OMNI IDL save file!!! Returning ...")
        return

    # try:
    #     omni_path = omnifile.resolve():
    # except FileNotFoundError:
    #     # doesn't exist
    #     print("Couldn't get OMNI IDL save file!!! Returning ...")
    #     return
    # else:
    #     # exists
    #     print("Opening " + file + ' ...')
    #     omni = idlsave.read(dir+file)
    #     return omni
