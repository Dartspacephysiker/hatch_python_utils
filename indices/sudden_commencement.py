import pandas as pd
import numpy as np

def load_sc_database():
    scdir = '/SPENCEdata/Research/database/sudden_commencement/'
    # output info
    outfile = 'SuddenCommencements__1999-2019.parquet'
    return pd.read_parquet(scdir+outfile)
