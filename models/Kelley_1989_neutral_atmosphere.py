# 2019/07/31
# def Kelley_1989_neutral_atmosphere():

import numpy as np
import pandas as pd


def load_sunspot_max_model():
    outdir = '/SPENCEdata/Research/database/'
    outfile = 'Kelley_1989__TableB1_neutral_atmosph.csv'
    return pd.read_csv(outdir+outfile)
