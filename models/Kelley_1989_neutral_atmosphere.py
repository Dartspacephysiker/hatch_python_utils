# 2019/07/31

import numpy as np
import pandas as pd


def load_model(*args, units='cm', sunspot_max=True,
               verbose=False):
    """
    load_model([z_km])

    z_km (array_like): Provide if you'd like the model interpolated to a new set of altitudes
    """
    outdir = '/SPENCEdata/Research/database/'
    if sunspot_max:
        print("Loading sunspot max model ...")
        outfile = 'Kelley_1989__TableB1_neutral_atmosph.csv'
    else:
        print("Loading sunspot min model ...")
        outfile = 'Kelley_1989__TableB2_neutral_atmosph.csv'

    kelley = pd.read_csv(outdir+outfile)

    if units == 'm':
        kelley['n'] *= 1.e6

    if len(args) > 0:
        if verbose:
            print("Interpolating to new heights ...")

        z_km = args[0]

        # interpolate

        # Do lin-log interp
        kelley['n'] = np.log10(kelley['n'])

        kelley.set_index('h', inplace=True)
        kelley = pd.concat([kelley, kelley.reindex(z_km)]).sort_index()
        kelley = kelley[~kelley.index.duplicated(keep='first')]

        # kelley = kelley.interpolate(method='polynomial', order=2)
        kelley = kelley.interpolate(method='linear')

        kelley['n'] = 10**(kelley['n'])
        for col in ['T', 'M', 'n']:
            kelley.loc[:, [col]] = kelley.loc[:, [col]].rolling(
                21, center=True, min_periods=1).mean()

        kelley = kelley.reset_index()

    kelley.set_index('h', inplace=True)

    # else:
    return kelley
