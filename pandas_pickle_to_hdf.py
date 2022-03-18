# ;; This buffer is for text that is not saved, and for Lisp evaluation.
# ;; To create a file, visit it with C-x C-f and enter text in its buffer.

import pandas as pd
import glob
import os
def pandas_pickle_to_hdf(fil,remove=False):
    # print(f"Opening {fil}")
    new = fil.replace('.pkl','.h5')
    try:
        df = pd.read_pickle(fil)
        print(f"Writing to {new}",end='...')
        df.to_hdf(new,'data')
        if remove:
            print("Removing orig pickle …")
            os.remove(fil)
        print('Done!')
    except Exception as e:
        print(e)
        print("Skipping!!")

if __name__ == "__main__":
    files = glob.glob('*.pkl')
    print("Found these files:")
    for i,f in enumerate(files):
        print(f"{i:3d}: {f}")
    print()
    cont = False
    while not cont:
        answer = input("Should I convert them to .h5? (y/n)")
        if answer.lower().startswith('y'):
            cont = True
        elif answer.lower().startswith('n'):
            print("OK, exiting ...")
            exit
        else:
            print("What? Do y/n")

    cont = False
    remove = False
    while not cont:
        answer = input("Should I remove the original pickles after successful .h5 conversion? (y/n)")
        if answer.lower().startswith('y'):
            cont = True
            remove = True
        elif answer.lower().startswith('n'):
            remove = False
            print("OK, not removing ...")
            exit
        else:
            print("What? Do y/n")


    for i,f in enumerate(files):
        print(f"{i:3d}: {f}")
        pandas_pickle_to_hdf(f,remove=remove)
