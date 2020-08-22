import pandas as pd
import numpy as np
def load_dst_db(verbose=False):
    indir = '/SPENCEdata/Research/database/dst/'
    infile = 'Dst__2000-202006.parquet'

    print(f"Loading {infile} ...")
    return pd.read_parquet(indir+infile)


def brett_storm_identifier(verbose=False):

    dfdst = load_dst_db()
    
    nDsts = dfdst.shape[0]
    
    storms = np.zeros(nDsts,dtype=np.int64)-1
    
    dst_drop_HH2   = np.zeros(nDsts,dtype=np.float64)-9999.
    dst_drop_HH2L  = np.zeros(nDsts,dtype=np.float64)-9999.
    
    HH1 = 16                        #Defines range of time in hours for: DST a minimum in +/- HH1 hours
    HH2 = 12                        #Defines range of time in hours for: DST drop in previous HH2 hours (NOTE: HH2 must be = or < HH1)
    setHH2L = 16
    HH2L = setHH2L                  #Defines range of time in hours for: DST drop in previous HH2L hours (NOTE: HH2L must be = or < HH1)
    DD1 = 27                        #Defines the DST drop required for SMALL storms
    setDD2 = 55
    DD2 = setDD2                    #Defines the DST drop required for LARGE storms
    
    # The Brett
    
    data_dst = dfdst['Dst'].values
    for ii in np.arange(HH1,nDsts-HH1+1,dtype=np.int64):
        temp = 1
    
        curdst = data_dst[ii]
    
        if verbose:
            if (ii % 1000) == 0:
                print("ii: ",ii,"/",nDsts)
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 2		Indicates DST value is between -50 and -20 (inclusive)
    
        if (curdst <= -20) and (curdst >= -50):
            temp *= 2
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 3		Indicates DST value is less than -50
        if (curdst < -50):
            temp *= 3
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 5		Indicates DST value is a minumum in a time period of plus/minus "HH1" hours
        if (curdst == np.min(data_dst[ii-HH1:ii+HH1+1])):
            temp *= 5
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 7		Indicates DST value dropped at least DD1 nT in the previous "HH2" hours
        if (np.max(data_dst[ii-HH2:ii])-curdst) >= DD1:
            temp *= 7
    
        dst_drop_HH2[ii]  = np.max(data_dst[ii-HH2:ii+1]) - data_dst[ii]
        dst_drop_HH2L[ii] = np.max(data_dst[ii-HH2L:ii+1]) - data_dst[ii]
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 11		Indicates DST value is NOT a repeat of a DST value in the previous 12 hours
        if ( (data_dst[ii] != data_dst[ii-1]) and \
             (data_dst[ii] != data_dst[ii-2]) and \
             (data_dst[ii] != data_dst[ii-3]) and \
             (data_dst[ii] != data_dst[ii-4]) and \
             (data_dst[ii] != data_dst[ii-5]) and \
             (data_dst[ii] != data_dst[ii-6]) and \
             (data_dst[ii] != data_dst[ii-7]) and \
             (data_dst[ii] != data_dst[ii-8]) and \
             (data_dst[ii] != data_dst[ii-9]) and \
             (data_dst[ii] != data_dst[ii-10]) and \
             (data_dst[ii] != data_dst[ii-11]) and \
             (data_dst[ii] != data_dst[ii-12]) ):
            temp *= 11
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 13		Indicates DST value dropped at least DD2 nT in the previous "HH2L" hours
        if (np.max(data_dst[ii-HH2L:ii+1])-data_dst[ii] >= DD2):
            temp *= 13 #should be GE DD2
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 17		Indicates DST value is NOT a repeat of a DST value in the previous 12 hours
    
        if ( (data_dst[ii] != data_dst[ii-1]) and \
             (data_dst[ii] != data_dst[ii-2]) and \
             (data_dst[ii] != data_dst[ii-3]) and \
             (data_dst[ii] != data_dst[ii-4]) and \
             (data_dst[ii] != data_dst[ii-5]) and \
             (data_dst[ii] != data_dst[ii-6]) and \
             (data_dst[ii] != data_dst[ii-7]) and \
             (data_dst[ii] != data_dst[ii-8]) and \
             (data_dst[ii] != data_dst[ii-9]) and \
             (data_dst[ii] != data_dst[ii-10]) and \
             (data_dst[ii] != data_dst[ii-11]) and \
             (data_dst[ii] != data_dst[ii-12]) ):
            temp *= 17
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        #   ;;Prime Factor: 19		Indicates the time is between the specified startdate/enddate
        #   ;;							(required since I load data from/through 7 days before/after the specified startdate/enddate)
        #      IF ( (data_jd[ii] GE jd_stadate) AND $
        #           (data_jd[ii] LE jd_enddate) ) THEN temp=temp*19
        # if ( (data_jd[ii] >= jd_stadate) and \
        #      (data_jd[ii] <= jd_enddate) ):
        temp *= 19
    
        #   ;;-----------------------------------------------------------------------------------------------------------------------------
        storms[ii] = temp		#store the temp index in the actual "storms" array
    
    
    # "STorm TImes SMall" = STTI_SM (Brett's nomenclature)
    stti_sm = ((storms %   2) == 0) & \
        ((storms %   5) == 0) & \
        ((storms %   7) == 0) & \
        ((storms %  11) == 0) & \
        ((storms %  19) == 0)
    # "STorm TImes LarGe" = STTI_LG (Brett's nomenclature)
    stti_lg = ((storms %   3) == 0) & \
        ((storms %   5) == 0) & \
        ((storms %  13) == 0) & \
        ((storms %  17) == 0) & \
        ((storms %  19) == 0)

    return stti_sm,stti_lg
