# 2019/09/10
########################################
# DATAFRAMES
# DF0     (20190910): Apply quantile as agg function
########################################

########################################
# DATAFRAMES
########################################
# DF0     (20190910): Apply quantile as agg function

dfCTcull['qdlat'].agg({'min': np.min,
                       'max': np.max,
                       'q50': lambda x: x.quantile(0.5)
                       'q80': lambda x: x.quantile(0.8)})
