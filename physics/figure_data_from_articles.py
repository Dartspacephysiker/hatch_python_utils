# 2019/08/27
# def figure_data_from_articles():


# (20190827) Johnson 1961, "Satellite Data Handbook", Ion collision frequencies
# (20190827) Lysak 1997, "", Alfvén speed profile (Figure 1)
# (20190827) Lotko and Zhang 2018, "", Alfvén speed, sigma0, sigmaP, etaPar profiles (Figure 1)


########################################
# (20190827) Johnson 1961, "Satellite Data Handbook", Ion collision frequencies
# Used in journal__20190813__Lotko_Zhang__survey_FORTRAN.ipynb

Johnsondir = '/SPENCEdata/Research/database/From_articles/Johnson - 1961/'

ion_collfreq_daysunmax = 'Johnson_1961__day_sunspot_max_ion_collision_freq.csv'
ion_collfreq_daysunmin = 'Johnson_1961__day_sunspot_min_ion_collision_freq.csv'

daysunmax = pd.read_csv(Johnsondir+ion_collfreq_daysunmax,
                        header=None, names=['nu', 'alt'])
daysunmin = pd.read_csv(Johnsondir+ion_collfreq_daysunmin,
                        header=None, names=['nu', 'alt'])

daysunmax.set_index('alt', inplace=True)
daysunmin.set_index('alt', inplace=True)

########################################
# (20190827) Lysak 1997, "", Alfvén speed profile (Figure 1)
# Used in journal__20190813__Lotko_Zhang__survey_FORTRAN.ipynb

Lysak1997Fig1AlfSpeedFile = "/SPENCEdata/Research/database/From_articles/Lysak - 1997/Lysak19997_alfSpeed.csv"
LysAlf = pd.read_csv(Lysak1997Fig1AlfSpeedFile, header=None)
LysAlf.sort_values(by=[0], inplace=True)
LysAlf.loc[:, 0] = LysAlf[0].rolling(5, center=True, min_periods=1).mean()
junk = axes[0].plot(LysAlf[1].values,
                    LysAlf[0].values,
                    alpha=0.5,
                    label='Lysak 1997, Fig 1')

########################################
# (20190827) Lotko and Zhang 2018, "", Alfvén speed, sigma0, sigmaP, etaPar profiles (Figure 1)
# Used in journal__20190813__Lotko_Zhang__survey_FORTRAN.ipynb
Lotko2018Fig1aAlfSpeedFile = "/SPENCEdata/Research/database/From_articles/Lotko and Zhang - 2018/LZ_Vprof.csv"
LotkAlf = pd.read_csv(Lotko2018Fig1aAlfSpeedFile, header=None)
LotkAlf.sort_values(by=[1], inplace=True)
# LotkAlf.loc[:,0] = LotkAlf[0].rolling(5,center=True, min_periods=1).mean()
junk = axes[0].plot(LotkAlf[0].values,
                    LotkAlf[1].values,
                    alpha=0.5,
                    label='Lotko and Zhang 2018, Fig 1')
