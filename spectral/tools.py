# 2019/09/03
# def tools():

import copy
import numpy as np
from scipy.optimize import curve_fit
import random

import matplotlib.pyplot as plt


def func_powerlaw2(x, m, c):
    return x**m * c


def func_loggedpowerlaw2(log10x, m, log10c):
    return 10**(m * log10x + log10c)


def loguniform(low=1e-10, high=1, size=None):
    return np.exp(np.random.uniform(np.log(low), np.log(high), size))


def f_resampler(freqs, fResampleNormedStdDevs,
                minFreq=0,
                maxFreq=4):
    f_resampled = np.sort(np.random.normal(
        loc=freqs,
        scale=freqs*fResampleNormedStdDevs))

    badFreqs = (f_resampled > maxFreq) | (
        f_resampled < minFreq)
    if np.where(badFreqs)[0].size > 0:
        f_resampled = f_resampled[~badFreqs]

    return f_resampled


class CalcMedianSpectrogram(object):

    def __init__(self,
                 saveDir='/home/spencerh/Research/sandbox_and_journals/saves_output_etc/',
                 fName=None,
                 frequencyBins=None,
                 doRescaleFreqs=None):

        if fName is None:
            print("Provide filename!")
            return None

        # self.saveDir = saveDir

        self.fName = None
        self.orbSetNav = None

        self.nOrbits = 0
        self.origSpecs = {}

        self.MonteCarloParams = None

        self.MCFitX = None
        self.MCFitY = None

        if frequencyBins is None:
            self.f = np.arange(0, 4.01, .01)
        else:
            self.f = frequencyBins

        if doRescaleFreqs is not None:
            self.f = np.arange(0, 2.01, .01)

        self.mednMaxN = 1000
        self.orbitIndex = np.zeros(self.mednMaxN, dtype=np.int16)

        self.doRescaleFreqs = doRescaleFreqs

        self.list_spectrogramDict = {}

        if isinstance(fName, list):
            for fNm in fName:
                self.add_specs_from_file(saveDir, fNm)
        else:
            self.add_specs_from_file(saveDir, fName)

        self.calcMedian(self.f)

    def add_specs_from_file(self, saveDir, fName):

        # try:
        npz = np.load(saveDir+fName, allow_pickle=True)
        # except:
        #     print("Couldn't load {:s}! Returning ...".format(fName))
        #     return None

        # if useQuant2Dict:
        # useDict = 'Quant2Dict'
        # else:
        useDict = 'Quant1Dict'

        QuantDict = npz[useDict].item()
        orbsToUse = npz['orbsToUse']
        orbSetNav = str(npz['orbSetNav'])
        refSpeed = npz['refSpeed']
        speedDict = npz['speedDict']
        useSpeedInd = npz['useSpeedInd']

        if self.fName == None:
            self.fName = [saveDir+fName]
        else:
            self.fName.append(saveDir+fName)

        if self.orbSetNav == None:
            self.orbSetNav = [orbSetNav]
        else:
            self.orbSetNav.append(orbSetNav)

        #
        self.list_spectrogramDict = {**self.list_spectrogramDict, **QuantDict}
        # self.orbits = copy.deepcopy(orbsToUse)

        foreOrbs = self.nOrbits

        # Laste inn alle origs
        for dude in orbsToUse:

            tmpspec = QuantDict[dude]

            if tmpspec is None:
                print("Orbit {0}: Nannies!".format(dude))
                continue

            # tmpspec[0]: frequencies
            # tmpspec[1]: PSD

            self.origSpecs[dude] = {'f': tmpspec[0].copy(),
                                    'spec': tmpspec[1].copy()}

            self.orbitIndex[self.nOrbits] = dude

            self.nOrbits += 1

        print("Added {:d} orbits from {:s} ...".format(
            self.nOrbits-foreOrbs, fName))

        # self.orbitIndex = self.orbitIndex[0:self.nOrbits]

    def calcMedian(self, interpfrequencies,
                   isMonteCarlo=False):

        # Initialize interpedSpecArr
        self.interpedSpecArr = np.zeros((interpfrequencies.size,
                                         self.nOrbits))

        for iOrb, dude in enumerate(self.orbitIndex[0:self.nOrbits]):

            # x : The x-coordinates at which to evaluate the interpolated values.
            # xp : The x-coordinates of the data points, must be increasing if argument period is not specified.
            #      Otherwise, xp is internally sorted after normalizing the periodic boundaries with xp = xp % period.
            # fp : The y-coordinates of the data points, same length as xp.
            # left : Value to return for x < xp[0], default is fp[0].
            # right : Value to return for x > xp[-1], default is fp[-1].
            self.interpedSpecArr[:, iOrb] = np.interp(
                interpfrequencies,
                self.origSpecs[dude]['f'],
                self.origSpecs[dude]['spec'],
                left=np.nan, right=np.nan, period=None)

        # Now median
        if isMonteCarlo:
            self.tmpMCmednSpec = np.median(self.interpedSpecArr, axis=1)
        else:
            self.mednSpec = np.median(self.interpedSpecArr, axis=1)

    def init_MonteCarlo(self, MonteCarloParams=None,
                        minFreq=None,
                        maxFreq=None,
                        do_logfits=True,
                        use_loguniform_for_coeff0=True):

        # MonteCarloParams = {'slope':(-4,-1),
        #                     'coeff':(0.5,10),
        #                     'N':1000,
        #                     'fResampleNormedStdDevs':MCresampleXStdDevs}

        if MonteCarloParams is None:
            print(
                "Please provide dict with 'slope','coeff','N', and 'fResampleNormedStdDevs' keys")
            return

        self.origFreq = copy.deepcopy(self.f)

        self.MC_minFreq = minFreq
        self.MC_maxFreq = maxFreq

        if self.MC_minFreq is not None:
            binStart_i = np.where(self.origFreq >= self.MC_minFreq)[0][0]
        else:
            binStart_i = 0

        if self.MC_maxFreq is not None:
            binStop_i = np.where(self.origFreq <= self.MC_maxFreq)[0][-1]
        else:
            binStop_i = self.origFreq.size-1

        self.origFitIndices = np.arange(binStart_i, binStop_i+1)

        self.MonteCarloParams = copy.deepcopy(MonteCarloParams)

        do_f_resampling = 'fResampleNormedStdDevs' in self.MonteCarloParams
        if do_f_resampling:
            # if self.fResampler is None:
            print("Init fResampler ...")

            self.fResampler = lambda: f_resampler(self.origFreq, self.MonteCarloParams['fResampleNormedStdDevs'],
                                                  minFreq=np.min(
                                                      self.origFreq),
                                                  maxFreq=np.max(self.origFreq))

        self.do_logfits = False
        if do_logfits:
            print("Doing fits to log of all quants")

            # self.baseFitFunc = func_loggedpowerlaw2
            self.baseFitFunc = lambda x, m, c: func_loggedpowerlaw2(
                np.log10(x), m, np.log10(c))

            self.do_logfits = True
        else:
            print("Doing fits to a straight powerlaw...")
            self.baseFitFunc = func_powerlaw2

        if use_loguniform_for_coeff0:
            print("Coeff inits med log-uniform!")
            self.coeff0func = loguniform
            self.use_loguniform_for_coeff0 = True
        else:
            print("Coeff inits med (not log) uniform!")
            self.coeff0func = np.random.uniform
            self.use_loguniform_for_coeff0 = False

    def execute_MonteCarlo(self,
                           verbose=True,
                           use_fixed_slope=False,
                           fixed_slope_val=-2,
                           maxfev=2000):
        """
        Use fResampler to randomize selection of resampling frequencies, regne ut medianspektrum; do it MonteCarloParams['N'] times. Then
        perform fit to each medianspektrum to get dists of fit params.
        """
        self.use_fixed_slope = use_fixed_slope
        self.fixed_slope_val = fixed_slope_val

        # Init all the neededses
        self.MC_medianspecs = []
        self.MC_medianfreqs = []
        self.MC_p0 = []

        if verbose:
            print("Regne ut median-spectra ...")
            print(
                "Monte Carlo-picking frequencies for resampling all spectra, and performing  ...")

        if self.use_fixed_slope:
            self.fitFunc = lambda x, c: self.baseFitFunc(
                x, self.fixed_slope_val, c)
        else:
            self.fitFunc = self.baseFitFunc

        MCFitParms = []
        MCFitCov = []

        # minFitFreq = self.MC_minFreq
        # maxFitFreq = self.MC_maxFreq
        # boundTypeStr = "MC min/max freqs"

        minFitFreq = np.min(self.origFreq)
        maxFitFreq = np.max(self.origFreq)
        boundTypeStr = "orig resampler freqs"
        print("Using {:s} ({:.3f} Hz, {:.3f} Hz) as bounds!".format(boundTypeStr,
                                                                    minFitFreq, maxFitFreq))

        for i in range(self.MonteCarloParams['N']):
            # if do_f_resampling:
            # self.medianSpec = np.interp(self.freqsResample, self.origFreq, self.origSpec,
            #                             left=np.nan, right=np.nan, period=None)

            self.MC_medianfreqs.append(copy.deepcopy(self.fResampler()))
            self.calcMedian(self.MC_medianfreqs[i], isMonteCarlo=True)
            self.MC_medianspecs.append(copy.deepcopy(self.tmpMCmednSpec))

            # NY

            coeffTmp = self.coeff0func(self.MonteCarloParams['coeff'][0],
                                       self.MonteCarloParams['coeff'][1])

            if self.use_fixed_slope:
                p0 = [coeffTmp]
            else:
                slopeTmp = np.random.uniform(self.MonteCarloParams['slope'][0],
                                             self.MonteCarloParams['slope'][1])

                p0 = [slopeTmp, coeffTmp]

            # print(i, p0)

            self.MC_p0.append(p0)

            xAll = self.MC_medianfreqs[i]

            # LIMIT TO FITLIMS
            if self.MC_minFreq is not None:
                binStart_i = np.where(xAll >= self.MC_minFreq)[0][0]
            else:
                binStart_i = 0

            if self.MC_maxFreq is not None:
                binStop_i = np.where(xAll <= self.MC_maxFreq)[0][-1]
            else:
                binStop_i = xAll.size-1

            xIndices = np.arange(binStart_i, binStop_i+1)

            x = xAll[xIndices]
            y = self.MC_medianspecs[i][xIndices]

            if x.size <= 2:
                print("Gonna run into trouble!")
                print("Might need to increase the number of points that we fit to ...")

            if np.where(~(np.isfinite(x) & np.isfinite(y)))[0].size > 0:
                breakpoint()

            sol1 = curve_fit(self.fitFunc, x, y,
                             maxfev=maxfev,
                             p0=p0)

            MCFitParms.append(sol1[0])
            MCFitCov.append(sol1[1])

        MCFitParms = np.array(MCFitParms)
        MCFitCov = np.array(MCFitCov)

        self.MCFitParms = MCFitParms
        self.MCFitCov = MCFitCov

        if len(MCFitParms.shape) > 1:
            self.MCfitParamsFinal = np.median(MCFitParms, axis=0)
            self.MCfitCovFinal = np.median(MCFitCov, axis=0)
        else:
            self.MCfitParamsFinal = np.median(MCFitParms)
            self.MCfitCovFinal = np.median(MCFitCov)

        self.MCFitX = self.origFreq[self.origFitIndices]
        self.MCFitY = self.fitFunc(self.MCFitX,
                                   *self.MCfitParamsFinal)

        self.final_MC_median()

    def final_MC_median(self):

        # Initialize interpedSpecArr
        tmpInterpedSpecArr = np.zeros((self.origFreq[self.origFitIndices].size,
                                       len(self.MC_medianfreqs)))

        for iMC in range(self.MonteCarloParams['N']):

            tmpInterpedSpecArr[:, iMC] = np.interp(
                self.origFreq[self.origFitIndices],
                self.MC_medianfreqs[iMC],
                self.MC_medianspecs[iMC],
                left=np.nan, right=np.nan, period=None)

        # Now median
        self.MCmednSpec = np.median(tmpInterpedSpecArr, axis=1)

    def powerLawFit(self, minFreq=None, maxFreq=None,
                    maxfev=2000, p0=None,
                    use_fixed_slope=False,
                    fixed_slope_val=-2,
                    do_logfits=True):  # ,
                    # MonteCarloParams=None):

        self.SingleMonteCarloParams = None

        self.use_fixed_slope = False
        self.fixed_slope_val = fixed_slope_val

        if minFreq is not None:
            binStart_i = np.where(self.f >= minFreq)[0][0]
        else:
            binStart_i = 0

        if maxFreq is not None:
            binStop_i = np.where(self.f <= maxFreq)[0][-1]
        else:
            binStop_i = self.f.size-1

        self.fitIndices = np.arange(binStart_i, binStop_i+1)

        x = self.f[self.fitIndices]
        y = self.mednSpec[self.fitIndices]

        if do_logfits:
            print("Doing fits to log of all quants")
            self.baseFitFunc = lambda x, m, c: func_loggedpowerlaw2(
                np.log10(x), m, np.log10(c))
        else:
            print("Doing fits to a straight powerlaw...")
            self.baseFitFunc = func_powerlaw2

        if self.use_fixed_slope:
            self.fitFunc = lambda x, c: self.baseFitFunc(
                x, self.fixed_slope_val, c)
        else:
            self.fitFunc = self.baseFitFunc

        if self.SingleMonteCarloParams is None:
            sol1 = curve_fit(self.fitFunc, x, y, maxfev=maxfev)
            # sol2 = curve_fit(func_powerlaw, x, y, p0=p0))

            self.fitParams = sol1[0]
            self.fitCov = sol1[1]
            self.fitX = x
            self.fitY = self.fitFunc(self.fitX, *self.fitParams)
            self.fitParmErrs = np.sqrt(np.diag(sol1[1]))

        elif self.SingleMonteCarloParams is not None:

            assert 2 < 0, "What was this for again? I'm using journal__20190903__Chrisfit__powerlawfits_to_MonteCarloed_median_spectra.ipynb, and it totally escapes me"

            if use_loguniform_for_coeff0:
                print("Coeff inits med log-uniform!")
                self.singlecoeff0func = loguniform
            else:
                print("Coeff inits med (not log) uniform!")
                self.singlecoeff0func = np.random.uniform

            do_f_resampling = 'fResampleNormedStdDevs' in self.SingleMonteCarloParams
            if do_f_resampling:
                print("Also resampling xs!")
                self.fResampler = lambda: (np.random.normal(
                    loc=self.f,
                    scale=self.f*self.SingleMonteCarloParams['fResampleNormedStdDevs']))

            MCFitParms = []
            MCFitCov = []
            for i in range(self.SingleMonteCarloParams['N']):
                if do_f_resampling:

                    self.fResample = self.fResampler()

                    self.fResample = np.sort(self.fResample)

                    badFreqs = (self.fResample > np.max(self.origFreq)) | (
                        self.fResample < np.min(self.origFreq))
                    if np.where(badFreqs)[0].size > 0:
                        self.fResample = self.fResample[~badFreqs]

                    self.mednSpec = np.interp(self.fResample, self.origFreq, self.origSpec,
                                              left=np.nan, right=np.nan, period=None)

                    if minFreq is not None:
                        binStart_i = np.where(
                            self.fResample >= minFreq)[0][0]
                    else:
                        binStart_i = 0

                    if maxFreq is not None:
                        binStop_i = np.where(
                            self.fResample <= maxFreq)[0][-1]
                    else:
                        binStop_i = self.fResample.size-1

                    self.tmpfitIndices = np.arange(binStart_i, binStop_i+1)

                    x = self.fResample[self.tmpfitIndices]
                    y = self.mednSpec[self.tmpfitIndices]

                # Now get p0
                coeffTmp = self.singlecoeff0func(self.SingleMonteCarloParams['coeff'][0],
                                                 self.SingleMonteCarloParams['coeff'][1])
                if self.use_fixed_slope:
                    p0 = [coeffTmp]
                else:
                    slopeTmp = np.random.uniform(self.SingleMonteCarloParams['slope'][0],
                                                 self.SingleMonteCarloParams['slope'][1])
                    p0 = [slopeTmp, coeffTmp]
                # print(i, p0)

                if np.where(~(np.isfinite(x) & np.isfinite(y)))[0].size > 0:
                    breakpoint()

                sol1 = curve_fit(self.fitFunc, x, y,
                                 maxfev=maxfev,
                                 p0=p0)

                MCFitParms.append(sol1[0])
                MCFitCov.append(sol1[1])

            MCFitParms = np.array(MCFitParms)
            MCFitCov = np.array(MCFitCov)

            # breakpoint()

            self.MCFitParms = MCFitParms
            self.MCFitCov = MCFitCov

            if len(MCFitParms.shape) > 1:
                self.fitParams = np.median(MCFitParms, axis=0)
                self.fitCov = np.median(MCFitCov, axis=0)
            else:
                self.fitParams = np.median(MCFitParms)
                self.fitCov = np.median(MCFitCov)

            self.fitY = self.fitFunc(self.f[self.fitIndices],
                                     *self.fitParams)

    def plot__show_median_median_spec(self, ax):
        return ax.plot(self.origFreq[self.origFitIndices], self.MCmednSpec,
                       label='Median median spec (N = {:d})'.format(self.MonteCarloParams['N']))

    def plot__show_fit_line(self, ax):

        fitLabel = r'P = $P_0 f^a$' + " (a = {:.2f}, $P_0$ = {:.2f})".format(self.MCfitParamsFinal[0],
                                                                             self.MCfitParamsFinal[1])

        # junk = ax.plot(self.origFreq[self.origFitIndices],self.MCmednSpec,
        #                label='Median median spec (N = {:d})'.format(self.MonteCarloParams['N']))
        junk = ax.plot(self.MCFitX, self.MCFitY,
                       color='orange',
                       label=fitLabel)

        return junk

    def plot_psds(self,
                  ylim=(1e-2, 7e7),
                  show_fit=True,
                  show_median_median_spec=False,
                  indivPSDPlotOpts=None,
                  medianPSDPlotOpts=None):

        if show_fit and self.MCFitX is None:
            print("Can't show fit!")
            # return None, None

        fig, ax = plt.subplots(1, 1)

        ax.grid()

        junk = ax.set_yscale('log')
        junk = ax.set_xscale('log')

        junk = ax.set_ylabel(r'PSD (nT$^2$/Hz)')
        junk = ax.set_xlabel('Frequency (Hz)')

        if indivPSDPlotOpts is None:

            indivPSDPlotOpts = dict(color='black',
                                    alpha=0.1)

        if medianPSDPlotOpts is None:

            medianPSDPlotOpts = dict(color='blue',
                                     alpha=1.0,
                                     linewidth=2.5)

        for iOrb, dude in enumerate(self.orbitIndex[0:self.nOrbits]):

            junk = ax.plot(self.origSpecs[dude]['f'],
                           self.origSpecs[dude]['spec'],
                           **indivPSDPlotOpts)

        junk = ax.plot(self.f,
                       self.mednSpec,
                       **medianPSDPlotOpts,
                       label="Orig median")

        junk = ax.set_ylim(ylim)
        # junk = ax.set_xlim('Frequency (Hz)')

        if show_fit and (self.MCFitX is not None):
            junk = self.plot__show_fit_line(ax)

        if show_median_median_spec:
            junk = self.plot__show_median_median_spec(ax)

        junk = ax.legend()

        return fig, ax

    def plot_MC_medianspecs(self,
                            Nlim=10000,
                            show_fit=True,
                            show_median_median_spec=False,
                            N_random_sample=None,
                            verbose=False,
                            plotOpts=None):

        if self.MonteCarloParams is None:
            print("Haven't yet initialized/executed Monte Carlo fits! Returning ... ")
            return None, None

        if plotOpts == None:
            plotOpts = dict(color='black',
                            alpha=0.01)

        fig, ax = plt.subplots(1, 1)

        ax.grid()

        junk = ax.set_yscale('log')
        junk = ax.set_xscale('log')

        junk = ax.set_ylabel(r'PSD (nT$^2$/Hz)')
        junk = ax.set_xlabel('Frequency (Hz)')

        junk = fig.suptitle(
            "N = {:d} Monte Carlo specs".format(len(self.MC_medianspecs)))

        inds = list(np.arange(len(self.MC_medianfreqs)))
        if N_random_sample is not None:
            if verbose:
                print("{:d}-sample of MC median spectra ...".format(N_random_sample))

            inds = random.sample(inds, k=N_random_sample)

        # for i, (freq, spec) in enumerate(zip(self.MC_medianfreqs, self.MC_medianspecs)):
        for iSpec in inds:

            freq = self.MC_medianfreqs[iSpec]
            spec = self.MC_medianspecs[iSpec]
            junk = ax.plot(freq,
                           spec,
                           **plotOpts)

        if show_fit and (self.MCFitX is not None):
            junk = self.plot__show_fit_line(ax)

        if show_median_median_spec:
            junk = self.plot__show_median_median_spec(ax)

        junk = ax.legend()

        return fig, ax

    def plot_MC_param_histos(self,
                             show_p0=False):
        if self.MonteCarloParams is None:
            print("Haven't yet initialized/executed Monte Carlo fits! Returning ... ")
            return None, None

        titleSuff = ''
        if show_p0:
            if self.use_fixed_slope:
                x0 = np.array([self.fixed_slope_val]*len(self.MC_p0))
                x1 = np.array(self.MC_p0)
            else:
                x0 = np.array([p0[0] for p0 in self.MC_p0])
                x1 = np.array([p0[1] for p0 in self.MC_p0])

            x0final = np.median(x0)
            x1final = np.median(x1)

            titleSuff = ' (p0s)'

        else:
            x0 = self.MCFitParms[:, 0]
            x1 = self.MCFitParms[:, 1]
            x0final = self.MCfitParamsFinal[0]
            x1final = self.MCfitParamsFinal[1]

        P0label = r'$P_0$'
        if self.use_loguniform_for_coeff0:
            x1 = np.log10(x1)
            x1final = np.log10(x1final)
            P0label = r'log($P_0$)'

        fig2, ax2 = plt.subplots(1, 2)

        junk = fig2.suptitle(",".join(self.orbSetNav)+titleSuff)

        junk = ax2[0].hist(x0)
        junk = ax2[0].axvline(x0final,
                              color='red',
                              label='med(a) = {:.3f}'.format(x0final))

        junk = ax2[1].hist(x1)
        junk = ax2[1].axvline(x1final,
                              color='red',
                              label=r'med($P_0$) = {:.3f}'.format(x1final))

        # if self.use_loguniform_for_coeff0:
        #     junk = ax2[1].set_xscale('log')

        junk = ax2[0].set_xlabel('a')
        junk = ax2[1].set_xlabel(P0label)

        junka = ax2[0].legend()
        junkb = ax2[1].legend()

        return fig2, ax2
