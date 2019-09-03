# 2019/09/03
# def tools():

import copy
import numpy as np
from scipy.optimize import curve_fit


def func_powerlaw2(x, m, c):
    return x**m * c


def f_resampler(freqs, fResampleStdDevs,
                minFreq=0,
                maxFreq=4):
    f_resampled = np.sort(np.random.normal(
        loc=freqs,
        scale=freqs*fResampleStdDevs))

    badFreqs = (f_resampled > maxFreq) | (
        f_resampled < minFreq)
    if np.where(badFreqs)[0].size > 0:
        f_resampled = f_resampled[~badFreqs]

    return f_resampled


class CalcMedianSpectrogram(object):

    def __init__(self,
                 saveDir='/home/spencerh/Research/sandbox_and_journals/saves_output_etc/',
                 fName=None,
                 # Quant1Dict,
                 # orbsToUse,
                 frequencyBins=None,
                 doRescaleFreqs=None):

        if fName is None:
            print("Provide filename!")
            return None

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

        self.saveDir = saveDir
        self.fName = fName

        self.orbSetNav = orbSetNav

        if frequencyBins is None:
            self.f = np.arange(0, 4.01, .01)
        else:
            self.f = frequencyBins

        if doRescaleFreqs is not None:
            self.f = np.arange(0, 2.01, .01)

        self.mednMaxN = 1000
        self.orbitIndex = np.zeros(self.mednMaxN, dtype=np.int16)

        #
        self.list_spectrogramDict = copy.deepcopy(QuantDict)
        # self.orbits = copy.deepcopy(orbsToUse)

        self.doRescaleFreqs = doRescaleFreqs

        self.nOrbits = 0
        self.origSpecs = {}

        # Laste inn alle origs
        for dude in orbsToUse:

            tmpspec = copy.deepcopy(self.list_spectrogramDict[dude])

            if tmpspec is None:
                print("Orbit {0}: Nannies!".format(dude))
                continue

            # tmpspec[0]: frequencies
            # tmpspec[1]: PSD

            self.origSpecs[dude] = {'f': tmpspec[0].copy(),
                                    'spec': tmpspec[1].copy()}

            self.orbitIndex[self.nOrbits] = dude

            self.nOrbits += 1

        self.orbitIndex = self.orbitIndex[0:self.nOrbits]

        self.calcMedian(self.f)

    def calcMedian(self, interpfrequencies):

        # Initialize interpedSpecArr
        self.interpedSpecArr = np.zeros((interpfrequencies.size,
                                         self.nOrbits))

        for iOrb, dude in enumerate(self.orbitIndex):

            # tmpspec = copy.deepcopy(self.list_spectrogramDict[dude])

            # if tmpspec is None:
            #     print("Orbit {0}: Nannies!".format(dude))
            #     continue

            # tmpspec[0]: frequencies
            # tmpspec[1]: PSD

            # tmpFreqs = tmpspec[0].copy()

            # Krever refSpeed, speedDict, useSpeedInd
            # if self.doRescaleFreqs is not None:
            #     tmpFreqs = tmpFreqs*refSpeed / \
            #         np.abs(speedDict[dude][useSpeedInd])

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
        self.mednSpec = np.median(self.interpedSpecArr, axis=1)

    def init_MonteCarlo(self, MonteCarloParams=None,
                        minFreq=None,
                        maxFreq=None):

        # MonteCarloParams = {'slope':(-4,-1),
        #                     'coeff':(0.5,10),
        #                     'N':1000,
        #                     'fResampleStdDevs':MCresampleXStdDevs}

        if MonteCarloParams is None:
            print(
                "Please provide dict with 'slope','coeff','N', and 'fResampleStdDevs' keys")
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

        do_f_resampling = 'fResampleStdDevs' in self.MonteCarloParams
        if do_f_resampling:
            # if self.fResampler is None:
            print("Init fResampler ...")

            self.fResampler = lambda: f_resampler(self.origFreq, self.MonteCarloParams['fResampleStdDevs'],
                                                  minFreq=np.min(
                                                      self.origFreq),
                                                  maxFreq=np.max(self.origFreq))

    def execute_MonteCarlo(self,
                           verbose=True,
                           use_fixed_slope=False,
                           fixed_slope_val=-2,
                           maxfev=2000):

        self.use_fixed_slope = use_fixed_slope
        self.fixed_slope_val = fixed_slope_val

        # Init all the neededses
        self.MC_medianspecs = []
        self.MC_medianfreqs = []

        if self.use_fixed_slope:
            self.fitFunc = lambda x, c: func_powerlaw2(
                x, self.fixed_slope_val, c)
        else:
            self.fitFunc = func_powerlaw2

        if verbose:
            print("Regne ut median-spectra ...")

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

            # breakpoint()

            self.MC_medianfreqs.append(copy.deepcopy(self.fResampler()))
            self.calcMedian(self.MC_medianfreqs[i])
            self.MC_medianspecs.append(copy.deepcopy(self.mednSpec))

            # NY

            slopeTmp = np.random.uniform(self.MonteCarloParams['slope'][0],
                                         self.MonteCarloParams['slope'][1])
            if self.use_fixed_slope:
                p0 = [slopeTmp]
            else:
                coeffTmp = np.random.uniform(self.MonteCarloParams['coeff'][0],
                                             self.MonteCarloParams['coeff'][1])
                p0 = [slopeTmp, coeffTmp]
            # print(i, p0)

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
                    fixed_slope_val=-2):  # ,
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

        if self.use_fixed_slope:
            self.fitFunc = lambda x, c: func_powerlaw2(
                x, self.fixed_slope_val, c)
        else:
            self.fitFunc = func_powerlaw2

        if self.SingleMonteCarloParams is None:
            sol1 = curve_fit(self.fitFunc, x, y, maxfev=maxfev)
            # sol2 = curve_fit(func_powerlaw, x, y, p0=p0))

            self.fitParams = sol1[0]
            self.fitCov = sol1[1]
            self.fitX = x
            self.fitY = self.fitFunc(self.fitX, *self.fitParams)
            self.fitParmErrs = np.sqrt(np.diag(sol1[1]))

        elif self.SingleMonteCarloParams is not None:

            do_f_resampling = 'fResampleStdDevs' in self.SingleMonteCarloParams
            if do_f_resampling:
                print("Also resampling xs!")
                self.fResampler = lambda: (np.random.normal(
                    loc=self.f,
                    scale=self.f*self.SingleMonteCarloParams['fResampleStdDevs']))

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

                slopeTmp = np.random.uniform(self.SingleMonteCarloParams['slope'][0],
                                             self.SingleMonteCarloParams['slope'][1])
                if self.use_fixed_slope:
                    p0 = [slopeTmp]
                else:
                    coeffTmp = np.random.uniform(self.SingleMonteCarloParams['coeff'][0],
                                                 self.SingleMonteCarloParams['coeff'][1])
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
