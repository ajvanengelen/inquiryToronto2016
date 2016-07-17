

import transitTools
import numpy as np
import matplotlib.pyplot as plt
import pickle
def twodlist(rows, cols):
    a=[]
    for row in xrange(rows): a += [[None]*cols]
    return a
doPlot = True

durationInMin = 30 * 60

timeInMin = np.arange(0, durationInMin, 1)

numTeams = 6
numSignalCurvesPerTeam = 3

periodsInDays = np.array([[4.3, 7.1, 34.6], \
                     [2.3, 9.4, 26.8], \
                     [1.8, 13.7, 42.9], \
                     [3.7, 10.4, 46.3], \
                     [5.2, 14.6, 46.3], \
                     [5.2, 15.6, 31.7], \
                     [2.4, 17.3, 39.4]])

depthsInPercents = np.array([[1.2, 0.4, 1.9], \
                        [1.6, 0.7, 1.1], \
                        [1.3, 0.6, 2.1], \
                        [1.4, 0.5, 1.6], \
                        [1.5, 0.5, 1.1], \
                        [1.4, 0.6, 1.2]])
depths = depthsInPercents / 100.

impactParameters = np.array([[0.27, 0.41, 0.63], \
                        [0.39, 0.71, 0.23], \
                        [0.53, 0.31, 0.62], \
                        [0.67, 0.52, 0.21], \
                        [0.18, 0.36, 0.43], \
                        [0.68, 0.26, 0.54]])
T0sInMinutes = np.array([[  900.,  1000.,  1050.],
                      [  750.,   850.,   500.],
                      [  600.,   650.,   700.],
                      [ 1250.,  1100.,   750.],
                      [  500.,  1250.,  1000.],
                      [  550.,  1100.,   850.]]) / 2.

rhos = np.ones((numTeams , numSignalCurvesPerTeam))

signalCurves = twodlist(numTeams, numSignalCurvesPerTeam)

for t in range(numTeams):
    for s in range(numSignalCurvesPerTeam):
        cn = numSignalCurvesPerTeam * t + s

        signalCurves[t][s] =  \
            transitTools.transit( timeInMin,    #array of time values.
                                  rho = rhos[t, s],          # density of main star (solar density)
                                  ld1=0.2,            # limb darkening coefficients (2 coeff's for quadratic limb darkening)
                                  ld2=0.4,
                                  period = periodsInDays[t,s], # in days
                                  impact = impactParameters[t,s],  # between 0 and 1 yields a transit
                                  rprs = np.sqrt(depths[t,s]),     # radii ratio
                                  T0InMin = T0sInMinutes[t,s])



if doPlot:
    plt.figure("signal light curves", figsize = (20, 20))
    plt.clf()

    for t in range(numTeams):
        for s in range(numSignalCurvesPerTeam):
            plt.subplot(numTeams , numSignalCurvesPerTeam,t *  numSignalCurvesPerTeam +  s + 1)
            plt.plot(signalCurves[t][s], \
                         label = 'team %i, curve %i, P = %4.2f, d = %4.2f, i = %4.2f' %(t,s,periodsInDays[t,s], depthsInPercents[t,s], impactParameters[t,s], ))
            plt.ylim([.96, 1.02])
            plt.legend(loc = 'lower right')
    # plt.show()
    plt.savefig('../plot/signalCurves.pdf')  

for t in range(numTeams):
    for s in range(numSignalCurvesPerTeam):
        np.savetxt('../data/signalcurve_%i_%i.txt' %(t, s), signalCurves[t][s], header = 'data in one-minute intervals')
        np.savetxt('../data/signalperiod_%i_%i.txt' %(t, s), [periodsInDays[t,s]], header = 'period in days')

        

#########################################################################################
numNoiseCurvesPerTeam = 3

powerLawIndices = [0,1,2]

noiseCurves = twodlist(numTeams, numNoiseCurvesPerTeam)
# import pdb
# pdb.set_trace()

for t in range(numTeams):
    for s in range(numNoiseCurvesPerTeam):
        noiseCurves[t][s] = transitTools.noise(totalObservingTimeInMinutes = durationInMin, \
                                             samplingRatePerMinute = 1., \
                                             amplitude = 0.01, \
                                             powerLawIndex = powerLawIndices[s])



if doPlot:
    plt.figure("noise light curves", figsize = (20, 20))
    plt.clf()

    for t in range(numTeams):
        for s in range(numNoiseCurvesPerTeam):
            plt.subplot(numTeams , numNoiseCurvesPerTeam,t *  numNoiseCurvesPerTeam +  s + 1)
            plt.plot(noiseCurves[t][s])

            # plt.ylim([.96, 1.04])

        
plt.savefig('../plot/noiseCurves.pdf')  

for t in range(numTeams):
    for s in range(numNoiseCurvesPerTeam):
        np.savetxt('../data/noisecurve_%i_%i.txt' %(t, s), noiseCurves[t][s], header = 'data in one-minute intervals')
        
###########################################################################################

numGroups = 12
numSignalPlusNoiseCurves = 20

pRange = [2,30]

depthHighRange = np.array([1,1.5])
#depthMidRange = np.array([0.4, 0.7])
depthMidRange = np.array([0.7, 0.9])
depthLowRange = np.array([0.1, 0.3])

signalRanges = np.array([depthHighRange, depthMidRange, depthLowRange])
numSignalTypes = len(signalRanges)

noiseAmplitudes =  np.array([0.01, 0.02, 0.1])

impactRange = np.array([0, 0.8])
T0Range = np.array([500., 1300.]) / 2  #factor of 2 not understood.
# params = twodlist(numGroups, numSignalTypes)
signals = twodlist(numGroups, numSignalTypes)
noises = twodlist(numGroups, numSignalTypes)
signalPlusNoiseCurves = twodlist(numGroups, numSignalTypes)
signal_to_noise = twodlist(numGroups, numSignalTypes)
signal_to_noise_chi = twodlist(numGroups, numSignalTypes)

plt.figure('signal plus noise curves', figsize = (20, 40))
# plt.figure("noise light curves")
plt.clf()


for g in range(numGroups):
    for s in range(numSignalTypes):

        params = {'depth' : np.random.uniform(low = signalRanges[s,0], high = signalRanges[s,1]), 
                  'impact' : np.random.uniform(low = impactRange[0], high = impactRange[1]), 
                  'period' : np.random.uniform(low = pRange[0], high = pRange[1]) , 
                  'T0sInMinutes' : np.random.uniform(low = T0Range[0] , high = T0Range[1]), 
                  'noiseAmplitude' : noiseAmplitudes[s]}

        signals[g][s] = transitTools.transit( timeInMin,    #array of time values.
                                              rho = 1.,          # density of main star (solar density)
                                              ld1=0.2,            # limb darkening coefficients (2 coeff's for quadratic limb darkening)
                                              ld2=0.4,
                                              period = params['period'], # in days
                                              impact = params['impact'],  # between 0 and 1 yields a transit
                                              rprs = np.sqrt(params['depth'] / 100),    # radii ratio
                                              T0InMin = params['T0sInMinutes'])

        noises[g][s] = noiseAmplitudes[s] * np.random.randn(len(timeInMin))

        signalPlusNoiseCurves[g][s] = signals[g][s] + noises[g][s]

	#signal_to_noise[g][s] = transitTools.signal2noise_rootN(noises[g][s], params['depth'], transit_duration) # what to do for duration?) 

	#tophat model
	tophat = transitTools.tophat_model(T0InMin, params['depth'], transit_duration)
	trapezoid = transitTools.trapezoid_model(T0InMin, params['depth'], transit_duration, inner_duration) #how to calc durations?
	signal_to_noise_chi[g][s] = transitTools.signal2noise_chi(noises[g][s], signalPlusNoiseCurves[g][s], tophat)

        #add S/N to params


        plt.subplot(numGroups , numSignalTypes,g *  numSignalTypes +  s + 1)
        plt.plot(signalPlusNoiseCurves[g][s], \
                     label = 'team %i, curve %i, P = %4.2f, d = %4.2f, i = %4.2f' %(g,s,params['period'], params['depth'], params['impact'] ), color = 'b')

        plt.plot(signals[g][s], color = 'r')

        plt.legend(loc = 'lower right')


        pickle.dump(params, open('../data/params_%02i_%i.pkl' %(g, s), 'w'))
        
        np.savetxt('../data/signalPlusNoiseCurve_%02i_%i.txt' %(g, s), signalPlusNoiseCurves[g][s],  header = 'data in one-minute intervals')
        np.savetxt('../data/signalPlusNoisePeriod_%02i_%i.txt' %(g, s), [params['period']], header = 'signal period in days')

        


plt.savefig('../plot/signalPlusNoiseCurves.pdf')  




    
