


from numpy import *


import matplotlib.pyplot as plt

import pdb

def timeSeries(N = 300, samplingRatePerDay = 1., amplitude = 1., kneeFreq = 1, powerLawIndex = 1, doPlots = False):


    lengthOfObservation = 1. / samplingRatePerDay

    duration = N * lengthOfObservation

    uncorrelated = random.randn(N)  #white/uncorrelated Gaussian noise 
   
    fftUncorrelated = fft.fft(uncorrelated)

    dOmega = 2 * pi / duration

    fftFreqs = fft.fftfreq(N, d = lengthOfObservation)

    powerSpec = amplitude *(1 +  abs( (fftFreqs / kneeFreq)**(-1 * powerLawIndex)))
    powerSpec[0] = 0.

    fftCorrelated = fftUncorrelated * sqrt(abs(powerSpec))


    if doPlots:
        plt.figure('power')
        plt.clf()
        plt.subplot(2,1,1)
        # plt.plot(fftFreqs, powerSpec)
        
        # plt.subplot(2,1,2)
        # plt.plot(fftFreqs, abs(fftCorrelated))
        plt.plot(powerSpec)
        
        plt.subplot(2,1,2)
        plt.plot(abs(fftCorrelated))
        
        plt.show()
        
        
        # pdb.set_trace()

    

    correlated = real(fft.ifft(nan_to_num(fftCorrelated)))

    # return (correlated, uncorrelated, [fftFreqs, powerSpec])

    return (lengthOfObservation * arange(float(N)), correlated)
