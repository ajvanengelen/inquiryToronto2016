


from numpy import *


import matplotlib.pyplot as plt

import pdb

def timeSeries(N = 300, samplingRatePerDay = 1., amplitude = 1., kneeFreq = 1, powerLawIndex = 1, doPlots = False):

    #time duration of each sample.
    lengthOfObservation = 1. / samplingRatePerDay

    #duration of the whold dataset
    duration = N * lengthOfObservation

    #random, uncorrelated (white) noiwe
    uncorrelated = random.randn(N)  #white/uncorrelated Gaussian noise 

    #fourier transform it.
    fftUncorrelated = fft.fft(uncorrelated)

    #get frequencies corresponding to the Fourier transform.
    fftFreqs = fft.fftfreq(N, d = lengthOfObservation)

    #power spectrum white white and 1/f noise components.
    powerSpec = amplitude *(1 +  abs( (fftFreqs / kneeFreq)**(-1 * powerLawIndex)))
    powerSpec[0] = 0.

    #multiply by square root of the power spectrum.
    fftCorrelated = fftUncorrelated * sqrt(abs(powerSpec))


    #make plots if requested..
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

    
        #fourier transform back to real space.
    correlated = real(fft.ifft(nan_to_num(fftCorrelated)))

    # return (correlated, uncorrelated, [fftFreqs, powerSpec])

    return (lengthOfObservation * arange(float(N)), correlated)
