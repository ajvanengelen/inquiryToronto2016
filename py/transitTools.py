


from numpy import *


import matplotlib.pyplot as plt

import pdb

def noise(totalObservingTimeInMinutes = 1800,  samplingRatePerMinute = 1., amplitude = 0.01, kneeFreq = 100., powerLawIndex = 0, doPlots = False):

    #time duration of each sample.
    lengthOfObservationInMinutes = 1. / samplingRatePerMinute

    N = int(np.round(totalObservingTimeInMinutes / lengthOfObservationInMinutes))

    #random, uncorrelated (white) noiwe
    uncorrelated = random.randn(N)  #white/uncorrelated Gaussian noise 

    #fourier transform it.
    fftUncorrelated = fft.fft(uncorrelated)

    #get frequencies corresponding to the Fourier transform.

    fftFreqs = fft.fftfreq(N, d = lengthOfObservationInMinutes)

    
    powerSpec = fftFreqs.copy()
    powerSpec[:] = 0
    #power spectrum white white and 1/f noise components.
    powerSpec[1:] = amplitude**2 *(1 +  abs( (fftFreqs[1:] / kneeFreq)**(-1 * powerLawIndex)))
    powerSpec[0] = amplitude**2

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

    # return (lengthOfObservationInMinutes * arange(float(N)), correlated)
    return correlated


import ktransit
import matplotlib.pyplot as plt
import numpy as np
from scipy import constants,optimize

def phased(lst,per=75.577294,T0=0.):
    return (lst.copy()-T0)/per%1

def center(lst):
    for i in range(len(lst)):
        if lst[i] > 0.5:
            lst[i] -= 1.
    return lst

# Define tophat model
def tophat_model(time, t0=900,depth=0.01,dur=2*60.): # top hat model function
    modelSignal = ones(1800)
    modelSignal[(time >t0-(dur/2.)) & (time < t0+(dur/2.))]=1-depth	    
    return modelSignal

# Define trapezoid function
def trapezoid_model(time, t0=900,depth=0.01,dur=2*60.,inner_dur=2*60/2.):
    modelSignal = np.ones(1800)
    modelSignal[(time > t0-(inner_dur/2.)) & (time < t0+(inner_dur/2.))]=1-depth
    angle_sides = np.linspace(start=1, stop=1-depth, num=(dur-inner_dur)/2., endpoint=False)
    modelSignal[(time > t0-(dur/2.)) & (time <= t0-(inner_dur/2.))]=angle_sides
    modelSignal[(time < t0+(dur/2.)) & (time >= t0+(inner_dur/2.))]=angle_sides[::-1]
    return modelSignal

def chi_squared(noise,SignalPlusNoiseCurve,modelSignal):
    chi_squared_value = np.sum(((SignalPlusNoiseCurve-modelSignal)/noise)**2)
    return chi_squared_value

def signal2noise_chi(noise,SignalPlusNoiseCurve,modelSignal):
    # chi square from tophat
    chi_squared_value = chi_squared(noise,SignalPlusNoiseCurve,modelSignal)  
    
    # chi square for the null hypothesis
    chi_squared_value_null = chi_squared(noise,SignalPlusNoiseCurve,ones(1800))  

    # calc S/N
    signal_to_noise = np.sqrt(np.abs(chi_squared_value_null - chi_squared_value))
    return signal_to_noise

def signal2noise_rootN(noise,depth,duration):
    signal_to_noise = (depth/noise)*np.sqrt(duration)
    return signal_to_noise

def transit( timeInMin,    #array of time values.
             rho=0.459,          # density of main star (solar density)
             ld1=0.2,            # limb darkening coefficients (2 coeff's for quadratic limb darkening)
             ld2=0.4,
             period = 1000, # in days
             impact=0,  # between 0 and 1 yields a transit
             rprs=0.4,    # radii ratio
             T0InMin = 1800. / 2., 
             ): #mid-transit time



    #
    #offset = 0
    #rvamp = 13.

    # Parameters
    G = 6.67259e-8     # (cgs)

    ld3=0.
    ld4=0.
    dil=0.
    zpt=1.             # zero point of transit data
    # veloffset=offset   # in same units as rv data

    ecosw=0.
    esinw=0.
    occ=0.
    # rvamp=rvamp        # in same units as rv data

    Rsun = 6.9598e10                   # radius of the sun (cgs)
    Msun = 1.989e33                    # mass of the sun (cgs)
    rhosun = Msun/(4./3*np.pi*Rsun**3) # density of the sun

    T0 = T0InMin/ (24. * 30.) 
    # Model
    M = ktransit.LCModel()

    M.add_star(
        rho=rho,
        ld1=ld1,
        ld2=ld2,
        ld3=ld3,
        ld4=ld4,
        dil=dil,
        zpt=zpt)
    #veloffset=offset)

    M.add_planet(
        T0=T0,
        period=period,
        impact=impact,
        rprs=rprs,
        ecosw=ecosw,
        esinw=esinw,
        occ=occ)
    #rvamp=rvamp)

    M.add_data(time = timeInMin.astype(float)/(60.*24.))
    #time=np.arange(0,1e3,1/24./60))

    tmod = M.transitmodel
    # (tmodtime, tmod)  = np.array(zip(*sorted(zip(center(phased(M.time)),tmod))))


    return tmod#(tmodtime, tmod)
