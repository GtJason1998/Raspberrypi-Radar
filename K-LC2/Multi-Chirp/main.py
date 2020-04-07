#!/usr/bin/python3
# -*- coding:utf-8 -*-

import ADS1256
import DAC8532
import RPi.GPIO as GPIO
from pylab import fft,fftshift,datetime
import matplotlib.pyplot as plt
import numpy as np

N = 32
M = 10
Ivec = np.zeros(N)
Qvec = np.zeros(N)
ramp = np.linspace(0,1,N)

# Radar parameters
c0 = 299792458.              # Speed of light in vacuum
fs = 30000                   # Sample frequency
f1 = 24.05e9                 # Start frequency
f2 = 24.16e9                 # Stop frequency, using external amplifier
f0 = (f1 + f2)/2             # Center frequency
lambda0 = c0/f0              # Center wavelength
B = f2 - f1                  # Absolute bandwidth
Npad = 4                     # Factor to expand data vector with zeros
Mpad = 10                    # Factor to expand data vector with zeros in Doppler-FFT
V0 = 2                       # Maximum modulation voltage
windowFunction1 = np.hanning(N) # windowFunction
windowFunction2 = np.hanning(M) # windowFunction
Rres = c0/(2*B*Npad)

measureData = np.zeros((M,N),dtype=complex)
sigRWin = np.zeros((M,N),dtype=complex)
sigDWin = np.zeros((M,N*Npad),dtype=complex)
sigRfft = np.zeros((M,N*Npad),dtype=complex)
sigDfft = np.zeros((M*Mpad,N*Npad),dtype=complex)
analyzeSignal1 = np.zeros(N,dtype=complex)

# Define a convenient signal processing function
def EstimateAmplitudeOffset(f, f0):
    """Estimate amplitude and offset of f0 to minimize distance to f."""
    A = sum(f0*f0)
    B = sum(f0)
    C = sum(f0*f)
    D = sum(f)
    N = len(f)
    amplitude = 1/(A*N - B**2)*(N*C - B*D)
    offset = 1/(A*N - B**2)*(-B*C + A*D)
    return(amplitude, offset)

try:
    ADC = ADS1256.ADS1256()
    DAC = DAC8532.DAC8532()
    ADC.ADS1256_init()
    DAC.DAC8532_Out_Voltage(DAC8532.channel_A, 0)
    DAC.DAC8532_Out_Voltage(DAC8532.channel_B, 0)
    
    j = 0
    while(1):
        j = j + 1
        T1 = 0                       # Up_Chirp Time
        # Up-chirp pulse
        for i in range(0,M):
            for n in np.arange(0,N):
                StartTime1 = datetime.datetime.now()         # Get start time 
                DAC.DAC8532_Out_Voltage(DAC8532.channel_A,V0*n/N)
                EndTime1 = datetime.datetime.now()           # Get end time
                T1 = T1 + (EndTime1 - StartTime1).total_seconds() # Compute total time
                ADC_Value = ADC.ADS1256_GetIQ()
                Ivec[n] = ADC_Value[0]*5.0/0x7fffff
                Qvec[n] = ADC_Value[1]*5.0/0x7fffff
            #print('i = ',i)
            analyzeSignal1 = Ivec + 1j*Qvec
            amplitude, offset = EstimateAmplitudeOffset(analyzeSignal1, ramp)
            analyzeSignal1 = analyzeSignal1 - (amplitude*ramp + offset)
            measureData[i,:] = analyzeSignal1
        T1 = T1/M
        print('T1 = ',T1)
        Vres = c0/(2*f0*T1*M*Mpad)

        for ii in range(M):
            sigRWin[ii,:] = measureData[ii,:]*windowFunction1
        
        for ii in range(M):
            sigRfft[ii,:] = fft(sigRWin[ii,:],N*Npad)
        
        for ii in range(N*Npad):
            sigDWin[:,ii] = sigRfft[:,ii]*windowFunction2
        
        for ii in range(N*Npad):
            sigDfft[:,ii] = fftshift(fft(sigDWin[:,ii],M*Mpad))
        '''
        plt.figure(1)
        plt.contourf(Rres*np.arange(1,N*Npad+1),Vres*(np.arange(1,M*Mpad+1) - M*Mpad/2),np.abs(sigDfft))
        plt.xlabel('Range/m')
        plt.ylabel('Velocity/mps')
        plt.pause(0.000001)
        plt.cla()
        '''
        '''
        np.savetxt('./Rres/Rres'+str(j)+'.txt',Rres*np.arange(1,N*Npad+1))
        np.savetxt('./Vres/Vres'+str(j)+'.txt',Vres*(np.arange(1,M*Mpad+1) - M*Mpad/2))
        np.savetxt('./sigDfft/sigDfft'+str(j)+'.txt',np.abs(sigDfft))
        '''
        
        plt.figure(1)
        plt.contourf(Rres*np.arange(1,N*Npad+1),np.arange(1,M+1),np.abs(sigRfft))
        plt.xlabel('Range/m')
        plt.ylabel('Frame')
        plt.pause(0.000001)
        plt.cla()
        
except Exception as e:
    print(e)
    GPIO.cleanup()
    print ("\r\nProgram end     ")
    exit()
