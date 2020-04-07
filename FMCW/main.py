#!/usr/bin/python3
# -*- coding:utf-8 -*-

import time 
import ADS1256
import DAC8532
import RPi.GPIO as GPIO
from pylab import *
import matplotlib.pyplot as plt
import numpy as np

N = 128
Ivec = np.zeros(N)
Qvec = np.zeros(N)
ramp = np.linspace(0,1,N)

# Radar parameters
c0 = 299792458.              # Speed of light in vacuum
fs = 7500                    # Sample frequency
F1 = 24.05e9                 # Start frequency
F2 = 24.45e9                 # Stop frequency, using external amplifier
f0 = (F1 + F2)/2             # Center frequency
lambda0 = c0/f0              # Center wavelength
B = F2 - F1                  # Absolute bandwidth
Nave = 1                     # Number of averaging pulses
Npad = 4                     # Factor to expand data vector with zeros
V0 = 2                       # Maximum modulation voltage
windowFunction = np.hamming(N) #windowFunction
f1 = fftshift(fftfreq(Npad*N, d=1/fs)) #Up-Chirp XLabel in frequency domain
f2 = fftshift(fftfreq(Npad*N, d=1/fs)) #Down-Chirp XLabel in frequency domain

#CFAR
P_fa = 10**(-6)
n_CFAR = 16
windowSize = n_CFAR + 1
shieldCell = 2
referenceCell = n_CFAR - shieldCell
factor_T = P_fa**(-1/referenceCell) - 1
CFAR_xLabel = f1[int(n_CFAR/2)-1:f1.shape[0]-int(n_CFAR/2)].copy()
CFAR_Level = np.zeros(CFAR_xLabel.shape[0],dtype=float)
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

    while(1):
        T1 = 0                       # Up_Chirp Time
        T2 = 0                       # Down_Chirp Time
        analyzeSignal1 = np.zeros(N,dtype=complex)
        analyzeSignal2 = np.zeros(N,dtype=complex)
        for m in np.arange(0,Nave):
             # Up-chirp pulse
            for n in np.arange(0,N):
                StartTime1 = datetime.datetime.now()         # Get start time 
                DAC.DAC8532_Out_Voltage(DAC8532.channel_A,V0*n/N)
                EndTime1 = datetime.datetime.now()           # Get end time
                T1 = T1 + (EndTime1 - StartTime1).total_seconds() # Compute total time
                ADC_Value = ADC.ADS1256_GetIQ()
                Ivec[n] = ADC_Value[0]*5.0/0x7fffff
                Qvec[n] = ADC_Value[1]*5.0/0x7fffff
            analyzeSignal1 = analyzeSignal1 + (Ivec + 1j*Qvec)/Nave
            

            # Down-chirp pulse
            for n in np.arange(0,N):
                StartTime2 = datetime.datetime.now()         # Get start time
                DAC.DAC8532_Out_Voltage(DAC8532.channel_A,V0*(N - n - 1)/N)
                EndTime2 = datetime.datetime.now()           # Get end time
                T2 = T2 + (EndTime2 - StartTime2).total_seconds() # Compute total time
                ADC_Value = ADC.ADS1256_GetIQ()
                Ivec[n] = ADC_Value[0]*5.0/0x7fffff
                Qvec[n] = ADC_Value[1]*5.0/0x7fffff
            analyzeSignal2 = analyzeSignal2 + (Ivec + 1j*Qvec)/Nave
            #T_average = (T1 + T2)/2
            #print("T_average = ",T_average)

        # Estimate the amplitude and offset of the ramp in up-chirp and down-chirp,
        # then subtract the result
        a10 = analyzeSignal1  # Save raw analytical signal for later plotting
        a20 = analyzeSignal2  # Save raw analytical signal for later plotting
        amplitude, offset = EstimateAmplitudeOffset(analyzeSignal1, ramp)
        analyzeSignal1 = analyzeSignal1 - (amplitude*ramp + offset)
        amplitude, offset = EstimateAmplitudeOffset(analyzeSignal2, ramp)
        analyzeSignal2 = analyzeSignal2 - (amplitude*ramp + offset)

        analyzeSignal1 = 5*analyzeSignal1*windowFunction
        analyzeSignal2 = 5*analyzeSignal2*windowFunction

        # Convert to frequency domain and analyze
        # The 'Npad' factor enables zero-padding the data for interpolation and nicer plots
        a1_f = fftshift(fft(analyzeSignal1, n=Npad*N))
        a2_f = fftshift(fft(analyzeSignal2, n=Npad*N))
        A1_f = abs(a1_f)
        A2_f = abs(a2_f)
        
        #CFAR
        for j in np.arange(CFAR_xLabel.shape[0]):
            #print("j+windowSize = ",j+windowSize)
            CFAR_Level[j] = np.sum(A1_f[j:j+windowSize].copy())
            for k in np.arange(1,int(shieldCell/2)+1):
                CFAR_Level[j] = CFAR_Level[j] - A1_f[j+int(n_CFAR/2)-k] - A1_f[j+int(n_CFAR/2)+k]
            CFAR_Level[j] = CFAR_Level[j] - A1_f[j+int(n_CFAR/2)]
            CFAR_Level[j] = CFAR_Level[j]/referenceCell
        CFAR_Level = CFAR_Level*factor_T
        # Compute the delta frequencies corresponding to the maximum peaks 
        # in up-chirp and down-chirp
        temp = f1[int(n_CFAR/2)-1:f1.shape[0]-int(n_CFAR/2)]
        #targetNum = temp[A1_f[int(n_CFAR/2)-1:f1.shape[0]-int(n_CFAR/2)] > CFAR_Level].shape[0]
        #print("target_num = ",targetNum)
        
        df1 = f1[abs(a1_f) == max(abs(a1_f))][0]
        df2 = f2[abs(a2_f) == max(abs(a2_f))][0]
        #print("f_b = ",abs(df1 - df2)/2)
        #print("f_d = ",(df1 + df2)/2)
        #print("f1 = ",df1)
        #print("f2 = ",df2)
        
        plt.figure(1)
        plt.plot(f1,(abs(a1_f)))
        plt.plot(temp,(CFAR_Level))
        #plt.plot(real(analyzeSignal1))
        plt.pause(0.001)
        plt.cla()
        plt.figure(2)
        plt.plot(f2,abs(a2_f))
        plt.pause(0.001)
        plt.cla()
        
        R = ((c0)/(4*B))*abs(T2*df2 + T1*df1)
        v = (c0/(4*f0))*(df1 - df2)
        

        print("R = ",R)
        print("v = ",v)


except Exception as e:
    print(e)
    GPIO.cleanup()
    print ("\r\nProgram end     ")
    exit()
