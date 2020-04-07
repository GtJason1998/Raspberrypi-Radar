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
analyzeSignal = zeros(N,dtype=complex)
C = 3*10**8
f_vco = 24.15*10**9
windowFunction = np.hamming(N)

try:
    ADC = ADS1256.ADS1256()
    DAC = DAC8532.DAC8532()
    ADC.ADS1256_init()

    DAC.DAC8532_Out_Voltage(DAC8532.channel_A,1)
    DAC.DAC8532_Out_Voltage(DAC8532.channel_B,1)

    while(1):
        #starttime = datetime.datetime.now()
        for i in range(0,N):
            ADC_Value = ADC.ADS1256_GetIQ()
            I = ADC_Value[0]*5.0/0x7fffff
            Q = ADC_Value[1]*5.0/0x7fffff
            analyzeSignal[i] = I + 1j*Q
            #print("0 ADC = %lf"%(I))
            #print("1 ADC = %lf"%(Q))
            #print("\33[9A")
        #endtime = datetime.datetime.now()
        #T = (endtime - starttime).total_seconds()
        analyzeSignal = analyzeSignal*2
        m = mean(analyzeSignal)
        analyzeSignal = (analyzeSignal - m)*windowFunction
        #print("mean = ",m)
        #print("fs = ",N/T," Hz")
        a_f = fftshift(fft(analyzeSignal,4*N))
        fvec = fftshift(fftfreq(4*N,d = 1/30000))
        f = fvec[abs(a_f) == max(abs(a_f))][0]
        v = (f*C)/(2*f_vco)
        print("Doppler frequency = {0:10.2f} Hz".format(f))
        print("Velocity          = {0:10.2f} m/s".format(v))
        plt.plot(fvec,10*np.log10(abs(a_f)))
        #plt.plot(real(analyzeSignal))
        plt.grid()
        plt.pause(0.001)
        plt.cla()
        #np.savetxt('./amplitude/Doppler.txt',abs(a_f))
        #np.savetxt('./realPart/inPhase.txt',real(analyzeSignal))
        #np.savetxt('./imaginaryPart/Quadrature.txt',imag(analyzeSignal))

except :
    #plt.ioff()
    #plt.show()
    #np.savetxt('fvec.txt',fvec)
    GPIO.cleanup()
    print ("\r\nProgram end     ")
    exit()
