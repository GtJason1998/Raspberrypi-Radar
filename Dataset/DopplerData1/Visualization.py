import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style('whitegrid')
fvec = np.loadtxt('./fvec.txt')
T = 1/7500
t = np.arange(0,fvec.shape[0]) * T

N = 35
for i in range(1,N):
    Doppler = np.loadtxt('./amplitude/Doppler'+str(i)+'.txt')
    inPhase = np.loadtxt('./realPart/inPhase'+str(i)+'.txt')
    Quadrature = np.loadtxt('./imaginaryPart/Quadrature'+str(i)+'.txt')

    plt.plot(fvec/1000,Doppler)
    plt.title('Amplitude'+str(i))
    plt.xlabel('frequency/KHz')
    plt.ylabel('Amplitude')
    plt.savefig('./pic/amplitude/Doppler'+str(i)+'.png')
    plt.close()

    plt.plot(t*1000,inPhase)
    plt.title('inPhase'+str(i))
    plt.xlabel('t/ms')
    plt.ylabel('Amplitude/v')
    plt.savefig('./pic/inPhase/inPhase'+str(i)+'.png')
    plt.close()

    plt.plot(t*1000,Quadrature)
    plt.title('Quadrature'+str(i))
    plt.xlabel('t/ms')
    plt.ylabel('Amplitude/v')
    plt.savefig('./pic/Quadrature/Quadrature'+str(i)+'.png')
    plt.close()