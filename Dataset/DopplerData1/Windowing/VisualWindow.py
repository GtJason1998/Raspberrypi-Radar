import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()
fvec = np.loadtxt('../fvec.txt')
T = 1/7500
t = np.arange(0,fvec.shape[0]) * T
P_fa = 10**(-6)     # false alert rate
n = 16
windowSize = n + 1
shieldCell = 4
referenceCell = n - shieldCell
factor_T = P_fa**(-1/referenceCell) - 1   
print(factor_T)
CFAR_xLabel = fvec[int(n/2)-1:fvec.shape[0]-int(n/2)].copy()
CFAR_Level = np.zeros(CFAR_xLabel.shape[0],dtype=float)

N = 35
for i in range(1,N):
    Doppler = np.loadtxt('../ampWithWindow/Doppler'+str(i)+'.txt')
    for j in range(CFAR_xLabel.shape[0]):
        CFAR_Level[j] = np.sum(Doppler[j:j+windowSize].copy())
        for k in range(1,int(shieldCell/2)+1):
            CFAR_Level[j] = CFAR_Level[j] - Doppler[j+int(n/2)-k] - Doppler[j+int(n/2)+k]
        CFAR_Level[j] = CFAR_Level[j] - Doppler[j+int(n/2)]
        CFAR_Level[j] = CFAR_Level[j]/referenceCell
    #inPhase = np.loadtxt('./realPart/inPhase'+str(i)+'.txt')
    #Quadrature = np.loadtxt('./imaginaryPart/Quadrature'+str(i)+'.txt')
    CFAR_Level = CFAR_Level*factor_T
    plt.plot(fvec/1000,Doppler)
    plt.plot(CFAR_xLabel/1000,CFAR_Level)
    plt.title('Amplitude'+str(i))
    plt.legend(['Amplitude','CA-CFAR'])
    plt.xlabel('frequency/KHz')
    plt.ylabel('Amplitude')
    plt.savefig('../pic/ampWithWindow_CFAR/Doppler'+str(i)+'.png')
    plt.close()

    plt.plot(fvec/1000,Doppler)
    plt.title('Amplitude'+str(i))
    plt.xlabel('frequency/KHz')
    plt.ylabel('Amplitude')
    plt.savefig('../pic/ampWithWindow/Doppler'+str(i)+'.png')
    plt.close()