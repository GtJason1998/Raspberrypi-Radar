import numpy as np
import matplotlib.pyplot as plt

N = 55
for i in range(1,N+1):
    Rres = np.loadtxt('./Rres/Rres'+str(i)+'.txt')
    Vres = np.loadtxt('./Vres/Vres'+str(i)+'.txt')
    sigDfft = np.loadtxt('./sigDfft/sigDfft'+str(i)+'.txt')
    plt.contourf(Rres,Vres,sigDfft)
    plt.title('Range-Doppler Map')
    plt.xlabel('Range/m')
    plt.ylabel('Velocity/mps')
    plt.savefig('./pic/jpeg/RDM'+str(i)+'.jpeg')
    plt.savefig('./pic/svg/RDM'+str(i)+'.svg')
    plt.close()