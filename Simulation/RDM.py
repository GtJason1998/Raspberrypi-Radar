import numpy as np
import matplotlib.pyplot as plt
from pylab import fft,fftshift

tarR = np.array([10,20,30])
tarV = np.array([0,-5,10])
c = 3*10**8
f0 = 24.05*10**9
T = 0.0002          #chirp Sweep Time
B = 400*10**6
L = 24              #slow-time dimension,num of chirps
N = 128             #fast-time dimension,num of samples
Npad = 10           #padding in order to improve measure precision
Lpad = 60           #padding in order to improve measure precision

##generate receive signal
S1 = np.zeros((L,N),dtype = complex)
for l in range(L):
    for n in range(N):
        S1[l,n] = 2*np.exp(1j*2*np.pi*((2*B*(tarR[0]+tarV[0]*T*(l+1))/(c*T)+(2*f0*tarV[0]/c))*T/N*(n+1)+((2*f0)*(tarR[0]+tarV[0]*T*(l+1)))/c))
        
S2 = np.zeros((L,N),dtype = complex)
for l in range(L):
    for n in range(N):
        S2[l,n] = 2*np.exp(1j*2*np.pi*((2*B*(tarR[1]+tarV[1]*T*(l+1))/(c*T)+(2*f0*tarV[1]/c))*T/N*(n+1)+((2*f0)*(tarR[1]+tarV[1]*T*(l+1)))/c))

S3 = np.zeros((L,N),dtype = complex)
for l in range(L):
    for n in range(N):
        S3[l,n] = 2*np.exp(1j*2*np.pi*((2*B*(tarR[2]+tarV[2]*T*(l+1))/(c*T)+(2*f0*tarV[2]/c))*T/N*(n+1)+((2*f0)*(tarR[2]+tarV[2]*T*(l+1)))/c))

sigReceive = S1 + S2 + S3
sigRWin = np.zeros((L,N),dtype=complex)
sigDWin = np.zeros((L,N*Npad),dtype=complex)
sigRfft = np.zeros((L,N*Npad),dtype=complex)
sigDfft = np.zeros((L*Lpad,N*Npad),dtype=complex)

windowFunction1 = np.hanning(N) # windowFunction
windowFunction2 = np.hanning(L) # windowFunction

for ii in range(L):
    sigRWin[ii,:] = sigReceive[ii,:]*windowFunction1
        
for ii in range(L):
    sigRfft[ii,:] = fft(sigRWin[ii,:],N*Npad)
        
for ii in range(N*Npad):
    sigDWin[:,ii] = sigRfft[:,ii]*windowFunction2
        
for ii in range(N*Npad):
    sigDfft[:,ii] = fftshift(fft(sigDWin[:,ii],L*Lpad))

Rres = c/(2*B*Npad)
Vres = c/(2*24.25*10**9*T*L*Lpad)
plt.contourf(Rres*np.arange(1,N*Npad+1),Vres*(np.arange(1,L*Lpad+1) - L*Lpad/2),np.abs(sigDfft))
plt.xlabel('Range/m')
plt.ylabel('Velocity/mps')
plt.title('Range-Doppler Map')
plt.minorticks_on()
plt.show()