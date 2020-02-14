import numpy as np
import windows

N = 512
Frame = 35
analyzeSignal = np.zeros(N,dtype=complex)
InPhase = np.zeros(N)
Quadrature = np.zeros(N)
windowFunction = windows.window(N,name="Hanning")
for i in range(1,Frame):
    InPhase = np.loadtxt('../realPart/inPhase'+str(i)+'.txt')
    Quadrature = np.loadtxt('../imaginaryPart/Quadrature'+str(i)+'.txt')
    analyzeSignal = InPhase + Quadrature*1j
    analyzeSignal = analyzeSignal*windowFunction
    m = np.mean(analyzeSignal)
    a_f = np.fft.fftshift(np.fft.fft(analyzeSignal - m))
    np.savetxt('../ampWithWindow/Doppler'+str(i)+'.txt',np.abs(a_f))