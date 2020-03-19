import numpy as np

N = 512
Frame = 35
analyzeSignal = np.zeros(N,dtype=complex)
InPhase = np.zeros(N)
Quadrature = np.zeros(N)
windowFunction = np.hamming(N)

for i in range(1,Frame):
    InPhase = np.loadtxt('../realPart/inPhase'+str(i)+'.txt')
    Quadrature = np.loadtxt('../imaginaryPart/Quadrature'+str(i)+'.txt')
    analyzeSignal = InPhase + Quadrature*1j
    m = np.mean(analyzeSignal)
    analyzeSignal = (analyzeSignal - m)*windowFunction
    a_f = np.fft.fftshift(np.fft.fft(analyzeSignal))
    np.savetxt('../ampWithWindow/Doppler'+str(i)+'.txt',np.abs(a_f))