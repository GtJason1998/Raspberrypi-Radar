import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#目标距离和速度
R = 20
v = 10

f0 = 24.25*10**9    #发射信号频率
fs = 7500           #采样频率
A0 = 1              #发射信号幅值
C = 3*10**8         #光速
#添加高斯白噪声
def wgn(x, snr):
    P_signal = np.sum(abs(x)**2)/len(x)
    P_noise = P_signal/10**(snr/10.0)
    return np.random.randn(len(x)) * np.sqrt(P_noise)
#构建IQ分析信号
def analySignal(t):
    return np.array([(A0/2)*np.cos(2*np.pi*f0*2*v*t/C - 2*np.pi*f0*2*R/C) + 
    1j*(A0/2)*np.cos(2*np.pi*f0*2*v*t/C - 2*np.pi*f0*2*R/C - np.pi/2)])

t = np.arange(0,2,1/fs)
s_receive = analySignal(t)[0,:]
s_receive = s_receive + wgn(s_receive,0.1)
fvec = np.fft.fftshift(np.fft.fftfreq(s_receive.shape[0],1/fs))
S_k = np.fft.fftshift(np.fft.fft(s_receive))
sns.set(style='ticks')
plt.plot(fvec,abs(S_k))
plt.grid()
plt.xlabel('frequency/Hz')
plt.ylabel('Amplitude')
plt.text(1616.5,abs(S_k)[10733],(1616.5,abs(S_k)[10733]),ha='center',va='bottom',fontsize=10)
plt.show()