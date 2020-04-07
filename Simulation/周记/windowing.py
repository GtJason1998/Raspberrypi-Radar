import numpy as np
import matplotlib.pyplot as plt
from pylab import mpl

mpl.rcParams['font.sans-serif'] = ['SimHei'] 
mpl.rcParams['axes.unicode_minus'] = False  
N = 128
fig0,ax0 = plt.subplots(1,1)
rectWindow = np.ones(N)
ax0.plot(np.concatenate((np.zeros(10),rectWindow,np.zeros(10))))
rectWindow = np.concatenate((np.zeros(1000),rectWindow,np.zeros(1000)))
rectWindow_w = np.fft.fftshift(np.fft.fft(rectWindow))
fig1,ax1 = plt.subplots(1,1)
ax1.plot(abs(rectWindow_w))

hanningWindow = np.hanning(N)
ax0.plot(np.concatenate((np.zeros(10),hanningWindow,np.zeros(10))))
hanningWindow = np.concatenate((np.zeros(1000),hanningWindow,np.zeros(1000)))
hanningWindow_w = np.fft.fftshift(np.fft.fft(hanningWindow))
ax1.plot(abs(hanningWindow_w))
hammingWindow = np.hamming(N)
ax0.plot(np.concatenate((np.zeros(10),hammingWindow,np.zeros(10))))
hammingWindow = np.concatenate((np.zeros(1000),hammingWindow,np.zeros(1000)))
hammingWindow_w = np.fft.fftshift(np.fft.fft(hammingWindow))
ax1.plot(abs(hammingWindow_w))
blackmanWindow = np.blackman(N)
ax0.plot(np.concatenate((np.zeros(10),blackmanWindow,np.zeros(10))))
blackmanWindow = np.concatenate((np.zeros(1000),blackmanWindow,np.zeros(1000)))
blackmanWindow_w = np.fft.fftshift(np.fft.fft(blackmanWindow))
ax1.plot(abs(blackmanWindow_w))
ax1.legend(['矩形窗','汉宁窗','汉明窗','布莱克曼窗'])
ax1.minorticks_on()
ax1.grid()
ax1.set_title('F(w)')

ax0.legend(['矩形窗','汉宁窗','汉明窗','布莱克曼窗'])
ax0.minorticks_on()
ax0.grid()
ax0.set_ylim(0,1.35)
ax0.set_title('f(t)')
plt.show()