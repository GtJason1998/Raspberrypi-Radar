'''
Create Window Function
including Rectangle,Triangle,Hanning,Hamming,Blackman
'''
import numpy as np

def window(N,name = "Rectangle"):
    '''
    window(N,name = "Rectangle/Triangle/Hanning/Hamming/Blackman")
    '''
    if name == "Rectangle":
        window = np.ones(N)
    if name == "Triangle":
        temp = np.array([2*i/(N - 1) for i in np.arange(int((N+1)/2))])
        window = np.append(temp,np.array([2-(2*i)/(N-1) for i in np.arange(int((N+1)/2),N)]))
    if name == "Hanning":
        window = np.array([0.5*(1 - np.cos((2*np.pi*i)/(N-1))) for i in np.arange(N)])
    if name == "Hamming":
        window = np.array([0.54 - 0.46*np.cos((2*np.pi*i)/(N-1)) for i in np.arange(N)])
    if name == "Blackman":
        window = np.array([0.42 - 0.5*np.cos((2*np.pi*i)/(N-1)) + 0.08*np.cos((4*np.pi*i)/(N-1)) for i in np.arange(N)])
    return window