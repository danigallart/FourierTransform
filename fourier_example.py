import numpy as np
import matplotlib.pyplot as plt
import fourier as fou

#Boundary conditions 
rangt = 1000                 #number of points for time domain
rangw = 1000                 #number of points for frequency domain
wi = -6                       #initial frequency of the frequency domain
wf = 6                        #final frequency of the freqeucny domain
ti = -5                       #initial time of the time domain
tf = 5                        #final time of the time domain

            
t = np.linspace(ti,tf,rangt)  #time
f = np.cos(2.0*np.pi*(3.0*t))*np.exp(-np.pi*t**2) #function to be transformed
#f = np.sin(t)
real_function = np.zeros(rangt) #in principle we should recover f if it was real
im_function = np.zeros(rangt) #imaginary part of the antitransform

real_fourier = np.zeros(rangw) #real part of the Fourier transform
im_fourier = np.zeros(rangw)   #imaginary part of the Fourier transform
w = np.linspace(wi,wf,rangw)   #frequency


#Initiate the Fourier transform

for i in range(rangw):
    real_fourier[i], im_fourier[i] = fou.ftransform(f,t,rangt,rangw,w,i) 

for i in range(rangt):
    real_function[i], im_function[i] = fou.anti_ftransform(real_fourier, im_fourier, rangw, rangt, w, t, i)

fou.plot_fourier(f, t, w, real_fourier, im_fourier, real_function, im_function)