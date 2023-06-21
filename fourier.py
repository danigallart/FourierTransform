import numpy as np
import matplotlib.pyplot as plt

#Boundary conditions 
rangt = 500                 #number of points for time domain
rangw = 500                 #number of points for frequency domain
wi = -2                       #initial frequency of the frequency domain
wf = 2                        #final frequency of the freqeucny domain
ti = -100                       #initial time of the time domain
tf = 100                        #final time of the time domain

            
t = np.linspace(ti,tf,rangt)  #time
#f = np.cos(2.0*np.pi*(3.0*t))*np.exp(-np.pi*t**2) #function to be transformed
f = np.sin(t)
real_function = np.zeros(rangt) #in principle we should recover f if it was real
im_function = np.zeros(rangt) #imaginary part of the antitransform

real_fourier = np.zeros(rangw) #real part of the Fourier transform
im_fourier = np.zeros(rangw)   #imaginary part of the Fourier transform
w = np.linspace(wi,wf,rangw)   #frequency

def ftransform(f,t,rangt,rangw,w,j):
    reftransform = 0.0
    #dt = np.diff(t)
    dt = (tf-ti)/rangt

    for i in range(1,rangt):
        reftransform += f[i]*np.cos(2.*np.pi*w[j]*t[i])*dt
    reftransform += dt*(f[0]*np.cos(2.*np.pi*w[j]*t[i])/2. + f[rangt-1]*np.cos(2.*np.pi*w[j]*t[i])/2.0)

    imftransform = 0.0

    for i in range(1,rangt):
        imftransform += f[i]*np.sin(-2.*np.pi*w[j]*t[i])*dt
    imftransform += dt*(f[0]*np.sin(-2.*np.pi*w[j]*t[i])/2. + f[rangt-1]*np.sin(-2.*np.pi*w[j]*t[i])/2.0)

    return reftransform, imftransform

def anti_ftransform(real_fourier,im_fourier,rangw,rangt,w,t,j):
    anti_reftransform = 0.0
    dw = (wf-wi)/rangw

    for i in range(1,rangw):
        anti_reftransform += real_fourier[i]*np.cos(2.*np.pi*w[i]*t[j])*dw + -im_fourier[i]*np.sin(2.*np.pi*w[i]*t[j])*dw
    
    anti_reftransform += (real_fourier[0]*np.cos(2.*np.pi*w[0]*t[j])/2. + real_fourier[rangw-1]*np.cos(2.*np.pi*w[rangw-1]*t[j])/2.)*dw + (-im_fourier[0]*np.sin(2.*np.pi*w[0]*t[j])/2. - im_fourier[rangw-1]*np.sin(2.*np.pi*w[rangw-1]*t[j])/2.)*dw

    anti_imftransform = 0.0

    for i in range(1,rangw):
        anti_imftransform += (real_fourier[i]*np.sin(2.*np.pi*w[i]*t[j]))*dw + (im_fourier[i]*np.cos(2.*np.pi*w[i]*t[j]))*dw

    anti_imftransform += (real_fourier[0]*np.sin(2.*np.pi*w[0]*t[j])/2.0 + real_fourier[rangw-1]*np.sin(2.*np.pi*w[rangw-1]*t[j])/2.)*dw + (im_fourier[0]*np.cos(2.*np.pi*w[0]*t[j])/2. + im_fourier[rangw-1]*np.cos(2.*np.pi*w[rangw-1]*t[j])/2.0)*dw

    return anti_reftransform, anti_imftransform



#Initiate the Fourier transform

for i in range(rangw):
    real_fourier[i], im_fourier[i] = ftransform(f,t,rangt,rangw,w,i) 

for i in range(rangt):
    real_function[i], im_function[i] = anti_ftransform(real_fourier, im_fourier, rangw, rangt, w, t, i)
