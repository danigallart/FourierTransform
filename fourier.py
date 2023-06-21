import numpy as np
import matplotlib.pyplot as plt


def ftransform(f,t,rangt,rangw,w,j):
    reftransform = 0.0
    #dt = np.diff(t)
    dt = (t[-1]-t[0])/rangt

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
    dw = (w[-1]-w[0])/rangw

    for i in range(1,rangw):
        anti_reftransform += real_fourier[i]*np.cos(2.*np.pi*w[i]*t[j])*dw + -im_fourier[i]*np.sin(2.*np.pi*w[i]*t[j])*dw
    
    anti_reftransform += (real_fourier[0]*np.cos(2.*np.pi*w[0]*t[j])/2. + real_fourier[rangw-1]*np.cos(2.*np.pi*w[rangw-1]*t[j])/2.)*dw + (-im_fourier[0]*np.sin(2.*np.pi*w[0]*t[j])/2. - im_fourier[rangw-1]*np.sin(2.*np.pi*w[rangw-1]*t[j])/2.)*dw

    anti_imftransform = 0.0

    for i in range(1,rangw):
        anti_imftransform += (real_fourier[i]*np.sin(2.*np.pi*w[i]*t[j]))*dw + (im_fourier[i]*np.cos(2.*np.pi*w[i]*t[j]))*dw

    anti_imftransform += (real_fourier[0]*np.sin(2.*np.pi*w[0]*t[j])/2.0 + real_fourier[rangw-1]*np.sin(2.*np.pi*w[rangw-1]*t[j])/2.)*dw + (im_fourier[0]*np.cos(2.*np.pi*w[0]*t[j])/2. + im_fourier[rangw-1]*np.cos(2.*np.pi*w[rangw-1]*t[j])/2.0)*dw

    return anti_reftransform, anti_imftransform

def plot_fourier(f, t, w, real_fourier, im_fourier, real_function, im_function):
    fig, axes = plt.subplots(2,3,figsize=(20,20))
    axes[0,0].plot(t,f)
    axes[0,0].set_title('Function to transform')
    axes[0,0].set_xlabel('time (s)')

    axes[0,1].plot(2.*np.pi*w,real_fourier)
    axes[0,1].set_title('Re[Fourier]')
    axes[0,1].set_xlabel('w (rad/s)')

    axes[0,2].plot(2.*np.pi*w,im_fourier)
    axes[0,2].set_title('Im[Fourier]')
    axes[0,2].set_xlabel('w (rad/s)')

    axes[1,0].plot(t,real_function)
    axes[1,0].set_title('Re[AntiFourier]')
    axes[1,0].set_xlabel('time (s)')

    axes[1,1].plot(t,im_function)
    axes[1,1].set_title('Im[AntiFourier]')
    axes[1,1].set_xlabel('time (s)')

    axes[1,2].plot(t,f)
    axes[1,2].plot(t,real_function)
    axes[1,2].plot(t,im_function)
    axes[1,2].set_title('Comparison Anti vs Function')
    axes[1,2].set_xlabel('time (s)')


    plt.show()

