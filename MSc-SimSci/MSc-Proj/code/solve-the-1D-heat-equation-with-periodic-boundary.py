import numpy as np


T_final = 5

# Physical parameters:

L=1
D=1

# *************************************************************************
# Numerical parameters:

N=64
dx=L/float(N)


def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step
        
x=[i for i in frange(0,L,dx)]
dt=0.01
N_timesteps=int(T_final/dt)


u0 = np.array([1 if x[i] >= .25 and x[i] <= .75 else 0 for i in range(0,N) ])

M=np.zeros((N_timesteps+1,N))

k_vec=np.array([ i*(2*pi/L) for i in frange(-N/2,(N/2),1)])



uhat=np.fft.fftshift(np.fft.fft(u0))

for i in range(1,N_timesteps):    
    uhat_new = uhat/(1+D*i*dt*k_vec*k_vec)
    uhat = uhat_new    
    un=np.fft.ifft(np.fft.ifftshift(uhat))
    M[i+1,:]=un
    