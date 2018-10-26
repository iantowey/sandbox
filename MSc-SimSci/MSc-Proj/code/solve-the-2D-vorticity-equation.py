import numpy as np

def frange(start, stop, step):
    i = start
    while i < stop:
        yield i
        i += step
        

#MSc project students to solve the 2D vorticity equation.

T_final = 10

# *************************************************************************
# Physical parameters:

Lx=1
Ly=1
nu_p=5.9e-30
nu_0=0.05

p=8

# *************************************************************************
# Numerical parameters:

Nx=256
dx=Lx/float(Nx)
x=np.array([i for i in frange(0,Lx,dx)])

Ny=256
dy=Ly/float(Ny)
y=np.array([i for i in frange(0,Ly,dy)])

dt=0.01
N_timesteps=int(T_final/dt)

# *************************************************************************
# Parameters to do with the Q-function:
A0=1
R=0.9
kmax=9
kmin=7

omega0=np.zeros((1,Nx))

kx_vec=np.array([ i*(2*pi/Lx) for i in frange(-Nx/2,(Nx/2),1)])
ky_vec=np.array([ i*(2*pi/Ly) for i in frange(-Ny/2,(Ny/2),1)])

Ksq=np.zeros((Nx,Ny))
Kx=np.zeros((Nx,Ny))
Ky=np.zeros((Nx,Ny))

for i in range(1,Nx):
    for j in range(1,Ny):
        Ksq[i,j]=kx_vec[i]**2+ky_vec[j]**2
        Kx[i,j]=kx_vec[i]
        Ky[i,j]=ky_vec[j]

omega0=np.zeros((Nx,Ny))

def get_Q(t):
    

for t_ctr in range(1,N_timesteps):    
   omega_hat=np.fft.fftshift(np.fft.fft(omega0))
   
   psi_hat=omega_hat/Ksq # make sure not to divide by zero here.
   
   dx_psi_hat=sqrt(-1)*Kx*psi_hat
   dy_psi_hat=sqrt(-1)*Ky*psi_hat
   
   dx_psi=ifft2(ifftshift(dx_psi_hat))
   dy_psi=ifft2(ifftshift(dy_psi_hat))
   
   u=dy_psi
   v=-dx_psi
   
   dx_omega_hat=sqrt(-1)*Kx*omega_hat
   dy_omega_hat=sqrt(-1)*Ky*omega_hat
   
   dx_omega=np.fft.ifft(np.fft.ifftshift(dx_omega_hat))
   dy_omega=np.fft.ifft(np.fft.ifftshift(dy_omega_hat))
   
   u_dot_grad_omega=u*dx_omega+v*dy_omega
   
   #*****************************************************
   # TO DO: WRITE A FUNCTION TO DEFINE Q
   Q=get_Q(t_ctr)
   Q_hat=np.fft.fftshift(np.fft.fft(Q))
   
   u_dot_grad_omega_hat=np.fft.fftshift(np.fft.fft(u_dot_grad_omega))
   
   #*****************************************************
   # Implement Equation (8) on the board (01/06):
   
   numerator=omega_hat-dt*u_dot_grad_omega_hat-0.5*dt*nup*(Ksq**p)*omega_hat+dt*(Q_hat-nu0*omega_hat)
   denominator=1+0.5*dt*nup*(Ksq**p)
   
   omega_hat_new=numerator/denominator
   
   # TO DO: DE-ALIAS
   # omega_hat_new=dealias*(numerator/denominator)
   
   # Updated value of omega:
   omega=np.fft.ifft(np.fft.ifftshift(omega_hat_new))
