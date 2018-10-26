#import random as ran
import numpy as np
import matplotlib.pyplot as plt


T_final=1000
Lx=1
Ly=1
nu_p=5.9e-30
nu_0=0.05
p=8
im=1j

Nx=256
dx=Lx/float(Nx)
x=np.array([v*dx for v in range(0,Nx)])

Ny=256
dy=Ly/float(Ny)
y=np.array([v*dx for v in range(0,Ny)])

(X,Y) = np.meshgrid(x,y)

dt=1e-3
N_timesteps=np.floor(T_final/dt)

kx_vec=[v*2*np.pi/Lx for v in range(-Nx/2,Nx/2)]
ky_vec=[v*2*np.pi/Ly for v in range(-Ny/2,Ny/2)]


Ksq=np.zeros((Nx,Ny))
Kx=np.zeros((Nx,Ny))
Ky=np.zeros((Nx,Ny))

for i in range(0,Nx):
    kx_val=kx_vec[i]
    for j in range(0,Ny):
        ky_val=ky_vec[j]    
        Ksq[i,j]=(kx_val**2)+(ky_val**2)
        Kx[i,j]=kx_val
        Ky[i,j]=ky_val
        

Ksq_laplace=Ksq
Ksq_laplace[(Nx/2)+1,(Ny/2)+1]=1

A0=1
R=0.9

kmax=9*(2*np.pi/Lx)
kmin=7*(2*np.pi/Lx)

Mx_scale=np.zeros((Nx,Ny))

for i in range(0,Nx):
    for j in range(0,Ny):
        k_radius=np.sqrt(kx_vec[i]**2+ky_vec[j]**2)
        if k_radius>kmin  and k_radius<kmax :
            Mx_scale[i,j]=1
        

A0_num=A0*Nx*Ny/np.sqrt(np.pi*(kmax*kmax-kmin*kmin))*(np.sqrt(2)/1)

omega0=np.zeros((Nx,Ny))

dealias=np.ones((Nx,Ny))

abs_kx_max=kx_vec[Nx-1]
abs_ky_max=ky_vec[Ny-1]

for i in range(0,Nx):
    for j in range(0,Ny):        
        abs_kx=abs(kx_vec[i])
        abs_ky=abs(ky_vec[j])
        if ((abs_kx>((2/float(3))*abs_kx_max)) and (abs_ky>((2/float(3))*abs_ky_max)) ):
            dealias[i,j]=0

omega=omega0
forcing_hat=np.zeros((Nx,Ny))
omega_arr=[omega]
np.random.seed(0)
for t_ctr in range(1,50):
    omega_hat=np.fft.fftshift(np.fft.fft2(omega))

    psi_hat=np.zeros((Nx,Ny),'complex')
    for i in range(0,Nx):
        for j in range(0,Ny):        
            if Ksq_laplace[i,j] > 0:
                psi_hat[i,j]=omega_hat[i,j]/complex(Ksq_laplace[i,j]) # make sure not to divide by zero here.
    
    dx_psi_hat=im*Kx*psi_hat
    dy_psi_hat=im*Ky*psi_hat
    
    dx_psi=np.real(np.fft.ifft2(np.fft.ifftshift(dx_psi_hat)))
    dy_psi=np.real(np.fft.ifft2(np.fft.ifftshift(dy_psi_hat)))

    u=dy_psi
    v=-dx_psi
    
    dx_omega_hat=im*Kx*omega_hat
    dy_omega_hat=im*Ky*omega_hat

    dx_omega=np.real(np.fft.ifft2(np.fft.ifftshift(dx_omega_hat)))
    dy_omega=np.real(np.fft.ifft2(np.fft.ifftshift(dy_omega_hat)))
    
    u_dot_grad_omega=u*dx_omega+v*dy_omega

    Mx_random=np.random.random((Nx,Ny))
    random_part_hat=Mx_scale*np.exp(2*np.pi*im*Mx_random)
    forcing_hat=np.sqrt(1-R*R)*A0_num*random_part_hat+R*forcing_hat

    u_dot_grad_omega_hat=np.fft.fftshift(np.fft.fft2(u_dot_grad_omega))

    numerator=omega_hat*(1-0.5*dt*nu_p*(Ksq**p))-dt*u_dot_grad_omega_hat+dt*(forcing_hat-nu_0*omega_hat)
    denominator=1+0.5*dt*nu_p*(Ksq**p)

    omega_hat_new=dealias*(numerator/denominator)

    omega=np.real(np.fft.ifft2(np.fft.ifftshift(omega_hat_new)))
    omega_arr.append(omega)
    print 'iter  ' + str(t_ctr) + ' ' + str(sum(omega))
    



#fig = plt.figure(figsize=(10, 10))
#ax = fig.add_subplot(111)
#ax.set_title('<<figure initialised>>')
#im_plt = plt.imshow(omega, extent=(0,Lx,0,Ly))
#im_plt
#im_plt.set_data(omega)
#plt.draw()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
for i in range(0,len(omega_arr)):
    im = plt.imshow(omega_arr[i], animated=True)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)
#ani.save('dynamic_images.mp4')

plt.show()

---------


fig = plt.figure()

# ims is a list of lists, each row is a list of artists to draw in the
# current frame; here we are just animating one artist, the image, in
# each frame
ims = []
for i in range(0,len(omega_arr)):
    im = plt.contourf(omega_arr[i], animated=True)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=300, blit=True, repeat_delay=1000)
#ani.save('dynamic_images.mp4')

#plt.show()