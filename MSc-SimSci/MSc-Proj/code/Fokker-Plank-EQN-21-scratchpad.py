import numpy as np
import json

Tf = 0.1
w = 0.5
tau=0.01
delta=1
k_corr = 0.5
x = 0
tauy=tau
tauz=tau
Dy = 1
Dz = 1
gamma_val = 0.5
Lx=2*np.pi
Lz=75
Ly=Lz
Nx=64
Ny=32
Nz=32
dx=Lx/float(Nx)
dy=Ly/float(Ny)
dz=Lz/float(Nz)
Vmax=2*(Lz/2)/gamma_val
dx_min=min([dx,dy,dz])
dt=0.1*dx_min/Vmax
im=1j
t_vec=np.arange(0,Tf,dt)
n_timesteps=len(t_vec)
ck=0.95
x=np.arange(-Lx/float(2),(Lx/float(2)),dx)
y=np.arange(-Ly/float(2),(Ly/float(2)),dy)
z=np.arange(-Lz/float(2),(Lz/float(2)),dz)

p0=np.zeros((Nx,Ny,Nz))

w_y=2*Dy/tauy
w_z=2*Dz*(1-k_corr*k_corr)/tauz

for i in range(0,Nx):
    for j in range(0,Ny):
        for k in range(0,Nz):
            p0[i,j,k] = (1/float(Lx)) * (1-0.5*np.cos(2*x[i])) * np.exp(-(y[j]**2)/w_y) * np.exp(-(z[k]**2)/w_z)


nrm=sum(p0)*dx*dy*dz
p0=p0/nrm

kx=(([i+1 for i in range(0,Nx)] -ceil(Nx/2+1)) % Nx)-floor(Nx/2)
ky=(([i+1 for i in range(0,Ny)] -ceil(Ny/2+1)) % Ny)-floor(Ny/2)
kz=(([i+1 for i in range(0,Nz)] -ceil(Nz/2+1)) % Nz)-floor(Nz/2)

Kx=np.zeros(shape=(Nx,Ny,Nz))
Ky=np.zeros(shape=(Nx,Ny,Nz))
Kz=np.zeros(shape=(Nx,Ny,Nz))

for i in range(0,Nx):
    for j in range(0,Ny):
        for k in range(0,Nz):
            Kx[i,j,k]=kx[i]
            Ky[i,j,k]=ky[j]
            Kz[i,j,k]=kz[k]

Y=np.zeros(shape=(Nx,Ny,Nz))
Z=np.zeros(shape=(Nx,Ny,Nz))

for i in range(0,Nx):
    for j in range(0,Ny):
        for k in range(0,Nz):
            Y[i,j,k]=y[j]
            Z[i,j,k]=z[k]


Kx=(2*np.pi/float(Lx))*Kx
Ky=(2*np.pi/float(Ly))*Ky
Kz=(2*np.pi/float(Lz))*Kz

Kx=im*Kx
Ky=im*Ky
Kz=im*Kz

small=0
ksquare_lap=small*(Kx**2)+(Dy/(tauy**2))*(Ky**2)+(Dz*(1-k_corr*k_corr)/(tauz**2))*(Kz**2)

Vx=np.zeros(shape=(Nx,Ny,Nz))

for i in range(0,Nx):
    for j in range(0,Ny):
        for k in range(0,Nz):
            x_val=x[i]
            y_val=y[j]
            z_val=z[k] 
            prefac=1     
            Vx[i,j,k]=prefac*((w/gamma_val)+(1/gamma_val)*(-np.cos(x_val)+k_corr*np.sqrt(delta))*y_val+(z_val/gamma_val))

n_subdiv = 50

residual_vec=np.zeros(n_timesteps)
av_vec=residual_vec

p_hat_old=np.fft.fftn(p0)
Vp_hat=np.fft.fftn(Vx*p0)
py_hat=np.fft.fftn(p0*Y)
pz_hat=np.fft.fftn(p0*Z)
dy_py_hat= Ky*py_hat
dz_pz_hat= Kz*pz_hat

src0_hat=-Kx*Vp_hat
src1_hat=(1/tauy)*dy_py_hat+(1/tauz)*dz_pz_hat
src_hat=src0_hat+src1_hat

p_hat=((1+(1-ck)*dt*ksquare_lap)/(1-ck*dt*ksquare_lap))*p_hat_old+dt*(src_hat/(1-ck*dt*ksquare_lap))
p=np.real(np.fft.ifftn(p_hat))

residual_val=(abs(p-p0)).max()
residual_vec[0]=residual_val

p_x=0*x
for i in range(0,len(x)):
    p_x[i]=sum(p[i,:,:])

x1=np.concatenate([x, [np.pi]])
p1=np.concatenate([p_x,[p_x[0]]])
av_vec[0]=sum(p1*x1)/sum(p1)


for t_ctr in range(1,len(t_vec)):
    print t_ctr
    Vp_hat=np.fft.fftn(Vx*p)
    py_hat=np.fft.fftn(p*Y)
    pz_hat=np.fft.fftn(p*Z)
    dy_py_hat= Ky*py_hat
    dz_pz_hat= Kz*pz_hat
    
    src0_hat=-Kx*Vp_hat
    src1_hat=(1/tauy)*dy_py_hat+(1/tauz)*dz_pz_hat
    
    src_hat=src0_hat+src1_hat

    p_hat_new=(1/(1-2*dt*ksquare_lap))*p_hat_old+2*dt*(src_hat/(1-2*dt*ksquare_lap))
    
    p_hat_old=p_hat
    p_old=np.real(np.fft.ifftn(p_hat_old))
    
    p_hat=p_hat_new    
    p=np.real(np.fft.ifftn(p_hat))

    residual_val=(abs(p_old-p)).max()
    residual_vec[t_ctr]=residual_val
    for i in range(1,len(x)):
        p_x[i]=sum(p[i,:,:])

    x1=np.concatenate([x, [np.pi]])
    p1=np.concatenate([p_x,[p_x[0]]])
    av_vec[t_ctr]=sum(p1*x1)/sum(p1)
    
    #if t_ctr % n_subdiv == 0:
    #    p_x=p_x=0*x
    #    
    #    for i in range(0, len(x)):
    #        p_x[i]=sum(p[i,:,:])*dy*dz
    #    
    #    nrm1=sum(p_x)*dx;
    #    
    #    x_aug=[x,pi]
    #    p_aug=[p_x,p_x(0)]
    #    
    #    p_y=y*0
    #    
    #    for i in range(0, len(y)):
    #        p_y[i]=sum(p[:,i,:])*dx*dz
    #    
    #    nrm2=sum(p_y)*dy
    #    
    #    p_z=z*0
    #    
    #    for i in range(0, len(z)):
    #        p_z[i]=sum(p[:,:,i])*dx*dy
    #    
    #    nrm3=sum(p_z)*dz


nrm=sum(p)*dx*dy*dz
p_joint=p/nrm

#**************************************************************************
# Compute marginal distributions

p_x=np.zeros(len(x))
p_y=np.zeros(len(y))
p_z=np.zeros(len(z))

for i in range(0,len(x)):
    p_x[i]=sum(p_joint[i,:,:])*dy*dz

plt.plot(x,p_x) #PDF of X


for j in range(0,len(y)):
    p_y[j]=sum(p_joint[:,j,:])*dx*dz

import matplotlib.pyplot as plt
plt.plot(x,p_x) #PDF of X

for k in range(0,len(z)):
    p_z[k]=sum(p_joint[:,:,k])*dy*dz

p_xy=np.zeros((len(x),len(y)))
for i in range(0,len(x)):
    for j in range(0,len(y)):
        p_xy[i,j]=sum(p_joint[i,j,:])*dz
plt.imshow(p_xy)

p_yz=np.zeros((len(y),len(z)))
for j in range(0,len(y)):
    for k in range(0,len(z)):
        p_yz[j,k]=sum(p_joint[:,j,k])*dx
    
nrm=sum(p_yz)*dy*dz
p_yz=p_yz/nrm

p_yz_anal=np.zeros((len(y),len(z)))

w_y=2*Dy/tauy
w_z=2*Dz*(1-k_corr*k_corr)/tauz

for j in range(0,len(y)):
    for k in range(0,len(z)):
        p_yz_anal[j,k]=np.exp(-y[j]**2/w_y)*np.exp(-z[k]**2/w_z)
    

nrm=sum(p_yz_anal)*dy*dz
p_yz_anal=p_yz_anal/nrm


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
XX,YY = np.meshgrid(x,y)
XX = XX
YY = YY
ZZ = p_xy.T
colors = cm.viridis(ZZ)
rcount, ccount, _ = colors.shape
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(XX, YY, ZZ, rcount=rcount, ccount=ccount,
                       facecolors=colors, shade=False)
ax.set_xlim(min(x), max(x))
ax.set_ylim(min(y), max(y))
#ax.set_zlim(min(z), max(z))
ax.set_xlabel('x',fontsize=20)
ax.set_ylabel('y', fontsize=20)
ax.set_zlabel('z', fontsize=20)
fig.colorbar(surf, shrink=0.5, aspect=5)

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
YY,ZZ = np.meshgrid(y,z)
P_YZ = p_yz.T
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(YY, ZZ, P_YZ, rstride=1, cstride=1)
plt.show()

#**************************************************************************
# Compute lambda distribution

lambda_range=np.arange(-30,30,0.01)
p_lambda=0*lambda_range

for i in range(0,len(lambda_range)):
    sum_val=0
    lambda_val=lambda_range[i]
    
    if lambda_val==0:
        sum_val=0
    else:
        for ii in range(0,len(x)):
            g_val=-2*np.sin(x[ii])
    
            if((abs(lambda_val)/(abs(g_val)+1e-8))>=(Ly/2)):
                summand=0
            else:
                ll_map=1+(1/dy)*((lambda_val/g_val)+(Ly/2))
                ll_map=np.floor(ll_map)-1
                
                if(abs(ll_map)>Ny):
                    summand=0
                else:
                    summand=p_xy[ii,ll_map]*(dx/abs(g_val))
            sum_val=sum_val+summand
    p_lambda[i]=sum_val


ix=np.argmin(abs(lambda_range))

p_lambda[ix]=(p_lambda[ix-1]+p_lambda[ix+1])/2

plt.plot(lambda_range,p_lambda)

np.savetxt("/home/ian/Desktop/MSc Proj/fokker-plank/p.csv", p, delimiter=",")

with open('/home/ian/Desktop/MSc Proj/fokker-plank/p.csv', 'w') as f: 
    f.write(json.dumps(p, default=lambda x: list(x), indent=4))

with open('/home/ian/Desktop/MSc Proj/fokker-plank/x.csv', 'w') as f: 
    f.write(json.dumps(x, default=lambda x: list(x), indent=4))
with open('/home/ian/Desktop/MSc Proj/fokker-plank/y.csv', 'w') as f: 
    f.write(json.dumps(y, default=lambda x: list(x), indent=4))
with open('/home/ian/Desktop/MSc Proj/fokker-plank/z.csv', 'w') as f: 
    f.write(json.dumps(z, default=lambda x: list(x), indent=4))

with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_x.csv', 'w') as f: 
    f.write(json.dumps(p_x, default=lambda x: list(x), indent=4))
with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_y.csv', 'w') as f: 
    f.write(json.dumps(p_y, default=lambda x: list(x), indent=4))
with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_z.csv', 'w') as f: 
    f.write(json.dumps(p_z, default=lambda x: list(x), indent=4))

with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_xy.csv', 'w') as f: 
    f.write(json.dumps(p_xy, default=lambda x: list(x), indent=4))
with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_yz.csv', 'w') as f: 
    f.write(json.dumps(p_yz, default=lambda x: list(x), indent=4))

with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_yz_anal.csv', 'w') as f: 
    f.write(json.dumps(p_yz_anal, default=lambda x: list(x), indent=4))

with open('/home/ian/Desktop/MSc Proj/fokker-plank/p_lambda.csv', 'w') as f: 
    f.write(json.dumps(p_lambda, default=lambda x: list(x), indent=4))


#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator, FormatStrFormatter
#import numpy as np
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#
#Y_, Z_ = np.meshgrid(y, z)
#P_YZ = p_yz
#
## Plot the surface.
#surf = ax.plot_surface(Y_, Z_, P_YZ, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#
## Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#
## Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
#
#plt.show()

