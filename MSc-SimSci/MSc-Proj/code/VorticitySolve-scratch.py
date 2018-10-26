
# Parameters:
# Kill off energy piling up at large scales:
# take nu approx 0.1 omega_rms
nu_damping=0.05;

Nx=256;     # resolution in x
Ny=256;     # resolution in y
dt=1e-3;
Tf=100.0;     # final time
n_subdiv=5/dt;

im=1j;
dx=2*pi/Nx;
dy=2*pi/Ny;

nu=5.9e-30;
nuP=nu;

A0=1;
R=0.9;
kmax=9;
kmin=7;

# Crank-Nicholson parameter:
ck=0.5;
# Hyperdiffusivity
p=8;

# init='from_zero';
init='from_files';

xx=numpy.array([x*dx for x in range(0,Nx)]);
yy=numpy.array([y*dy for y in range(0,Ny)]);

t=0;

x01=pi-(pi/4);
x02=pi+(pi/4);

y01=pi;
y02=pi;

x03=x02;
y03=y02+(pi/(2*sqrt(2)));

r0=1/pi;
r0_sq=r0*r0;

(i,j) = np.meshgrid(range(1,Nx+1),range(1,Ny+1))
w=exp(-((i*dx-x01)**2+(j*dy-y01)**2)/r0_sq)+exp(-((i*dx-x02)**2+(j*dy-y02)**2)/r0_sq)-0.5*exp(-((i*dx-x03)*2+(j*dy-y03)*2)/r0_sq);
w=(1/(pi*r0_sq))*w;

theta=0*w;
for ii in range(0,Ny):
    for jj in range(0,Nx):
        rad=sqrt((xx[ii]-pi)**2+(yy[jj]-pi)**2);
        theta[jj,ii]=0.75*exp(-rad*rad/1.0);


kx = np.matrix(np.transpose(np.ones((1,Nx)))*(mod(np.array(range(1,Nx+1))-ceil(Nx/2+1),Nx)-floor(Nx/2))) # matrix of wavenumbers in x direction 
ky=(mod(np.transpose(np.matrix(range(0,Ny)))-ceil(Ny/2+1),Ny)-floor(Ny/2))*np.ones(Ny) # matrix of wavenumbers in y direction 

k_mod=sqrt(np.square(kx)+np.square(ky))
Mx_scale=(k_mod<=kmax)&(k_mod>=kmin)

forcing_hat=Mx_scale*np.ones((Ny,Nx));
A0_num=A0*Nx*Ny/sqrt(pi*(kmax*kmax-kmin*kmin))*(sqrt(2)/1);

kx=im*kx
ky=im*ky


ksquare_viscous=np.square(kx)+np.square(ky)        # Laplacian in Fourier space
k_hyperdiff=-((-1)**p)*ksquare_viscous**p
ksquare_poisson=ksquare_viscous    
ksquare_poisson[1,1]=1             # fixed Laplacian in Fourier space for Poisson's equation


dealias=np.any([abs(kx)<(2/3)*(Nx/2),abs(ky)<(2/3)*(Ny/2)], axis=0);
w_hat=np.fft.fftshift(np.fft.fft(w))
theta_hat=np.fft.fftshift(np.fft.fft(theta))

nn=((Tf/dt)/n_subdiv)+1;
t_vec=np.zeros(21);
l2_f_vec=t_vec;
l2_w_vec=t_vec;

ctr1=1;
ctr2=1;

x0=pi;
y0=x0;

while t<Tf:
    psi_hat = -w_hat/ksquare_poisson;  # Solve Poisson's Equation
    u  =real(np.fft.fft( ky*psi_hat));      # Compute  y derivative of stream function ==> u
    v  =real(np.fft.fft(-kx*psi_hat));      # Compute -x derivative of stream function ==> v
    w_x=real(np.fft.fft( kx*w_hat  ));      # Compute  x derivative of vorticity
    w_y=real(np.fft.fft( ky*w_hat  ));      # Compute  y derivative of vorticity
    theta_x=real(np.fft.fft( kx*theta_hat  ));      # Compute  x derivative of vorticity
    theta_y=real(np.fft.fft( ky*theta_hat  ));      # Compute  y derivative of vorticity

    conv     = u*w_x + v*w_y;         # evaluate the convective derivative (u,v).grad(w)   
    conv_hat = np.fft.fft(conv);              # go back to Fourier space

