function []=new_vorticity_eqn(T_final)

% *************************************************************************
% Lennon O'Naraigh
% Matlab code snipipet to help the MSc project students to solve the 2D 
% vorticity equation.
%
% *************************************************************************
% Physical parameters:

Lx=1;
Ly=1;
nu_p=5.9e-30;
nu_0=0.05;

p=8;

% *************************************************************************
% Numerical parameters:

Nx=256;
dx=Lx/Nx;
x=0:dx:(Lx-dx);

Ny=256;
dy=Ly/Ny;
y=0:dy:(Ly-dy);

dt=0.01;
N_timesteps=floor(T_final/dt);

% *************************************************************************
% Parameters to do with the Q-function:

A0=1;
R=0.9;
kmax=9;
kmin=7;

% *************************************************************************
% Allocating some arrays:

% Initial data:
omega0=zeros(1,N);

% Set up wavevector:
kx_vec=(2*pi/Ly)*( -(Nx/2):1:(Nx/2)-1 );

% Set up wavevector:
ky_vec=(2*pi/Ly)*( -(Ny/2):1:(Ny/2)-1 );

Ksq=zeros(Nx,Ny);
Kx=zeros(Nx,Ny);
Ky=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        Ksq(i,j)=kx_vec(i)^2+ky_vec(j)^2;
        Kx(i,j)=kx_vec(i);
        Ky(i,j)=ky_vec(j);
    end
end

% Remark: If Ksq(i,j)==0 set Ksq(i,j) to 1.


% *************************************************************************
% Initialization:

% Here, you will have to initialize omega0.  You can actually initialize it
% to zero.

omega0=zeros(Nx,Ny);

% Remark: Make sure that omega0 has mean zero.

% *************************************************************************
% Time loop:


for t_ctr=1:N_timesteps
    
   omega_hat=fftshift(fft2(omega));
   
   psi_hat=omega_hat./Ksq; % make sure not to divide by zero here.
   
   dx_psi_hat=sqrt(-1)*Kx.*psi_hat;
   dy_psi_hat=sqrt(-1)*Ky.*psi_hat;
   
   dx_psi=ifft2(ifftshift(dx_psi_hat));
   dy_psi=ifft2(ifftshift(dy_psi_hat));
   
   u=dy_psi;
   v=-dx_psi;
   
   dx_omega_hat=sqrt(-1)*Kx.*omega_hat;
   dy_omega_hat=sqrt(-1)*Ky.*omega_hat;
   
   dx_omega=ifft2(ifftshift(dx_omega_hat));
   dy_omega=ifft2(ifftshift(dy_omega_hat));
   
   u_dot_grad_omega=u.*dx_omega+v.*dy_omega;
   
   % TO DO: WRITE A FUNCTION TO DEFINE Q
   %Q_hat=get_Q_hat(t_ctr);

    %generate Q_hat
    k_mod=sqrt(kx.^2+ky.^2);
    Mx_scale=(k_mod<=kmax)&(k_mod>=kmin);

    forcing_hat=Mx_scale.*(ones(Ny,Nx));
    A0_num=A0*Nx*Ny/sqrt(pi*(kmax*kmax-kmin*kmin))*(sqrt(2)/1);

    Mx_random=rand(Ny,Nx);
    random_part_hat=Mx_scale.*exp(2*pi*im*Mx_random);        
    forcing_hat=sqrt(1-R*R)*A0_num*random_part_hat+R*forcing_hat;
	
    Q_hat = forcing_hat;
   
   u_dot_grad_omega_hat=fftshift(fft2(u_dot_grad_omega));
   
   % Implement Equation (8) on the board (01/06):
   
   numerator=omega_hat-dt*u_dot_grad_omega_hat-0.5*dt*nup*(Ksq.^p).*omega_hat+dt*(Q_hat-nu0*omega_hat);
   denominator=1+0.5*dt*nup*(Ksq.^p);
   
   omega_hat_new=numerator./denominator;
   
   % TO DO: DE-ALIAS
   % omega_hat_new=dealias.*(numerator./denominator);
   
   % Updated value of omega:
   omega=ifft2(ifftshift(omega_hat_new));
   
   
   
   
   
   
end

end

