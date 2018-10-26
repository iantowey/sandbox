
T_final=100;

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

im=sqrt(-1);

% *************************************************************************
% Numerical parameters:

Nx=256;
dx=Lx/Nx;
x=0:dx:(Lx-dx);

Ny=256;
dy=Ly/Ny;
y=0:dy:(Ly-dy);

[X,Y]=meshgrid(x,y);

dt=1e-3;
N_timesteps=floor(T_final/dt);

% *************************************************************************
% Allocating some arrays:

% Set up wavevectors:
kx_vec=(2*pi/Lx)*( -(Nx/2):1:(Nx/2)-1 );
ky_vec=(2*pi/Ly)*( -(Ny/2):1:(Ny/2)-1 );

Ksq=zeros(Nx,Ny);
Kx=zeros(Nx,Ny);
Ky=zeros(Nx,Ny);

for i=1:Nx
    kx_val=kx_vec(i);
    for j=1:Ny
        ky_val=ky_vec(j);
        
        Ksq(i,j)=(kx_val^2)+(ky_val^2);
        Kx(i,j)=kx_val;
        Ky(i,j)=ky_val;
    end
end

% Remark: If Ksq(i,j)==0 set Ksq(i,j) to 1.

Ksq_laplace=Ksq;
Ksq_laplace((Nx/2)+1,(Ny/2)+1)=1;

% *************************************************************************
% Forcing function.  The idea for this forcing is taken from the paper of
% Lilly (1969).

A0=1;
R=0.9;

kmax=9*(2*pi/Lx);
kmin=7*(2*pi/Lx);

% The following lines makes the array forcing_hat zero, unless k_mod is
% within the specified range.

Mx_scale=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        
        k_radius=kx_vec(i)^2+ky_vec(j)^2;
        k_radius=sqrt(k_radius);
        
        if( (k_radius>kmin) && (k_radius<kmax) )
            Mx_scale(i,j)=1;
        end
    end
end


forcing_hat=zeros(Nx,Ny);
A0_num=A0*Nx*Ny/sqrt(pi*(kmax*kmax-kmin*kmin))*(sqrt(2)/1);

% *************************************************************************
% Initialization:

% Here, you will have to initialize omega0.  You can actually initialize it
% to zero.

omega0=zeros(Nx,Ny);

% *************************************************************************
    % De-aliasing matrix:

dealias=ones(Nx,Ny);

abs_kx_max=kx_vec(end);
abs_ky_max=ky_vec(end);

for i=1:Nx
    for j=1:Ny
        
        abs_kx=abs(kx_vec(i));
        abs_ky=abs(ky_vec(j));
        
        if( (abs_kx>(2/3)*abs_kx_max) && (abs_ky>(2/3)*abs_ky_max) )
            dealias(i,j)=0;
        end
    end
end

% *************************************************************************
% Time loop:

omega=omega0;
rand('twister', 0);

for t_ctr=1:50
    
   omega_hat=fftshift(fft2(omega));
   psi_hat=omega_hat./Ksq_laplace; % make sure not to divide by zero here.
   
   dx_psi_hat=sqrt(-1)*Kx.*psi_hat;
   dy_psi_hat=sqrt(-1)*Ky.*psi_hat;
   
   dx_psi=ifft2(ifftshift(dx_psi_hat),'symmetric');
   dy_psi=ifft2(ifftshift(dy_psi_hat),'symmetric');
   
   u=dy_psi;
   v=-dx_psi;
   
   dx_omega_hat=sqrt(-1)*Kx.*omega_hat;
   dy_omega_hat=sqrt(-1)*Ky.*omega_hat;
   
   dx_omega=ifft2(ifftshift(dx_omega_hat),'symmetric');
   dy_omega=ifft2(ifftshift(dy_omega_hat),'symmetric');
   
   u_dot_grad_omega=u.*dx_omega+v.*dy_omega;
   
   %***********************************************************************
   % Generate Q_hat
   % Q_hat=forcing_hat;
   Mx_random=rand(Nx,Ny);
   random_part_hat=Mx_scale.*exp(2*pi*im*Mx_random);        
   forcing_hat=sqrt(1-R*R)*A0_num*random_part_hat+R*forcing_hat;

   %***********************************************************************
   
   u_dot_grad_omega_hat=fftshift(fft2(u_dot_grad_omega));
   
   % Implement Equation (8) on the board (01/06):
   
   numerator=omega_hat.*(1-0.5*dt*nu_p*(Ksq.^p))-dt*u_dot_grad_omega_hat+dt*(forcing_hat-nu_0*omega_hat);
   denominator=1+0.5*dt*nu_p*(Ksq.^p);
   
   % omega_hat_new=numerator./denominator;
   num_div_den=numerator./denominator;
   omega_hat_new=dealias.*(numerator./denominator);
    
   % Updated value of omega:
   omega=ifft2(ifftshift(omega_hat_new),'symmetric');
   fprintf('iter %d %d %d\n',t_ctr,sum(sum(omega_hat_new)),sum(sum(omega)));
   
   if(mod(t_ctr,1)==0)
       set(pcolor(X,Y,omega),'edgecolor','none')
       axis equal
       xlim([0 Lx])
       ylim([0 Ly])
       colorbar
       title(strcat('vorticity at t=',num2str(t_ctr*dt)))
       drawnow
   end
end


