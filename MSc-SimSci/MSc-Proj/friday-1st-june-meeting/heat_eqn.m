function [x,M]=heat_eqn(T_final)

% *************************************************************************
% Lennon O'Naraigh and Zhou Lyu, 01/06/2018
% Matlab code to solve the 1D heat equation with periodic boundary
% conditions.
% 
% u_t=Du_{xx},   x\in (0,L), t>0
% u(x,t=0)=u0
%
% *************************************************************************
% Physical parameters:

T_final=10;

L=1;
D=1;

% *************************************************************************
% Numerical parameters:

N=64;
dx=L/N;sa
x=0:dx:(L-dx);

dt=0.01;
N_timesteps=floor(T_final/dt);

% *************************************************************************
% Allocating some arrays:

% Initial data:
u0=zeros(1,N);

% Array of spacetime data:
M=zeros(N_timesteps+1,N);

% Set up wavevector:
k_vec=(2*pi/L)*( -(N/2):1:(N/2)-1 );

% *************************************************************************
% Initialization:

for i=1:N
    
    if( (x(i)<1/4) || (x(i)>3/4))
        u0(i)=0;
    else
        u0(i)=1;
    end
end

M(1,:)=u0;
 
% Compute the Fourier Transform of u:
uhat=fftshift(fft(u0));

% *************************************************************************
% Time loop:

for t_ctr=1:N_timesteps
    
    uhat_new=uhat./(1+D*dt*k_vec.*k_vec);
    uhat=uhat_new;
    
    un=ifft(ifftshift(uhat));
    M(t_ctr+1,:)=un;
    
end



end

