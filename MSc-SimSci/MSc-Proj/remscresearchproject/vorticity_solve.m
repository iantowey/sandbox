% function [xx,yy,w_f,theta_f,t_vec,l2_f_vec,l2_w_vec]=vorticity_solve()
function [x0_vec,y0_vec,w_tot_vec,lambda_vec]=vorticity_solve()
clear all;

%**************************************************************************
%                                                                         
%  Dw/Dt = nu.Laplacian(w)+Q
%  
%  Laplacian(psi) = -w  
%  u = psi_y
%  v =-psi_x
%
%  Q_hat(n+1,k)=A0_num*sqrt(1-R*R)*random(n+1,k)+R*Q_hat(n,k)
%  
%**************************************************************************

global Nx
global Ny

% Parameters:
% Kill off energy piling up at large scales:
% take nu approx 0.1 omega_rms
nu_damping=0.05;
Nx=256;     
Ny=256;     
dt=1e-3;
Tf=100.0;     
n_subdiv=5/dt;
im=sqrt(-1);
dx=2*%pi/Nx;
dy=2*%pi/Ny;
nu=5.9e-30;
nuP=nu;
A0=1;
R=0.9;
kmax=9;
kmin=7;

% Crank-Nicholson parameter:
ck=0.5;
% Hyperdiffusivity
p=8;

% init='from_zero';
init='from_files';

%**************************************************************************
% Initialize code

xx=(0:(Nx-1))*dx;
yy=(0:(Ny-1))*dy;

t=0;

% Define initial vorticity distribution
switch lower(init)
    case {'from_zero'}
        disp('generating intial fields')
        x01=%pi-(%pi/4);
        x02=%pi+(%pi/4);

        y01=%pi;
        y02=%pi;

        x03=x02;
        y03=y02+(pi/(2*sqrt(2)));

        r0=1/pi;
        r0_sq=r0*r0;

        [i,j]=meshgrid(1:Nx,1:Ny);
        w=exp(-((i*dx-x01).^2+(j*dy-y01).^2)/r0_sq)+exp(-((i*dx-x02).^2+(j*dy-y02).^2)/r0_sq)-0.5*exp(-((i*dx-x03).^2+(j*dy-y03).^2)/r0_sq);
        w=(1/(pi*r0_sq))*w;
        
        theta=0*w;
        for ii=1:Ny
            for jj=1:Nx
                rad=sqrt((xx(ii)-pi)^2+(yy(jj)-pi)^2);
                theta(jj,ii)=0.75*exp(-rad*rad/1.0);
            end
        end
    case {'from_files'}
        struct=load('w_start.mat');
        w=struct.w;
        struct=load('theta_start.mat');
        theta=struct.theta;
    otherwise
        disp('unknown initialization');
        return
end
disp('initialization complete')


kx=ones(1,Ny)'*(mod((1:Nx)-ceil(Nx/2+1),Nx)-floor(Nx/2)); % matrix of wavenumbers in x direction 
ky=(mod((1:Ny)'-ceil(Ny/2+1),Ny)-floor(Ny/2))*ones(1,Nx); % matrix of wavenumbers in y direction 


%**************************************************************************
% Forcing

k_mod=sqrt(kx.^2+ky.^2);
Mx_scale=(k_mod<=kmax)&(k_mod>=kmin);

forcing_hat=Mx_scale.*(ones(Ny,Nx));
A0_num=A0*Nx*Ny/sqrt(pi*(kmax*kmax-kmin*kmin))*(sqrt(2)/1);

%**************************************************************************

kx=im*kx;
ky=im*ky;

ksquare_viscous=kx.^2+ky.^2;        % Laplacian in Fourier space
k_hyperdiff=-((-1)^p)*ksquare_viscous.^p;
ksquare_poisson=ksquare_viscous;    
ksquare_poisson(1,1)=1;             % fixed Laplacian in Fourier space for Poisson's equation

dealias=abs(kx)<(2/3)*(Nx/2)&abs(ky)<(2/3)*(Ny/2);

w_hat=fft2(w);
theta_hat=fft2(theta);

%**************************************************************************
% Update step

nn=((Tf/dt)/n_subdiv)+1;
t_vec=0*(1:nn);
l2_f_vec=t_vec;
l2_w_vec=t_vec;

ctr1=1;
ctr2=1;

x0=pi;
y0=x0;
while t<Tf
    % Compute the stream function and get the velocity and gradient of vorticity
    psi_hat = -w_hat./ksquare_poisson;  % Solve Poisson's Equation
    u  =real(ifft2( ky.*psi_hat));      % Compute  y derivative of stream function ==> u
    v  =real(ifft2(-kx.*psi_hat));      % Compute -x derivative of stream function ==> v
    w_x=real(ifft2( kx.*w_hat  ));      % Compute  x derivative of vorticity
    w_y=real(ifft2( ky.*w_hat  ));      % Compute  y derivative of vorticity
    theta_x=real(ifft2( kx.*theta_hat  ));      % Compute  x derivative of vorticity
    theta_y=real(ifft2( ky.*theta_hat  ));      % Compute  y derivative of vorticity

    conv     = u.*w_x + v.*w_y;         % evaluate the convective derivative (u,v).grad(w)   
    conv_hat = fft2(conv);              % go back to Fourier space
    
    conv_theta     = u.*theta_x + v.*theta_y;          
    conv_theta_hat = fft2(conv_theta);              

    conv_hat = dealias.*conv_hat;   % Perform spherical dealiasing 2/3 rule
    
    % No need to de-alias convective derivative of theta because it is
    % linear ???
    conv_theta_hat = dealias.*conv_theta_hat;

    Mx_random=rand(Ny,Nx);
    random_part_hat=Mx_scale.*exp(2*pi*im*Mx_random);        
    forcing_hat=sqrt(1-R*R)*A0_num*random_part_hat+R*forcing_hat;
    
    damping_hat=-nu_damping*w_hat;

    % Compute Solution at the next step
    w_hat_new = ((1 + (1-ck)*dt*nu*k_hyperdiff)./(1 - ck*dt*nu*k_hyperdiff)).*w_hat         + dt*((-conv_hat+forcing_hat+damping_hat )./(1 - ck*dt*nu*k_hyperdiff)); 
    theta_hat_new = ((1 + (1-ck)*dt*nuP*k_hyperdiff)./(1 - ck*dt*nuP*k_hyperdiff)).*theta_hat + dt*((-conv_theta_hat+0)./(1 - ck*dt*nuP*k_hyperdiff)); 

    t=t+dt;

    % Generating outputs
    if( mod(ctr1,n_subdiv)==0 )
        w=real(ifft2(w_hat_new));
        theta=real(ifft2(theta_hat_new));
        forcing=real(ifft2(forcing_hat));
        
        l2_nrm_f=(1/(Nx*Ny))*sum(sum(forcing.*forcing));
        l2_nrm_w=(1/(Nx*Ny))*sum(sum(w.*w));
        l2_f_vec(ctr2)=l2_nrm_f;
        l2_w_vec(ctr2)=l2_nrm_w;
        
        t_string=num2str(ctr2);
%         name_string_w     = 'output_vorticity';
%         name_string_theta = 'output_theta';
%         filename_w    = strcat(name_string_w,t_string);
%         filename_theta= strcat(name_string_theta,t_string);
%         save(filename_w,'w')
%         save(filename_theta,'theta')
        
        t_vec(ctr2)=t;
        ctr2=ctr2+1;
    end
        
        [x0_new,y0_new]=get_particle_states0(u,v,x0,y0,dt);
        x0_vec(ctr1)=x0_new;
        y0_vec(ctr1)=y0_new;
        
        s_hat=kx.*ky.*psi_hat;
        d_hat=ky.*ky.*psi_hat-kx.*kx.*psi_hat;
        d_hat=d_hat/2;

        d=real(ifft2(d_hat));
        s=real(ifft2(s_hat));
        
        [w_tot_x0,lambda_x0]=get_particle_states1(s,d,x0,y0);
        w_tot_vec(ctr1)=w_tot_x0;
        lambda_vec(ctr1)=lambda_x0;
        
        x0=x0_new;
        y0=y0_new;
    
        w_hat=w_hat_new;
        theta_hat=theta_hat_new;
        ctr1=ctr1+1;
    
        if(mod(ctr1,100)==0)
            w=real(ifft2(w_hat_new));
            max(max(w))
            t
        end
     
end

w_f=real(ifft2(w_hat_new));
theta_f=real(ifft2(theta_hat_new));

!pwd | mail -s finished LENNON.ONARAIGH@UCD.IE

    function [x0_new,y0_new]=get_particle_states0(u,v,x0,y0,dt)
       
        u0=interp2(u,x0,y0);
        v0=interp2(v,x0,y0);
        
        x0_new=x0+u0*dt;
        y0_new=y0+v0*dt;
    
    end


    function [w_tot_x0,lambda_x0]=get_particle_states1(s,d,x0,y0)
       
        d=d';
        s=s';
        
        phi=zeros(Nx,Ny);
        
        for i=1:Nx
            for j=1:Ny
                d_val=d(i,j);
                s_val=s(i,j);

                if(d_val==0)
                    phi(i,j)=0;
                else
                    alpha_val=s_val/d_val;
                    nrm=2*sqrt(alpha_val*alpha_val+1)*(-alpha_val+sqrt(alpha_val*alpha_val+1));
                if(nrm<1)
                    nrm=1;
                end
                    nrm=sqrt(nrm);
                    phi(i,j)=acos(1/nrm);
                end
            end
        end
        
        psi_angle=(pi/4)-phi';

        psi_angle_hat=fft2(psi_angle);
        psi_angle_x=real(ifft2(kx.*psi_angle_hat));
        psi_angle_y=real(ifft2(ky.*psi_angle_hat));
        conv=u.*psi_angle_x+v.*psi_angle_y;
        w_tot=(w/2)+conv;
        
        w_tot_x0=interp2(w_tot,x0,y0);
        
        
        s=s';
        d=d';
        lambda_unsigned=sign(d).*sqrt(d.^2+s.^2);
        lambda_x0=interp2(lambda_unsigned,x0,y0);
       
    end

end

