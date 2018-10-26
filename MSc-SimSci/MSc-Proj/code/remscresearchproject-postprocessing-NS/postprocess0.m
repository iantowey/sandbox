function [d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y]=postprocess0(w,theta)

global Nx
global Ny

[nNx,nNy]=size(w);
Nx=nNx;
Ny=nNy;

L=2*pi;
im=sqrt(-1);

kx=mod((1:Nx)-ceil(Nx/2+1),Nx)-floor(Nx/2); 
ky=mod((1:Ny)-ceil(Ny/2+1),Ny)-floor(Ny/2);

Kx=zeros(Nx,Ny);
Ky=zeros(Nx,Ny);

for i=1:Nx
    for j=1:Ny
        Kx(i,j)=kx(i);
        Ky(i,j)=ky(j);
    end
end

Kx=(2*pi/L)*im*Kx;
Ky=(2*pi/L)*im*Ky;

ksquare_viscous=Kx.^2+Ky.^2;        % Laplacian in Fourier space
ksquare_poisson=ksquare_viscous;    
ksquare_poisson(1,1)=1;               % fixed Laplacian in Fourier space for Poisson's equation

w_hat=fft2(w);
theta_hat=fft2(theta);

psi_hat = -w_hat./ksquare_poisson;

s_hat=Kx.*Ky.*psi_hat;
d_hat=Ky.*Ky.*psi_hat-Kx.*Kx.*psi_hat;
d_hat=d_hat/2;

thetax_hat=Kx.*theta_hat;
thetay_hat=Ky.*theta_hat;

d=real(ifft2(d_hat));
s=real(ifft2(s_hat));

u  =real(ifft2( Ky.*psi_hat));
v  =real(ifft2(-Kx.*psi_hat));

theta_x=real(ifft2(thetax_hat));
theta_y=real(ifft2(thetay_hat));

grad=sqrt(theta_x.^2+theta_y.^2);

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

beta_vec=zeros(Ny,Nx);

for i=1:Ny
    for j=1:Nx

        %******************************************************************
        
        nrm_theta=sqrt(theta_x(i,j)^2+theta_y(i,j)^2);

        if(nrm_theta==0)
            beta_vec(i,j)=0;
        else
            cosbeta=theta_x(i,j)/nrm_theta;
            sinbeta=theta_y(i,j)/nrm_theta;
                   
            if(cosbeta>=0)
                if(sinbeta>=0)
                    beta_val=acos(cosbeta);
                else
                    beta_val=2*pi-acos(cosbeta);
                end
            elseif(sinbeta>=0)
                    beta_val=pi-acos(abs(cosbeta));
            else
                beta_val=pi+acos(abs(cosbeta));
            end

            
            beta_vec(i,j)=beta_val;

        end
        
        %******************************************************************

    end
end

psi_angle=(pi/4)-phi;


psi_angle_hat=fft2(psi_angle);
psi_angle_x=real(ifft2(Kx.*psi_angle_hat));
psi_angle_y=real(ifft2(Ky.*psi_angle_hat));
conv=u.*psi_angle_x+v.*psi_angle_y;
w_tot=(w/2)+conv;

lambda_unsigned=sign(d).*sqrt(d.^2+s.^2);

sig1_w=sum(sum(w_tot))/(Nx*Ny);
sig2_w=sum(sum(w_tot.^2))/(Nx*Ny);

sig1_l=sum(sum(lambda_unsigned))/(Nx*Ny)
sig2_l=sum(sum(lambda_unsigned.^2))/(Nx*Ny)

sig_lw=sum(sum((w_tot-sig1_w).*(lambda_unsigned-sig1_l)))/(Nx*Ny);

Xangle=2*(psi_angle+beta_vec);

grad_av=sqrt(sum(sum(grad.^2))/(Nx*Ny));
ctr=1;
for i=1:Nx
    for j=1:Ny
        if(grad(i,j)>3*grad_av)
            ll(ctr)=-2*lambda_unsigned(i,j)*sin(Xangle(i,j));
            ctr=ctr+1;
        end
    end
end

prod=s.*(theta_y.^2-theta_x.^2)-2*d.*theta_x.*theta_y;
prod=prod./(grad.^2);


end