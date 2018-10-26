function [x,y,z,p_joint,p_x,p_y,p_z,p_yz,p_yz_anal,residual_vec,av_vec,lambda_range,p_lambda]=jointpde_solve3d_param(Tf,w,delta,k_corr,tau)

Tf=20;
w=0.5;
tau=0.01;
delta = 1;
k_corr=0.5
tauy=tau;
tauz=tau;

Dy=1;
Dz=delta;
gamma_val=1/2;

Lx=2*pi;
Lz=75;
Ly=Lz;
Nx=64;
Ny=32;
Nz=32;

dx=Lx/Nx;
dy=Ly/Ny;
dz=Lz/Nz;

% CFL condition on dt - leapfrog is stable if cfl is satisfied.
Vmax=2*(Lz/2)/gamma_val;
dx_min=min([dx,dy,dz]);
dt=0.1*dx_min/Vmax;
display(dt)
im=sqrt(-1);

t_vec=0:dt:Tf;
n_timesteps=length(t_vec);

ck=0.95;

%**************************************************************************
% Initialize code

x=-(Lx/2):dx:((Lx/2)-dx);
y=-(Ly/2):dy:((Ly/2)-dy);
z=-(Lz/2):dz:((Lz/2)-dz);

p0=zeros(Nx,Ny,Nz);

w_y=2*Dy/tauy;
w_z=2*Dz*(1-k_corr*k_corr)/tauz;

for i=1:1
    for j=1:1
        for k=1:1
            p0(i,j,k)=(1/Lx)*(1-0.5*cos(2*x(i)))*exp(-(y(j)^2)/w_y)*exp(-(z(k)^2)/w_z);
        end
    end
end
                
nrm=sum(sum(sum(p0)))*dx*dy*dz;
p0=p0/nrm;

kx=mod((1:Nx)-ceil(Nx/2+1),Nx)-floor(Nx/2);
ky=mod((1:Ny)-ceil(Ny/2+1),Ny)-floor(Ny/2);
kz=mod((1:Nz)-ceil(Nz/2+1),Nz)-floor(Nz/2);

Kx=zeros(Nx,Ny,Nz);
Ky=zeros(Nx,Ny,Nz);
Kz=zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            Kx(i,j,k)=kx(i);
            Ky(i,j,k)=ky(j);
            Kz(i,j,k)=kz(k);
        end
    end
end

Y=zeros(Nx,Ny,Nz);
Z=zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            Y(i,j,k)=y(j);
            Z(i,j,k)=z(k);
        end
    end
end

Kx=(2*pi/Lx)*Kx;
Ky=(2*pi/Ly)*Ky;
Kz=(2*pi/Lz)*Kz;

Kx=im*Kx;
Ky=im*Ky;
Kz=im*Kz;

small=0;
ksquare_lap=small*(Kx.^2)+(Dy/(tauy^2))*(Ky.^2)+(Dz*(1-k_corr*k_corr)/(tauz^2))*(Kz.^2);

Vx=zeros(Nx,Ny,Nz);

for i=1:Nx
    for j=1:Ny
        for k=1:Nz
        x_val=x(i);
        y_val=y(j);
        z_val=z(k); 
        prefac=1;         
        Vx(i,j,k)=prefac*((w/gamma_val)+(1/gamma_val)*(-cos(x_val)+k_corr*sqrt(delta))*y_val+(z_val/gamma_val));
        end
    end
end

%**************************************************************************
% % Movie code
n_subdiv = 50;

%**************************************************************************

residual_vec=0*(1:n_timesteps);
av_vec=residual_vec;

p_hat_old=fftn(p0);
Vp_hat=fftn(Vx.*p0);
py_hat=fftn(p0.*Y);
pz_hat=fftn(p0.*Z);
dy_py_hat= Ky.*py_hat;
dz_pz_hat= Kz.*pz_hat;

src0_hat=-Kx.*Vp_hat;
src1_hat=(1/tauy)*dy_py_hat+(1/tauz)*dz_pz_hat;
src_hat=src0_hat+src1_hat;

p_hat=((1+(1-ck)*dt*ksquare_lap)./(1-ck*dt*ksquare_lap)).*p_hat_old+dt*(src_hat./(1-ck*dt*ksquare_lap));
p=real(ifftn(p_hat));

residual_val=max(max(max(abs(p-p0))));
residual_vec(1)=residual_val;

p_x=0*x;
for i=1:length(x)
    p_x(i)=sum(sum(p(i,:,:)));
end
x1=[x,pi];
p1=[p_x,p_x(1)];
av_vec(1)=sum(p1.*x1)/sum(p1);


for t_ctr=2:length(t_vec)
        
    Vp_hat=fftn(Vx.*p);
    py_hat=fftn(p.*Y);
    pz_hat=fftn(p.*Z);
    dy_py_hat= Ky.*py_hat;
    dz_pz_hat= Kz.*pz_hat;
    
    src0_hat=-Kx.*Vp_hat;
    src1_hat=(1/tauy)*dy_py_hat+(1/tauz)*dz_pz_hat;
    
    src_hat=src0_hat+src1_hat;

    p_hat_new=(1./(1-2*dt*ksquare_lap)).*p_hat_old+2*dt*(src_hat./(1-2*dt*ksquare_lap));
    
    p_hat_old=p_hat;
    p_old=real(ifftn(p_hat_old));
    
    p_hat=p_hat_new;    
    p=real(ifftn(p_hat));

    residual_val=max(max(max(abs(p_old-p))));
    residual_vec(t_ctr)=residual_val;
    for i=1:length(x)
        p_x(i)=sum(sum(p(i,:,:)));
    end
    x1=[x,pi];
    p1=[p_x,p_x(1)];
    av_vec(t_ctr)=sum(p1.*x1)/sum(p1);
    
    if(mod(t_ctr,n_subdiv)==0)
        
        p_x=(1:length(x))*0;
        
        for i=1:length(x)
            p_x(i)=sum(sum(p(i,:,:)))*dy*dz;
        end
        
        nrm1=sum(p_x)*dx;
        
        x_aug=[x,pi];
        p_aug=[p_x,p_x(1)];
        display(sum(x_aug.*p_aug)/sum(p_aug))
        
        p_y=(1:length(y))*0;
        
        for i=1:length(y)
            p_y(i)=sum(sum(p(:,i,:)))*dx*dz;
        end
        
        nrm2=sum(p_y)*dy;
        
                p_z=(1:length(z))*0;
        
        for i=1:length(z)
            p_z(i)=sum(sum(p(:,:,i)))*dx*dy;
        end
        
        nrm3=sum(p_z)*dz;
        
        plot(x,p_x/nrm1,'color','blue','linewidth',2)
        hold on
        plot(y,p_y/nrm2,'color','green','linewidth',2,'linestyle','--')
        plot(z,p_z/nrm3,'color','red','linewidth',2,'linestyle','-.')
        title(num2str(t_ctr*dt));
        hold off
        drawnow
    end
    
    
end

nrm=sum(sum(sum(p)))*dx*dy*dz;
p_joint=p/nrm;

%**************************************************************************
% Compute marginal distributions

p_x=0*(1:length(x));
p_y=0*(1:length(y));
p_z=0*(1:length(z));

for i=1:length(x)
    p_x(i)=sum(sum(p_joint(i,:,:)))*dy*dz;
end

for j=1:length(y)
    p_y(j)=sum(sum(p_joint(:,j,:)))*dx*dz;
end

for k=1:length(z)
    p_z(k)=sum(sum((p_joint(:,:,k))))*dy*dz;
end

p_xy=zeros(length(x),length(y));

for i=1:length(x)
    for j=1:length(y)
        p_xy(i,j)=sum(p_joint(i,j,:))*dz;
    end
end

p_yz=zeros(length(y),length(z));

for j=1:length(y)
    for k=1:length(z)
        p_yz(j,k)=sum(p_joint(:,j,k))*dx;
    end
end

nrm=sum(sum(p_yz))*dy*dz;
p_yz=p_yz/nrm;

p_yz_anal=zeros(length(y),length(z));

w_y=2*Dy/tauy;
w_z=2*Dz*(1-k_corr*k_corr)/tauz;

for j=1:length(y)
    for k=1:length(z)
        p_yz_anal(j,k)=exp(-y(j)^2/w_y)*exp(-z(k)^2/w_z);
    end
end

nrm=sum(sum(p_yz_anal))*dy*dz;
p_yz_anal=p_yz_anal/nrm;

%**************************************************************************
% Compute lambda distribution

lambda_range=-30:0.01:30;
p_lambda=0*lambda_range;

for i=1:length(lambda_range)
    sum_val=0;
    lambda_val=lambda_range(i);
    
    if(lambda_val==0)
        sum_val=0;
    else
    
        for ii=1:length(x)  
            g_val=-2*sin(x(ii));

%             if((abs(lambda_val)/(abs(g_val)+1e-8))>=(Ly/2))
%                 summand=0;
%             else
%                 ll_map=1+(1/dy)*((lambda_val/g_val)+(Ly/2))
%                 ll_map=floor(ll_map);
%                 summand=p_xy(ii,ll_map)*(dx/abs(g_val));
%             end
%             sum_val=sum_val+summand;

            if((abs(lambda_val)/(abs(g_val)+1e-8))>=(Ly/2))
                summand=0;
            else
                ll_map=1+(1/dy)*((lambda_val/g_val)+(Ly/2));
                ll_map=floor(ll_map);
                
                if(abs(ll_map)>Ny)
                    summand=0;
                else
                    summand=p_xy(ii,ll_map)*(dx/abs(g_val));
                end
            end
            sum_val=sum_val+summand;
            
        end

    end
    
    p_lambda(i)=sum_val;
end

[mmm,ix]=min(abs(lambda_range));

p_lambda(ix)=(p_lambda(ix-1)+p_lambda(ix+1))/2;

plot(lambda_range,p_lambda,'color','blue','linewidth',2)
drawnow

%**************************************************************************
        
end

