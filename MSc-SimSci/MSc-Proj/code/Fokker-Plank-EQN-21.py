import numpy as np

class FokkerPlank:
    def __init__(self,Tf):
        self.solved=False
        self.Tf = Tf
        self.x = 0
        self.delta=1
        self.gamma_val = 0.5
        self.k_corr = 0.5
        self.x = 0
        self.tau=0.01
        self.tauy=self.tau
        self.tauz=self.tau
        self.w = 0.5
        self.k = 0.5
        self.Dy = 1
        self.Dz = 1
        self.Lx=2*np.pi
        self.Lz=75
        self.Ly=self.Lz
        self.Nx=64
        self.Ny=32
        self.Nz=32
        self.dx=self.Lx/float(self.Nx)
        self.dy=self.Ly/float(self.Ny)
        self.dz=self.Lz/float(self.Nz)
        self.Vmax=2*(self.Lz/2)/self.gamma_val
        self.dx_min=min([self.dx,self.dy,self.dz])
        self.dt=0.1*self.dx_min/self.Vmax
        self.im=1j
        self.t_vec=np.arange(0,self.Tf,self.dt)
        self.n_timesteps=len(self.t_vec)
        self.ck=0.95
        self.x=np.arange(-self.Lx/float(2),(self.Lx/float(2)),self.dx)
        self.y=np.arange(-self.Ly/float(2),(self.Ly/float(2)),self.dy)
        self.z=np.arange(-self.Lz/float(2),(self.Lz/float(2)),self.dz)
        self.w_y=2*self.Dy/self.tauy
        self.w_z=2*self.Dz*(1-self.k_corr*self.k_corr)/self.tauz
        self.p0=np.zeros((self.Nx,self.Ny,self.Nz))
        self.lambda_range=np.arange(-30,30,0.01)
        for i in range(0,self.Nx):
            for j in range(0,self.Ny):
                for k in range(0,self.Nz):
                    self.p0[i,j,k] = (1/float(self.Lx)) * (1-0.5*np.cos(2*self.x[i])) * np.exp(-(self.y[j]**2)/self.w_y) * np.exp(-(self.z[k]**2)/self.w_z)
        self.nrm=sum(self.p0)*self.dx*self.dy*self.dz
        self.p0=self.p0/self.nrm

        self.kx=(([i+1 for i in range(0,self.Nx)] -np.ceil(self.Nx/2+1)) % self.Nx)-np.floor(self.Nx/2)
        self.ky=(([i+1 for i in range(0,self.Ny)] -np.ceil(self.Ny/2+1)) % self.Ny)-np.floor(self.Ny/2)
        self.kz=(([i+1 for i in range(0,self.Nz)] -np.ceil(self.Nz/2+1)) % self.Nz)-np.floor(self.Nz/2)

        self.Kx=np.zeros(shape=(self.Nx,self.Ny,self.Nz))
        self.Ky=np.zeros(shape=(self.Nx,self.Ny,self.Nz))
        self.Kz=np.zeros(shape=(self.Nx,self.Ny,self.Nz))

        for i in range(0,self.Nx):
            for j in range(0,self.Ny):
                for k in range(0,self.Nz):
                    self.Kx[i,j,k]=self.kx[i]
                    self.Ky[i,j,k]=self.ky[j]
                    self.Kz[i,j,k]=self.kz[k]

        self.Y=np.zeros(shape=(self.Nx,self.Ny,self.Nz))
        self.Z=np.zeros(shape=(self.Nx,self.Ny,self.Nz))

        for i in range(0,self.Nx):
            for j in range(0,self.Ny):
                for k in range(0,self.Nz):
                    self.Y[i,j,k]=self.y[j]
                    self.Z[i,j,k]=self.z[k]

        self.Kx=self.im*((2*np.pi/float(self.Lx))*self.Kx)
        self.Ky=self.im*((2*np.pi/float(self.Ly))*self.Ky)
        self.Kz=self.im*((2*np.pi/float(self.Lz))*self.Kz)
        
        small=0
        self.ksquare_lap=small*(self.Kx**2)+(self.Dy/(self.tauy**2))*(self.Ky**2)+(self.Dz*(1-self.k_corr*self.k_corr)/(self.tauz**2))*(self.Kz**2)

        self.Vx=np.zeros(shape=(self.Nx,self.Ny,self.Nz))

        for i in range(0,self.Nx):
            for j in range(0,self.Ny):
                for k in range(0,self.Nz):
                    self.x_val=self.x[i]
                    self.y_val=self.y[j]
                    self.z_val=self.z[k] 
                    prefac=1     
                    self.Vx[i,j,k]=prefac*((self.w/self.gamma_val)+(1/self.gamma_val)*(-np.cos(self.x_val)+self.k_corr*np.sqrt(self.delta))*self.y_val+(self.z_val/self.gamma_val))

        self.n_subdiv = 50

        self.residual_vec=np.zeros(self.n_timesteps)
        self.av_vec=self.residual_vec
        
        self.p_hat_old=np.fft.fftn(self.p0)
        self.Vp_hat=np.fft.fftn(self.Vx*self.p0)
        self.py_hat=np.fft.fftn(self.p0*self.Y)
        self.pz_hat=np.fft.fftn(self.p0*self.Z)
        self.dy_py_hat= self.Ky*self.py_hat
        self.dz_pz_hat= self.Kz*self.pz_hat
        
        self.src0_hat=-self.Kx*self.Vp_hat
        self.src1_hat=(1/self.tauy)*self.dy_py_hat+(1/self.tauz)*self.dz_pz_hat
        self.src_hat=self.src0_hat+self.src1_hat
        
        self.p_hat=((1+(1-self.ck)*self.dt*self.ksquare_lap)/(1-self.ck*self.dt*self.ksquare_lap))*self.p_hat_old+self.dt*(self.src_hat/(1-self.ck*self.dt*self.ksquare_lap))
        self.p=np.real(np.fft.ifftn(self.p_hat))
        
        self.residual_val=(abs(self.p-self.p0)).max()
        self.residual_vec[0]=self.residual_val
        
        self.p_x=0*self.x
        for i in range(0,len(self.x)):
            self.p_x[i]=sum(self.p[i,:,:])
        
        self.x1=np.concatenate([self.x, [np.pi]])
        self.p1=np.concatenate([self.p_x,[self.p_x[0]]])
        self.av_vec[0]=sum(self.p1*self.x1)/sum(self.p1)

    def __compute_marginal_distributions(self):
        self.p_x=np.zeros(len(self.x))
        self.p_y=np.zeros(len(self.y))
        self.p_z=np.zeros(len(self.z))
    
        for i in range(0,len(self.x)):
            self.p_x[i]=sum(self.p_joint[i,:,:])*self.dy*self.dz
                
        for j in range(0,len(self.y)):
            self.p_y[j]=sum(self.p_joint[:,j,:])*self.dx*self.dz
    
        for k in range(0,len(self.z)):
            self.p_z[k]=sum(self.p_joint[:,:,k])*self.dy*self.dz
    
        self.p_xy=np.zeros((len(self.x),len(self.y)))
        for i in range(0,len(self.x)):
            for j in range(0,len(self.y)):
                self.p_xy[i,j]=sum(self.p_joint[i,j,:])*self.dz
    
        self.p_yz=np.zeros((len(self.y),len(self.z)))
        for j in range(0,len(self.y)):
            for k in range(0,len(self.z)):
                self.p_yz[j,k]=sum(self.p_joint[:,j,k])*self.dx
            
        nrm=sum(self.p_yz)*self.dy*self.dz
        self.p_yz=self.p_yz/nrm
        
        self.p_yz_anal=np.zeros((len(self.y),len(self.z)))
        
        w_y=2*self.Dy/self.tauy
        w_z=2*self.Dz*(1-self.k_corr*self.k_corr)/self.tauz
        
        for j in range(0,len(self.y)):
            for k in range(0,len(self.z)):
                self.p_yz_anal[j,k]=np.exp(-self.y[j]**2/w_y)*np.exp(-self.z[k]**2/self.w_z)
                    
        nrm=sum(self.p_yz_anal)*self.dy*self.dz
        self.p_yz_anal=self.p_yz_anal/nrm
        
    def __do_solve_iter(self, idx):
        Vp_hat=np.fft.fftn(self.Vx*self.p)
        py_hat=np.fft.fftn(self.p*self.Y)
        pz_hat=np.fft.fftn(self.p*self.Z)
        dy_py_hat= self.Ky*py_hat
        dz_pz_hat= self.Kz*pz_hat
        
        self.src0_hat=-self.Kx*Vp_hat
        self.src1_hat=(1/self.tauy)*dy_py_hat+(1/self.tauz)*dz_pz_hat
        
        self.src_hat=self.src0_hat+self.src1_hat
    
        self.p_hat_new=(1/(1-2*self.dt*self.ksquare_lap))*self.p_hat_old+2*self.dt*(self.src_hat/(1-2*self.dt*self.ksquare_lap))
        
        self.p_hat_old=self.p_hat
        self.p_old=np.real(np.fft.ifftn(self.p_hat_old))
        
        self.p_hat=self.p_hat_new    
        self.p=np.real(np.fft.ifftn(self.p_hat))
    
        self.residual_val=(abs(self.p_old-self.p)).max()
        self.residual_vec[idx]=self.residual_val
        for i in range(1,len(self.x)):
            self.p_x[i]=sum(self.p[i,:,:])
    
        self.x1=np.concatenate([self.x, [np.pi]])
        self.p1=np.concatenate([self.p_x,[self.p_x[0]]])
        self.av_vec[idx]=sum(self.p1*self.x1)/sum(self.p1)
        return None

    def solve(self):
        if not self.solved:
            self.solved = True
            for t_ctr in range(1,len(self.t_vec)):
                print t_ctr
                self.__do_solve_iter(t_ctr)
            #    self.Vp_hat=np.fft.fftn(self.Vx*self.p)
            #    self.py_hat=np.fft.fftn(self.p*self.Y)
            #    self.pz_hat=np.fft.fftn(self.p*self.Z)
            #    self.dy_py_hat= self.Ky*self.py_hat
            #    self.dz_pz_hat= self.Kz*self.pz_hat
            #    
            #    self.src0_hat=-self.Kx*self.Vp_hat
            #    self.src1_hat=(1/self.tauy)*self.dy_py_hat+(1/self.tauz)*self.dz_pz_hat
            #    
            #    self.src_hat=self.src0_hat+self.src1_hat
            #
            #    self.p_hat_new=(1/(1-2*self.dt*self.ksquare_lap))*self.p_hat_old+2*self.dt*(self.src_hat/(1-2*self.dt*self.ksquare_lap))
            #    
            #    self.p_hat_old=self.p_hat
            #    self.p_old=np.real(np.fft.ifftn(self.p_hat_old))
            #    
            #    self.p_hat=self.p_hat_new    
            #    self.p=np.real(np.fft.ifftn(self.p_hat))
            #
            #    self.residual_val=(abs(self.p_old-self.p)).max()
            #    self.residual_vec[t_ctr]=self.residual_val
            #    for i in range(1,len(self.x)):
            #        self.p_x[i]=sum(self.p[i,:,:])
            #
            #    self.x1=np.concatenate([self.x, [np.pi]])
            #    self.p1=np.concatenate([self.p_x,[self.p_x[0]]])
            #    self.av_vec[t_ctr]=sum(self.p1*self.x1)/sum(self.p1)
            nrm=sum(self.p)*self.dx*self.dy*self.dz
            self.p_joint=self.p/nrm
            self.__compute_marginal_distributions()
        else:
            print 'Fokker-Plank already solved for given parameters, use methods '

    def p(self):
        return self.p_joint
        
    def p_x(self):
        return self.p_x
        
    def p_y(self):
        return self.p_y
    
    def p_z(self):
        return self.p_z

    def p_xy(self):
        return self.p_xy

    def p_yz(self):
        return self.p_yz

    def p_lambda(self):
        
        self.p_lambda=0*self.lambda_range
        for i in range(0,len(self.lambda_range)):
            sum_val=0
            lambda_val=self.lambda_range[i]
            
            if lambda_val==0:
                sum_val=0
            else:
                for ii in range(0,len(fp.x)):
                    g_val=-2*np.sin(fp.x[ii])
            
                    if((abs(lambda_val)/(abs(g_val)+1e-8))>=(fp.Ly/2)):
                        summand=0
                    else:
                        ll_map=1+(1/fp.dy)*((lambda_val/g_val)+(fp.Ly/2))
                        ll_map=int(np.floor(ll_map)-1)
                        
                        if(abs(ll_map)>fp.Ny):
                            summand=0
                        else:
                            summand=fp.p_xy[ii,ll_map]*(fp.dx/abs(g_val))
                    sum_val=sum_val+summand
            self.p_lambda[i]=sum_val
        ix=np.argmin(abs(self.lambda_range))
        self.p_lambda[ix]=(self.p_lambda[ix-1] + self.p_lambda[ix+1])/2
        return self.p_lambda
    
        
##############################################

fp = FokkerPlank(.1)

fp.solve()

#plt.plot(fp.x,fp.p_x) #PDF of X
#plt.plot(fp.y,fp.p_y) #PDF of X
#plt.plot(fp.z,fp.p_z) #PDF of X
#plt.imshow(fp.p_xy)
#plt.imshow(fp.p_yz)
#plt.imshow(fp.p_yz_anal)
plt.plot(fp.lambda_range,fp.p_lambda)
 
#chartDataDir='/home/ian/Desktop/MSc Proj/chart-data/'
#np.savetxt(chartDataDir + "fp_p_x.csv" , fp.p_x, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "fp_p_y.csv" , fp.p_y, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "fp_p_z.csv" , fp.p_z, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "fp_p_xy.csv" , fp.p_xy, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "fp_p_yz.csv" , fp.p_yz, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "p_yz_anal.csv" , fp.p_yz_anal, delimiter=",",fmt='%1.10e')
#np.savetxt(chartDataDir + "p_lambda.csv" , fp.p_lambda, delimiter=",",fmt='%1.10e')


#plt.plot(lambda_range,p_lambda)
