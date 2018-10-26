
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class VorticitySolve:
    def __init__(self,T_final, T_init=None, omega_init_file=None, theta_init_file=None):
        self.T_final = T_final
        self.T_init = T_init if T_init != None else 0;
        self.Lx = 2*np.pi
        self.Ly = 2*np.pi
        self.nu_p = 5.9e-30
        self.nu_0 = 0.05
        self.p = 8
        self.dt = 1e-3
        self.im = 1j
        self.Nx = 256
        self.Ny = 256
        self.dx = self.Lx/float(self.Nx)
        self.dy = self.Ly/float(self.Ny)
        self.x = np.array([v*self.dx for v in range(0,self.Nx)])
        self.y = np.array([v*self.dx for v in range(0,self.Ny)])
        (self.X,self.Y) = np.meshgrid(self.x,self.y)
        self.N_timesteps = np.floor(self.T_final/self.dt)
        self.N_timesteps_init = np.floor(self.T_init/self.dt)
        self.kx_vec = [v*2*np.pi/self.Lx for v in range(-self.Nx/2,self.Nx/2)]
        self.ky_vec = [v*2*np.pi/self.Ly for v in range(-self.Ny/2,self.Ny/2)]
        self.Ksq = np.zeros((self.Nx,self.Ny))
        self.Kx = np.zeros((self.Nx,self.Ny))
        self.Ky = np.zeros((self.Nx,self.Ny))
        for i in range(0,self.Nx):
            self.kx_val = self.kx_vec[i]
            for j in range(0,self.Ny):
                self.ky_val = self.ky_vec[j]    
                self.Ksq[i,j] = (self.kx_val**2)+(self.ky_val**2)
                self.Kx[i,j] = self.kx_val
                self.Ky[i,j] = self.ky_val
        self.Ksq_laplace = self.Ksq
        self.Ksq_laplace[(self.Nx/2)+1,(self.Ny/2)+1] = 1
        self.A0 = 1
        self.R = 0.9
        self.kmax = 9*(2*np.pi/self.Lx)
        self.kmin = 7*(2*np.pi/self.Lx)
        self.Mx_scale = np.zeros((self.Nx,self.Ny))
        for i in range(0,self.Nx):
            for j in range(0,self.Ny):
                k_radius = np.sqrt(self.kx_vec[i]**2+self.ky_vec[j]**2)
                if k_radius>self.kmin  and k_radius<self.kmax :
                    self.Mx_scale[i,j] = 1
        self.A0_num = self.A0*self.Nx*self.Ny/np.sqrt(np.pi*(self.kmax*self.kmax-self.kmin*self.kmin))*(np.sqrt(2)/1)
        self.omega0 = np.zeros((self.Nx,self.Ny))
        self.sigma = 0.1*self.Lx
        self.exponenttheta = ((self.X -self.Lx/2)**2 + (self.Y- self.Ly/2)**2)/(2*self.sigma**2)
        self.theta0 = np.exp(-self.exponenttheta)
        self.dealias = np.ones((self.Nx,self.Ny))
        self.abs_kx_max = self.kx_vec[self.Nx-1]
        self.abs_ky_max = self.ky_vec[self.Ny-1]
        for i in range(0,self.Nx):
            for j in range(0,self.Ny):        
                abs_kx = abs(self.kx_vec[i])
                abs_ky = abs(self.ky_vec[j])
                if ((abs_kx>((2/float(3))*self.abs_kx_max)) and (abs_ky>((2/float(3))*self.abs_ky_max)) ):
                    self.dealias[i,j] = 0
        self.omega = self.omega0
        self.theta = self.theta0
        
        if omega_init_file != None and theta_init_file != None:
            self.theta=np.loadtxt(theta_init_file,dtype='float', delimiter=',')
            self.omega=np.loadtxt(omega_init_file,dtype='float', delimiter=',')
            
        self.forcing_hat = np.zeros((self.Nx,self.Ny))
        self.omega_arr = [self.omega]

    def solve(self):
        print '***************************************************************'
        print '**    Number of iterations to complete ' + str((self.N_timesteps+1)-self.N_timesteps_init)
        print '***************************************************************'
        for t_ctr in range(int(self.N_timesteps_init)+1,int(self.N_timesteps)+1):
            
            omega_hat=np.fft.fftshift(np.fft.fft2(self.omega))
            theta_hat=np.fft.fftshift(np.fft.fft2(self.theta))
            
            psi_hat=np.zeros((self.Nx,self.Ny),'complex')
            for i in range(0,self.Nx):
                for j in range(0,self.Ny):        
                    if abs(self.Ksq_laplace[i,j]) > 0:
                        psi_hat[i,j]=omega_hat[i,j]/complex(self.Ksq_laplace[i,j]) # make sure not to divide by zero here.
            
            dx_psi_hat=self.im*self.Kx*psi_hat
            dy_psi_hat=self.im*self.Ky*psi_hat
            dx_omega_hat=self.im*self.Kx*omega_hat
            dy_omega_hat=self.im*self.Ky*omega_hat
            dx_theta_hat=self.im*self.Kx*theta_hat  
            dy_theta_hat=self.im*self.Ky*theta_hat  

            dx_psi=np.real(np.fft.ifft2(np.fft.ifftshift(dx_psi_hat)))
            dy_psi=np.real(np.fft.ifft2(np.fft.ifftshift(dy_psi_hat)))
            dx_omega=np.real(np.fft.ifft2(np.fft.ifftshift(dx_omega_hat)))
            dy_omega=np.real(np.fft.ifft2(np.fft.ifftshift(dy_omega_hat)))
            dx_theta=np.real(np.fft.ifft2(np.fft.ifftshift(dx_theta_hat)))
            dy_theta=np.real(np.fft.ifft2(np.fft.ifftshift(dy_theta_hat)))
                    
            u=dy_psi
            v=-dx_psi
                                
            conv_theta     = u*dx_theta + v*dy_theta        
            conv_theta_hat = np.fft.fftshift(np.fft.fft2(conv_theta))
 
            u_dot_grad_omega = u*dx_omega + v*dy_omega
            u_dot_grad_omega_hat=np.fft.fftshift(np.fft.fft2(u_dot_grad_omega))
       
            Mx_random=np.random.random((self.Nx,self.Ny))
            random_part_hat=self.Mx_scale*np.exp(2*np.pi*self.im*Mx_random)
            self.forcing_hat=np.sqrt(1-self.R*self.R)*self.A0_num*random_part_hat+self.R*self.forcing_hat
                
            #vorticity            
            omega_hat_numerator=omega_hat*(1-0.5*self.dt*self.nu_p*(self.Ksq**self.p))-self.dt*u_dot_grad_omega_hat+self.dt*(self.forcing_hat-self.nu_0*omega_hat)
            omega_hat_denominator=1+0.5*self.dt*self.nu_p*(self.Ksq**self.p)
            omega_hat_new=self.dealias*(omega_hat_numerator/omega_hat_denominator)
            
            # invection / diffusion eqn
            numerator_theta = theta_hat - self.dt*conv_theta_hat;
            denominator_theta = (1 + self.dt*((-1)**self.p)*((-1)**self.p)*self.nu_p*((self.Ksq**self.p)));
            theta_hat_new = self.dealias*(numerator_theta/denominator_theta)

            self.theta=np.real(np.fft.ifft2(np.fft.ifftshift(theta_hat_new)))
            self.omega=np.real(np.fft.ifft2(np.fft.ifftshift(omega_hat_new)))
            print 'iter  ' + str(t_ctr) #+ ' ' + str(sum(self.omega)) + ' ' + str(self.theta[0,0])
            if t_ctr % 1000 == 0 :
                np.savetxt("/tmp/omega" + str(t_ctr) + ".csv", self.omega, delimiter=",",fmt='%1.64e')
                np.savetxt("/tmp/theta" + str(t_ctr) + ".csv", self.theta, delimiter=",",fmt='%1.64e')
        plt.imshow(self.omega)
        plt.imshow(self.theta)

    def save_data_frames(self):
        for i in range(0, len(self.omega_arr)):
            np.savetxt("/tmp/" + str(i) + "_omega.csv", self.omega_arr[i], delimiter=",")

    def view_animation(self,save=False):
        fig = plt.figure()

        # ims is a list of lists, each row is a list of artists to draw in the
        # current frame; here we are just animating one artist, the image, in
        # each frame
        ims = []
        for i in range(0,len(self.omega_arr)):
            im = plt.imshow(self.omega_arr[i], animated=True)
            ims.append([im])
        
        #create anaimation
        ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)

        #save animation to disk
        if save:
            ani.save('/tmp/Vorticity2D.mp4',bitrate=1000, dpi=100)

        plt.show()
        
                
########################
#init
########################

#vs = VorticitySolve(100)
    
########################
#Solve
########################
#vs.solve()

########################
#Save data to filesystem
########################
#vs.save_data_frames()

########################
#View animation
########################
#vs.view_animation(True)