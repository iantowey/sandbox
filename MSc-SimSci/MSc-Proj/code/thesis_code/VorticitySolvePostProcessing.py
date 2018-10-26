import numpy as np

class VorticitySolvePostProcessing:

    @staticmethod
    def extract_model_parameters(w,theta):
    
        (nNx,nNy)=w.shape
        Nx=nNx
        Ny=nNy
        
        L=2*np.pi
        im=1j
        
        kx=[v % Nx for v in (range(1,Nx+1)-np.ceil(Nx/2+1) ) ] - np.floor(Nx/2)
        ky=[v % Ny for v in (range(1,Ny+1)-np.ceil(Ny/2+1) ) ] - np.floor(Ny/2)
        
        Kx=np.zeros((Nx,Ny))
        Ky=np.zeros((Nx,Ny))
        
        for i in range(0,Nx):
            Kx[:,i]=kx
            Ky[i,:]=ky
        
        Kx=(2*np.pi/L)*im*Kx
        Ky=(2*np.pi/L)*im*Ky
        
        ksquare_viscous=Kx**2+Ky**2        # Laplacian in Fourier space
        ksquare_poisson=ksquare_viscous    
        ksquare_poisson[1,1]=1               # fixed Laplacian in Fourier space for Poisson's equation
        
        w_hat=np.fft.fft2(w)
        theta_hat=np.fft.fft2(theta)
        
        psi_hat = 0*w_hat
        for i in range(0,Nx):
            for j in range(0,Ny):        
                if abs(ksquare_poisson[i,j]) > 0:
                    psi_hat[i,j] = -w_hat[i,j]/complex(ksquare_poisson[i,j]) # make sure not to divide by zero here.
                    
                    
        s_hat=Kx*Ky*psi_hat
        d_hat=Ky*Ky*psi_hat-Kx*Kx*psi_hat
        d_hat=d_hat/2
        
        thetax_hat=Kx*theta_hat
        thetay_hat=Ky*theta_hat
        
        d=np.real(np.fft.ifft2(d_hat))
        s=np.real(np.fft.ifft2(s_hat))
        
        u = np.real(np.fft.ifft2( Ky*psi_hat))
        v = np.real(np.fft.ifft2(-Kx*psi_hat))
        
        theta_x=np.real(np.fft.ifft2(thetax_hat))
        theta_y=np.real(np.fft.ifft2(thetay_hat))
        
        grad=np.sqrt(theta_x**2+theta_y**2)
        
        phi=np.zeros((Nx,Ny))
        
        for i in range(0,Nx):
            for j in range(0,Ny):
                d_val=d[i,j]
       	        s_val=s[i,j]
        
           	if d_val==0:
            	   phi[i,j]=0
           	else:
           	    alpha_val=s_val/d_val
           	    nrm=2*np.sqrt(alpha_val*alpha_val+1)*(-alpha_val+np.sqrt(alpha_val*alpha_val+1))
                    if nrm < 1:
                        nrm = 1
                    nrm = np.sqrt(nrm)
                    phi[i,j] = np.arccos(1/nrm)
      		
        beta_vec=np.zeros((Ny,Nx))
        
        for i in range(0,Nx):
            for j in range(0,Ny):
        
                nrm_theta=np.sqrt(theta_x[i,j]**2+theta_y[i,j]**2)
        
                if nrm_theta==0:
                    beta_vec[i,j]=0
                else:
                    cosbeta=theta_x[i,j]/nrm_theta
                    sinbeta=theta_y[i,j]/nrm_theta
                        
                    if cosbeta >= 0:
                        if sinbeta>=0:
                            beta_val=np.arccos(cosbeta)
                        else:
                            beta_val=2*np.pi-np.arccos(cosbeta)
                    elif sinbeta >= 0:
                        beta_val=np.pi-np.arccos(abs(cosbeta))
                    else:
                        beta_val=np.pi+np.arccos(abs(cosbeta))
                    
                    beta_vec[i,j]=beta_val
        
        psi_angle=(np.pi/4)-phi
        
        psi_angle_hat=np.fft.fft2(psi_angle)
        psi_angle_x=np.real(np.fft.ifft2(Kx*psi_angle_hat))
        psi_angle_y=np.real(np.fft.ifft2(Ky*psi_angle_hat))
        conv=u*psi_angle_x+v*psi_angle_y
        w_tot=(w/2)+conv
        
        mu=np.sign(d)*np.sqrt(d**2+s**2)
        
        sig1_w=np.sum(w_tot)/(Nx*Ny)
        sig2_w=np.sum(w_tot**2)/(Nx*Ny)
        
        sig1_l=np.sum(mu)/(Nx*Ny)
        sig2_l=np.sum(mu**2)/(Nx*Ny)
        
        sig_lw=np.sum((w_tot-sig1_w)*(mu-sig1_l))/(Nx*Ny)
        
        Xangle=2*(psi_angle+beta_vec)
        
        grad_av=np.sqrt(np.sum(grad**2)/(Nx*Ny))
        
        #initialise LAMBDA array
        LAMBDA=np.zeros(np.sum(grad > 3*grad_av))
        k=0
        for i in range(0,Nx):
            for j in range(0,Ny):
                if grad[i,j] > 3*grad_av:
                    LAMBDA[k] =-2*mu[i,j]*np.sin(Xangle[i,j])
                    k = k + 1
        
        LAMBDA=np.array(LAMBDA)
        prod=s*(theta_y**2-theta_x**2)-2*d*theta_x*theta_y
        prod=prod/(grad**2)
        
        return (u,v,d,s, phi,Xangle,w_tot,mu,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,LAMBDA,prod,theta_x,theta_y)


    @staticmethod
    def histogram_avg(rootdir,nlow,nhigh,nres):
        
        hist_bin_width=0.01
        num_str=str(nlow)
        filename_theta ='theta' + num_str + '.csv'
        filename_omega ='omega' + num_str + '.csv'
    
        t=np.loadtxt(rootdir + filename_theta,dtype='float', delimiter=',')
        w=np.loadtxt(rootdir + filename_omega,dtype='float', delimiter=',')
    
        (u,v, d,s, phi,Xangle,w_tot,mu,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,LAMBDA,prod,theta_x,theta_y)=VorticitySolvePostProcessing.extract_model_parameters(w,t)
    
        (yw, xw) = np.histogram(LAMBDA,np.arange(-1,1,hist_bin_width))
        (Xyw, Xxw) = np.histogram(np.reshape(Xangle,Xangle.shape[0]*Xangle.shape[1],1),np.arange(-3,15,hist_bin_width))
        (Zyw, Zxw) = np.histogram(np.reshape(mu,mu.shape[0]*mu.shape[1],1),np.arange(-1,1,hist_bin_width))

        sum_vec_x=xw
        sum_vec_y=yw

        Xsum_vec_x=Xxw
        Xsum_vec_y=Xyw

        Zsum_vec_x=Zxw
        Zsum_vec_y=Zyw
    
        for k in np.arange(nlow,nhigh+1000,1000):
            print k
            filectr=str(k)
            num_str=str(filectr)
            filename_theta ='theta' + num_str + '.csv'
            filename_omega ='omega' + num_str + '.csv'
        
            t=np.loadtxt(rootdir + filename_theta,dtype='float', delimiter=',')
            w=np.loadtxt(rootdir + filename_omega,dtype='float', delimiter=',')
        
            (u, v, d,s, phi,Xangle,w_tot,mu,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,LAMBDA,prod,theta_x,theta_y)=VorticitySolvePostProcessing.extract_model_parameters(w,t)
        
            (yw, xw) = np.histogram(LAMBDA,np.arange(-1,1,hist_bin_width))
            (Xyw, Xxw) = np.histogram(np.reshape(Xangle,Xangle.shape[0]*Xangle.shape[1],1),np.arange(-3,15,hist_bin_width))
            (Zyw, Zxw) = np.histogram(np.reshape(mu,mu.shape[0]*mu.shape[1],1),np.arange(-1,1,hist_bin_width))

            sum_vec_x=sum_vec_x+xw
            sum_vec_y=sum_vec_y+yw

            Xsum_vec_x=Xsum_vec_x+Xxw
            Xsum_vec_y=Xsum_vec_y+Xyw    

            Zsum_vec_x=Zsum_vec_x+Zxw
            Zsum_vec_y=Zsum_vec_y+Zyw    
            
        xa=sum_vec_x/float((nhigh-nlow)/1000)
        ya=sum_vec_y/float((nhigh-nlow)/1000)
 
        Xxa=Xsum_vec_x/float((nhigh-nlow)/1000)
        Xya=Xsum_vec_y/float((nhigh-nlow)/1000)

        Zxa=Zsum_vec_x/float((nhigh-nlow)/1000)
        Zya=Zsum_vec_y/float((nhigh-nlow)/1000)
        
        return {'LAMBDA':(ya,xa),'Xangle':(Xya,Xxa),'mu':(Zya,Zxa)}


    @staticmethod
    def G(x, w):
        return np.exp(-x**2/(2*w*w))
   

    @staticmethod
    def gaussian_smoothing_loc(x,y,w):
    
        # Gaussian filter for a function on a periodic domain.
    
        x=x-2*np.pi
        L=x(len(x))-x[0]
    
        n=len(x)
        x_a=0*(range(0,3*n))
        y_a=x_a
    
        dx=abs(x[1]-x[0])
    
        for i in range(0,n):
            x_a[i]=x[i]-L
            y_a[i]=y[i]
    
        for i in range(0,n):
            x_a[i+n]=x[i]
            y_a[i+n]=y[i]
        
        for i in range(0,n):
            x_a[i+2*n]=x[i]+L
            y_a[i+2*n]=y[i]
        
        y_s=y_a
        
        for i in range(0,3*n):
            y_s[j]=sum(G(x_a-x_a[j],w)*y_a)/sum(G(x_a-x_a[j], w)) 
        
        x_ss=0*(range(0,n))
        y_ss=x_ss
    
        for i in range(0,n):
            x_ss[i] =x_a[i+n]
            y_ss[i] =y_s[i+n]
        
        for i in range(0,np.floor(n/2)):
            x_ss[i+np.floor(n/2)] = -x_ss[i]
            y_ss[i+np.floor(n/2)] = y_ss[i]
        
        (x_ss1,ix)=np.sort(x_ss)
        
        for i in range(0,len(x_ss1)):
            y_ss1[i]=y_ss[ix[i]]
        
        x_ss=x_ss1
        y_ss=y_ss1
        
        nrm=sum(y_ss)*dx
        y_ss=y_ss/nrm
        
        return (x_ss,y_ss,x_a,y_a)