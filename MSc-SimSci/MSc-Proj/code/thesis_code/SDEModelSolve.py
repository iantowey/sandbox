import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

"""

Solve the SDE model directly using Euler-Maruyama algorithm.

"""
class SDEModelSolve:
    """
    Constructor.
    
    Args:
        t_periods  : mandatory integer
        rng_seed_1 : optional Y rng seed
        rng_seed_2 : optional Z rng seed
    
    Returns:
    """
    def __init__(self,t_periods, rng_seed_1=None, rng_seed_2=None):
        assert(t_periods != None or t_periods > 0)
        self.rng_seed_1 = rng_seed_1
        self.rng_seed_2 = rng_seed_2
        self.w = 0.5         #from page 11
        self.k = 0.5         #from page 11
        self.tau = 2         #from page 11
        self.Dy = 2          #from page 11
        self.Dz = 2          #from page 11
        self.dt=1e-6         #from page 11
        self.gamma=0.5       #from page 9
        self.T=self.tau*t_periods         #"integrate over T_periods time constants"
        self.strength=2      #"Weiner process of strength 2"
        self.ts = np.arange(0, self.T, self.dt)
        self.dX = np.zeros(self.ts.size)
        self.dY = np.zeros(self.ts.size)
        self.dZ = np.zeros(self.ts.size)
        
    """
    X SDE function definition.
    
    Args:
        x: value of x at time t ; x(t)
        y: value of y at time t ; y(t) 
        z: value of z at time t ; z(t)  
    
    Returns:
        float: value of x at time t+1 ; x(t+1)
    """
    def dX_sde(self,x,y,z):
        return ((self.w +(-np.cos(x) + self.k*np.sqrt(self.Dz/self.Dy))*y + z)/self.gamma)*self.dt
        
    """
    Y SDE function definition.
    
    Args:
        y:     value of y at time t ; y(t) 
        Wy:    wiener process value at time t; Wy(t)  
    
    Returns:
        float: value of y at time t+1 ; y(t+1)
    """
    def dY_sde(self,y,Wy):
        return (-y/self.tau)*self.dt + (np.sqrt(self.Dy)/self.tau)*Wy
        
    """
    Z SDE function definition.
    
    Args:
        z:     value of z at time t ; z(t) 
        Wz:    wiener process value at time t; Wz(t)  
    
    Returns:
        float: value of z at time t+1 ; z(t+1)
    """
    def dZ_sde(self,z,Wz):
        return (-z/self.tau)*self.dt + (np.sqrt(self.Dz*(1-self.k**2))/self.tau)*Wz

    """
    Solve SDE model and plot lambda and angle X.
    
    Args:
    
    Returns:
    """
    def solve_n_plot(self, save_plot_dir=None):
        np.random.seed(self.rng_seed_1)
        dWy=np.random.normal(loc = 0, scale = np.sqrt(self.dt),size=self.ts.size)
        np.random.seed(self.rng_seed_2)
        dWz=np.random.normal(loc = 0, scale = np.sqrt(self.dt),size=self.ts.size)

        for i in range(1, self.ts.size):
    
            x_t_minus_1 = self.dX[i-1]
            y_t_minus_1 = self.dY[i-1]
            z_t_minus_1 = self.dZ[i-1]
        
            y_t = y_t_minus_1 + self.dY_sde(y_t_minus_1,dWy[i])
            z_t = z_t_minus_1 + self.dZ_sde(z_t_minus_1,dWz[i])
            x_t = x_t_minus_1 + self.dX_sde(x_t_minus_1,y_t_minus_1,z_t_minus_1)
        
            self.dX[i] = x_t  
            self.dY[i] = y_t
            self.dZ[i] = z_t  
            
        dX_scaled_mod_2pi = np.array([ x % ((-1 if x < 0 else 1)*2*np.pi) for x in self.dX])/(2*np.pi)
        #plot angle X
        plt.figure(1)
        plt.plot(self.ts, dX_scaled_mod_2pi)
        plt.ylabel('X')
        plt.xlabel('t')
        plt.xticks(np.arange(0, self.T+1, step=5))
        plt.show()
        if save_plot_dir != None:
            plt.savefig(save_plot_dir + '/x-angle.png')

        
        lambda_est = -2*np.array(self.dY)*np.sin(self.dX)
        #plot lamdba
        plt.figure(2)
        plt.plot(self.ts, lambda_est)
        plt.ylabel(r'$\Lambda$')
        plt.xlabel('t')
        plt.xticks(np.arange(0, self.T+1, step=5))
        plt.show()
        if save_plot_dir != None:
            plt.savefig(save_plot_dir + '/lambda.png')
