import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
        
L = 1
D = 1
x_min = 0
x_max = 1
u_min = 0
u_max = 1
t_min = 0
t_max = 10
N = 2**12
x_vec = np.arange(0,L,float(L)/N)
                 
Delta_t = float(t_max) / N
u0 = [ 1 if float(i)/N >= .25 and float(i)/N <= .75 else 0 for i in range(0,N) ]
k_vec = (2*pi/L)*np.arange(-N/2,N/2)

un = []
un.append(u0)
uhat = np.fft.fftshift(np.fft.fft(u0))
for i in range(1,N+1):
    uhat_new = uhat/(1+D*(i*Delta_t)*k_vec*k_vec)
    uhat = uhat_new
    un.append(np.fft.ifft(np.fft.ifftshift(uhat)))

plt.axes().set_xlim(-x_min - 0.05, x_max +.05)
plt.axes().set_ylim(-u_min - 0.05, u_max +.05)
plt.plot(x_vec,un[0])
for i in range(1,len(un)):
        plt.plot(x_vec,un[i])

#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************

fig, ax = plt.subplots()
line, = ax.plot(x_vec, u0, lw=2)
ax.grid()
xdata, ydata = [], []
time_label = ax.text(0.05, 0.90, 'time = 0.0', transform=ax.transAxes) # initialize the time label for the graph

def data_gen():
    for i in range(0, len(un)):
        yield i*Delta_t, x_vec, un[i]

def init():
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlim(-0.05, 1.05)
    del xdata[:]
    del ydata[:]
    line.set_data(xdata, ydata)
    return line,

def run(data):
    # update the data
    t, x, u = data
    time_label.set_text('time = %.3f' % t) # Display the current time to the accuracy of your liking.
    line.set_data(x, u)
    return line,time_label

ani = animation.FuncAnimation(fig, run, data_gen, blit=False, interval=1000, repeat=True, init_func=init)
plt.show()

#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************
#*****************************************************************************************************************************************

def SolveHeatEquationFFT(u0, k_vec, t, D):
    uhat = np.fft.fftshift(np.fft.fft(u0))
    uhat_new = uhat/(1+D*t*k_vec*k_vec)
    uhat = uhat_new
    un.append(np.fft.ifft(np.fft.ifftshift(uhat)))
















class FFTHeat1D_test(object):
    """
    Solves the 1D heat equation, dU/dt = d**2U/dx**2, using FFT.
    Inputs:
    - initial condition, a twice differentiable periodic function
    - dt, size of the time step.
    """

    def __init__(self, u, dt):
        self.ux = u
        self.uk = self.compute_to_k(u)
        self.n = len(self.uk)
        self.dt = dt
        self.e = self.get_derivConst()

    def get_derivConst(self):
        kVec = np.arange(-self.n/2,self.n/2)
        derivConst = -(2*np.pi*kVec)**2*self.dt
        return np.exp(derivConst)

    def compute_to_k(self, u):
        " Converts vector from spatial to frequency domain. "
        u = np.fft.fft(u)
        return np.fft.fftshift(u)
    
    def compute_to_x(self, uk):
        " Converts vector from frequency to spatial domain. "
        uk = np.fft.ifftshift(uk)
        return np.fft.ifft(uk)

    def time_step(self):
        " Takes a step dt in time using FFT. "
        self.uk = self.e * self.uk
        self.ux = np.real(self.compute_to_x(self.uk))