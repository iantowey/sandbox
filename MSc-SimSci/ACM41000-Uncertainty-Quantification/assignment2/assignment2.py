import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class HeatEquation1DFFTSolve:
    
    def __init__(self, L, D, N, u0, t_min, t_max):
        self.L = L
        self.D = D
        self.N = N
        self.u0 = u0
        self.t_min = t_min
        self.t_max = t_max
        self.x_vec = np.arange(0, L, float(L)/N)
        self.delta_t = float(t_max) / N
        self.k_vec = (2*np.pi/L) * np.arange(-N/2, N/2)
        self.un = [self.u0]
        self.solution_epsilon = 0.000001
        
    def solve(self):
        uhat = np.fft.fftshift(np.fft.fft(self.u0))
        for i in range(1,self.N+1):
            uhat_new = uhat/(1+self.D*(i*self.delta_t)*self.k_vec*self.k_vec)
            uhat = uhat_new
            self.un.append(np.fft.ifft(np.fft.ifftshift(uhat)))
        
    def plot_solution(self):
        plt.axes().set_xlim(-0.05, self.L +.05)
        plt.axes().set_ylim(-0.05, self.L +.05)
        plt.axes().grid()
        for i in range(0, len(self.un)):
            plt.plot(self.x_vec,self.un[i])
            if (max(abs(self.un[i]))-min(abs(self.un[i]))) < self.solution_epsilon:
                break

    def animate_solution(self):
        fig, ax = plt.subplots()
        line, = ax.plot(self.x_vec, self.u0, lw=2)
        ax.grid()
        xdata, ydata = [], []
        time_label = ax.text(0.05, 0.90, 'time = 0.0', transform=ax.transAxes) # initialize the time label for the graph

        def solution_data():
            for i in range(0, len(self.un)):
                yield i*self.delta_t, self.x_vec, self.un[i]

        def init():
            ax.set_ylim(-0.05, self.L+.05)
            ax.set_xlim(-0.05, self.L+.05)
            del xdata[:]
            del ydata[:]
            line.set_data(xdata, ydata)
            return line,

        def run(data):
            t, x, u = data
            time_label.set_text('time = %.3f' % t) # Display the current time to the accuracy of your liking.
            line.set_data(x, u)
            return line,time_label

        ani = animation.FuncAnimation(
            fig, 
            run, 
            solution_data, 
            blit=False, 
            interval=1000, 
            repeat=True, 
            init_func=init)
        return ani

########################
########################
##Init class
########################
########################

NN=2**14
he = HeatEquation1DFFTSolve(
    L = 1, 
    D = 1, 
    N = NN, 
    u0 = np.array([1 if float(i)/NN >= .25 and float(i)/NN <= .75 else 0 for i in range(0,NN) ]),
    t_min = 0, 
    t_max = 10)
    
########################
#Solve
########################
he.solve()

########################
#Solution Plot
########################
he.plot_solution()

########################
#Solution Animation
########################
he_ani = he.animate_solution()
