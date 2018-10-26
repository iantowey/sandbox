import matplotlib.pyplot as plt
from scipy import integrate

#function def
def du_dt(u,t,a,b):
    return a - b*u - exp(-u)

path = '/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignment1'
#population size
N = 10000
for param in [((1.5,0.9),'1', 3),((2,3.5),'2',3),((3.1,2.5),'3',3),((3.1,0.5),'4',3),((1.5,3.5),'5',3),((1.1,0.01),'6',50)]:
    a, b = param[0]
    file_idx = param[1]
    tau_max = param[2]
    u_step = linspace(0, tau_max,  10001)
    #solve ode
    U, infodict = integrate.odeint(du_dt, y0 = 0, t = u_step, args = (a,b), full_output=True)
    if infodict['message'] != 'Integration successful':
        #rescale x,y,z reletive to the population size
        x_x0 = exp(-U.T[0]) * N / a
        y_y0 = (a - b*U.T[0] - exp(-U.T[0])) * N / a
        z_z0 = (b*U.T[0]) * N / a
        #plot x_x0, y_y0, z_z0
        fig = plt.figure(facecolor='w')
        ax = fig.add_subplot(111, axis_bgcolor='#dddddd', axisbelow=True)
        ax.plot(u_step, x_x0, 'b', alpha=0.5, lw=2, label='Susceptible')
        ax.plot(u_step, y_y0, 'r', alpha=0.5, lw=2, label='Infected')
        ax.plot(u_step, z_z0, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
        ax.set_xlabel(r'$\tau$')
        ax.set_ylabel('Number People')
        ax.legend(loc='best')
        plt.title( 'Population N = 10000\n'+ r'$\frac{du}{d\tau} = ' + str(a) + ' - ' + str(b) + 'u -e^{-u}$')
        plt.legend(loc='best')
        plt.grid()
        plt.show()
        plt.savefig(path + '/q7_' + file_idx+ '.png')
    else:
        print 'error (a = {}, b = {})'.format(a, b)
        
