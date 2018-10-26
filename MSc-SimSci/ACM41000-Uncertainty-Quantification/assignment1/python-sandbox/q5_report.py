from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def f(u, a, b):
    return a - b*u[0] - exp(-u[0])
    
for param in [((1.5,0.9),'1', 10),((2,3.5),'2',10),((3.1,2.5),'3',10),((3.1,0.5),'4',10),((1.5,3.5),'5',10),]:
    a,b = param[0]
    file_idx = param[1]
    u_max = param[2]
    U = np.linspace(start = 0, stop = 5, num = 10000, endpoint=True)
    F_U = [f([u],a,b) for u in U]
    plt.plot(U, F_U)
    plt.grid(True)
    plt.text(0.25, -4.8, r'$f(u) = a - bu - e^{-u}$')
    #plt.title('f(u) vs u, (a=1.1,b=2.5)')
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.show()
plt.savefig('/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignment1/q5_a.png')

f(sol.x+0.1,2,2.5)


sol = optimize.root(fun = f,x0 = [0], args = (1.1,2.5), method = 'hybr')


sol.x + 0.1

step = 0.01
a_p = np.arange(1 + step, 10 + step, step)
b_p = np.arange(-4.5 + step, 4.5 + step, step)
a_ij = []
b_ij = []
u_ij = []

for i in range(0,len(a_p)):
    for j in range(0,len(b_p)):
        sol = optimize.root(fun = f,x0 = [0], args = (a_p[i],b_p[j]), method = 'hybr')
        a_ij.append(a_p[i])
        b_ij.append(b_p[j])
        u_ij.append(sol.x[0])

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(a_ij, b_ij, u_ij, lw=step)
ax.set_xlabel("a")
ax.set_ylabel("b")
ax.set_zlabel(r'$u_{*}$')
ax.set_title("Fixed point " + r'$u_{*}$' + " dependence on constants a,b")

plt.show()
plt.savefig('/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignment1/q5_b.png')

def du_dt(u,t,a,b):
    return a - b*u - exp(-u)

from scipy import integrate
u_step = linspace(0, 5,  10001)
X, infodict = integrate.odeint(du_dt, y0 = 0, t = u_step, args = (5, 2.5), full_output=True)
X1, infodict1 = integrate.odeint(du_dt, y0 = 0, t = u_step, args = (2, 3.5), full_output=True)
X2, infodict2 = integrate.odeint(du_dt, y0 = 0, t = u_step, args = (1.5, 2.5), full_output=True)
X3, infodict3 = integrate.odeint(du_dt, y0 = 0, t = u_step, args = (1.5, 3.5), full_output=True)
infodict['message']
infodict1['message']

X2.T[0]

aa = 4



import matplotlib.pyplot as plt


u_step = linspace(0, 5,  10001)
for params in [((5,1),'b'),((2,3.5),'g'),((1.5,2.5),'r'),((1.5,3.5),'y')]:
    X, infodict = integrate.odeint(du_dt, y0 = 0, t = u_step, args = params[0], full_output=True)
    if infodict['message'] != 'Integration successful':
        plt.plot(u_step, X.T[0], params[1], label=r'$\frac{du}{d\tau} = ' + str(params[0][0]) + ' - ' + str(params[0][1]) + 'u -e^{-u}$')
    else:
        print 'error (a = {}, b = {})'.format(params[0][0], params[0][1])
        
#plt.plot(u_step, X.T[0], 'b', label=r'$\frac{du}{d\tau} = a - bu -e^{-u}$')
#plt.plot(u_step, X1.T[0], 'g', label=r'$\frac{du}{d\tau} = 2 - 3.5u -e^{-u}$')
#plt.plot(u_step, X2.T[0], 'r', label=r'$\frac{du}{d\tau} = 1.5 - 2.5u -e^{-u}$')
#plt.plot(u_step, X3.T[0], 'y', label=r'$\frac{du}{d\tau} = 1.5 - 3.5u -e^{-u}$')

plt.title(r'$\frac{du}{d\tau} = a - bu -e^{-u}$')
plt.legend(loc='best')
plt.xlabel('u')
plt.ylabel('f(u)')
plt.grid()
plt.show()


        