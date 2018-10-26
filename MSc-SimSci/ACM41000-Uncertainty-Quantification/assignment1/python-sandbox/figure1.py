path = '/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignment1/'

from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def f(u, a, b):
    return a - b*u[0] - exp(-u[0])

step = 0.01
a_p = np.arange(1 + step, 10 + step, step)
b_p = np.arange(0 + step, 9 + step, step)
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
plt.savefig(path + '3D_plot.png')

