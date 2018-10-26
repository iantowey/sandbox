path = '/home/ian/Desktop/ACM41000-Uncertainty-Quantification/assignment1/'

from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def f(u, a, b):
    return a - b*u[0] - exp(-u[0])
    
U = np.linspace(start = 0, stop = 5, num = 10000, endpoint=True)
for param in [((1.5,0.9),'1', 10),((2,3.5),'2',10),((3.1,2.5),'3',10),((3.1,0.5),'4',10),((1.5,3.5),'5',10),]:
    a,b = param[0]
    file_idx = param[1]
    u_max = param[2]
    F_U = [f([u],a,b) for u in U]
    plt.plot(U, F_U)
    plt.grid(True)
    plt.text(0.25, -4.8, r'$f(u) = a - bu - e^{-u}$')
    plt.xlabel('u')
    plt.ylabel('f(u)')
    plt.show()
plt.savefig(path + 'fp_diff_ab.png')
