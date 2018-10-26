
import numpy as np
code_dir ='/home/ian/Desktop/MSc-Proj/code/thesis_code/'

theta=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/theta250000.csv',dtype='float', delimiter=',')
w=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/omega250000.csv',dtype='float', delimiter=',')
import imp
dns_solver = imp.load_source('VorticitySolvePostProcessing', code_dir +'VorticitySolvePostProcessing.py')
(d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y)=dns_solver.VorticitySolvePostProcessing.extract_model_parameters(w,theta)

#plt.plot(np.sum(theta_x,axis=0)/256,theta_y)

sum(grad > 3)

#https://matplotlib.org/examples/color/colormaps_reference.html

plt.imshow(Xangle, cmap="Spectral")
plt.imshow(phi, cmap="Spectral")
plt.imshow(theta)
plt.imshow(grad)
plt.title("Tracer Gradient Field  " + r'$\theta$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig('/home/ian/Desktop/tracer-field.png')


plt.imshow(w, cmap="afmhot")
plt.title("Vorticity Field  " + r'$\omega$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig('/home/ian/Desktop/MSc Proj/vorticity-field.png')


np.max(phi)

rootdir='/home/ian/Desktop/MSc-Proj/saved-datasets/'
hist_dict = VorticitySolvePostProcessing.histogram_avg(rootdir,20000,250000,-1)

(ya,xa) = hist_dict['LAMBDA']
(Xya,Xxa) = hist_dict['Xangle']
(Zya,Zxa) = hist_dict['mu']
plt.plot(xa[:-1],ya)
plt.plot(Xxa[:-1],Xya)
plt.plot(Zxa[:-1],Zya)

Xya[0]

##
#
#   variance of omega v t , check NS simulations settle down to steady state.
#
#
#
import numpy as np

rootdir='/home/ian/Desktop/MSc Proj/saved-datasets/'
t=[]
omega_var=[]
for i in np.arange(1000,251000,1000):
    filename_omega ='omega' + str(i) + '.csv'
    w=np.loadtxt(rootdir+filename_omega,dtype='float', delimiter=',')
    t.append(i)
    N=w.shape[0]
    omega_var.append((1/float(N**2))*sum(w**2))
    print rootdir + filename_omega

time=np.array(t)/1000
chartDataDir='/home/ian/Desktop/MSc-Proj/chart-data/'
np.savetxt(chartDataDir + "time.csv" , time, delimiter=",",fmt='%1.10e')
np.savetxt(chartDataDir + "omega_var.csv", omega_var, delimiter=",",fmt='%1.64e')

        
    
import matplotlib.pyplot as plt
time=np.loadtxt(chartDataDir + "time.csv", dtype='float', delimiter=',')
omega_variance=np.loadtxt(chartDataDir + "omega_var.csv", dtype='float', delimiter=',')

1/(0.1*np.sqrt(np.mean(omega_variance)))

1/0.05

plt.plot(time, omega_variance)
plt.ylabel(r'$|| \omega ||^{2}_{2}$',fontsize=10)
plt.xlabel('t', fontsize=10)
plt.show()



VorticitySolve(250,150,'/home/ian/Desktop/MSc Proj/saved-datasets/omega150000.csv','/home/ian/Desktop/MSc Proj/saved-datasets/theta150000.csv').solve()


import os
dir_path = os.path.dirname(os.path.realpath(__file__))
cwd = os.getcwd()
import sys
sys.path.append('/home/ian/Desktop/MSc\ Proj/code')
import thesis_code.SDEModelSolve
from thesis_code import SDEModelSolve


