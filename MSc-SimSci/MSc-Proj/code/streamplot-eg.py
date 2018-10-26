import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import numpy as np
code_dir ='/home/ian/Desktop/MSc-Proj/code/thesis_code/'

theta=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/theta250000.csv',dtype='float', delimiter=',')
w=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/omega250000.csv',dtype='float', delimiter=',')
theta=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/theta1000.csv',dtype='float', delimiter=',')
w=np.loadtxt('/home/ian/Desktop/MSc-Proj/saved-datasets/omega1000.csv',dtype='float', delimiter=',')
import imp
dns_solver = imp.load_source('VorticitySolvePostProcessing', code_dir +'VorticitySolvePostProcessing.py')
(u,v,d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y)=dns_solver.VorticitySolvePostProcessing.extract_model_parameters(w,theta)


fig = plt.figure(1)
Y, X = np.mgrid[0:256, 0:256]
U=u
V=v
#strm = plt.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
plt.streamplot(X, Y, U, V, density=[0.5, 1])
#fig.colorbar(strm.lines)

plt.set_title('Varying Line Width')
