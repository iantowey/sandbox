import numpy as np

w=np.loadtxt('/home/ian/Desktop/MSc Proj/saved-datasets/omega100000.csv',dtype='float', delimiter=',')
t=np.loadtxt('/home/ian/Desktop/MSc Proj/saved-datasets/theta100000.csv',dtype='float', delimiter=',')

(d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y) = VorticitySolvePostProcessing().postprocess0(w,t)


(d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y) = VorticitySolvePostProcessing().postprocess0(w,t)

(xa,ya) = VorticitySolvePostProcessing().postprocess3('/home/ian/Desktop/MSc Proj/saved-datasets/',20000,100000,-1)

plt.plot(ya[:-1,],xa)

ll.shape


(aaa,bbb) = np.histogram(ll,np.arange(-1,1,0.01))

plt.plot(bbb[:-1],aaa)
