import imp
import numpy as np
import matplotlib.pyplot as plt

code_dir ='/home/ian/Desktop/MSc-Proj/code/thesis_code/'
data_dir ='/home/ian/Desktop/MSc-Proj/saved-datasets/'
img_save_dir = '/home/ian/Desktop/MSc-Proj/thesistemplate/imgs/'

"""

Simulations and plot executions of direct solve of SDE model

"""
sde_solver = imp.load_source('SDEModelSolve', code_dir +'SDEModelSolve.py')
sde_solver.SDEModelSolve(10, 0, 2**32-1).solve_n_plot(img_save_dir)


"""

Vorticity Solve post processing 

"""
theta=np.loadtxt(data_dir + 'theta250000.csv',dtype='float', delimiter=',')
w=np.loadtxt(data_dir + 'omega250000.csv',dtype='float', delimiter=',')
dns_solver = imp.load_source('VorticitySolvePostProcessing', code_dir +'VorticitySolvePostProcessing.py')
(u,v, d,s, phi,Xangle,w_tot,lambda_unsigned,sig1_w,sig2_w,sig1_l,sig2_l,sig_lw,grad,ll,prod,theta_x,theta_y)=dns_solver.VorticitySolvePostProcessing.extract_model_parameters(w,theta)


"""

Plots

"""

plt.figure(1)
plt.imshow(w, cmap="afmhot")
plt.title("Vorticity Field  " + r'$\omega$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'vorticity-field-dns.png')

plt.figure(2)
plt.imshow(theta)
plt.title("Concentration Field  " + r'$\theta$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'concentration-field-dns.png')

plt.figure(3)
plt.imshow(lambda_unsigned)
plt.title(r'$\mu = sign(d)\sqrt{s^{2} + d^{2}}$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'mu-field-dns.png')

plt.figure(4)
plt.imshow(phi)
plt.title("The angle " + r'$\phi$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'phi-angle-dns.png')

plt.figure(5)
plt.imshow(Xangle, cmap="Spectral")
plt.title("The angle " + r'$X$')
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'X-angle-dns.png')

plt.figure(6)
plt.imshow(grad, cmap="Spectral")
plt.title("Concentration gradient")
plt.ylabel('y')
plt.xlabel('x')
plt.colorbar()
plt.savefig(img_save_dir+'theta-gradient-dns.png')

"""

Construct emperical PDF of the variables one the 

"""
hist_dict = dns_solver.VorticitySolvePostProcessing.histogram_avg(data_dir,20000,250000,-1)

(ya,xa) = hist_dict['LAMBDA']
(Xya,Xxa) = hist_dict['Xangle']
(Zya,Zxa) = hist_dict['mu']

plt.figure(1)
plt.title(r'$\Lambda$')
plt.plot(xa[:-1],ya/np.max(ya))
plt.savefig(img_save_dir+'lambda-hist-pdf.png')
plt.figure(2)
plt.title("The angle " + r'$X$')
Xxaa = Xxa[:-1]
plt.plot(Xxaa[600:1250]-6.2,Xya[600:1250]/np.max(Xya[600:1250]))
plt.savefig(img_save_dir+'angle-x-hist-pdf.png')
plt.figure(3)
plt.title('$\mu$')
plt.plot(Zxa[:-1],Zya)
plt.savefig(img_save_dir+'mu-hist-pdf.png')


"""

Run FokkerPlanck simulation

"""

import imp
import numpy as np
code_dir ='/home/ian/Desktop/MSc-Proj/code/thesis_code/'
data_dir ='/home/ian/Desktop/MSc-Proj/saved-datasets/'
img_save_dir = '/home/ian/Desktop/MSc-Proj/thesistemplate/imgs/'
chartDataDir='/home/ian/Desktop/MSc-Proj/chart-data/'

opt_params={'tau': 2}
opt_params_str=str(opt_params).replace("'","").replace(":","_").replace(" ","").replace("{","").replace("}","").replace(",","-")

chartDataDir=chartDataDir + str(opt_params_str) + '/'
import os
os.mkdir(chartDataDir)       

fp_solver = imp.load_source('FokkerPlank', code_dir +'FokkerPlankSolve.py')
fp = fp_solver.FokkerPlank(1,opt_params)
fp.solve()

#save to file

fp_quantities={'x':None,'y':None,'z':None,'p_x':None,'p_y':None,'p_z':None,'p_xy':None,'p_yz':None, 'p_yz_anal': None,'p_lambda':None,'lambda_range':None}
for q in fp_quantities.keys():
    file_name = 'fp_<QUANT>.csv'.replace('<QUANT>',q)
    np.savetxt(chartDataDir + file_name , getattr(fp,q), delimiter=",",fmt='%1.10e')
#load from file
for q in fp_quantities.keys():
    file_name = 'fp_<QUANT>.csv'.replace('<QUANT>',q)
    fp_quantities[q] = np.loadtxt(chartDataDir + file_name, dtype='float', delimiter=',')


def plot_p_xy_surface(x_param,y_param,p_xy_param,data_dir):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import numpy as np
    XX,YY = np.meshgrid(x_param,y_param)
    XX = XX
    YY = YY
    ZZ = p_xy_param.T
    colors = cm.hot((ZZ-ZZ.min())/(ZZ.max() - ZZ.min()))
    rcount, ccount, _ = colors.shape
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(XX, YY, ZZ, rcount=rcount, ccount=ccount,facecolors=colors, shade=True)
    surf.set_facecolor((0,0,0,0))
    ax.set_xlim(min(x_param), max(x_param))
    ax.set_ylim(min(y_param), max(y_param))
    ax.set_xlabel('x',fontsize=10)
    ax.set_ylabel('y', fontsize=10)
    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel(r'$P_{xy}$', fontsize=10)
    plt.show()
    plt.savefig(data_dir + 'fp_p_xy.png')

plot_p_xy_surface(fp_quantities['x'],fp_quantities['y'],fp_quantities['p_xy'],chartDataDir)

def plot_p_yz_surface(y_param, z_param, p_yz_param, data_dir):
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import cm
    import matplotlib.pyplot as plt
    YY,ZZ = np.meshgrid(y_param,z_param)
    YY = YY
    ZZ = ZZ
    P_YZ = p_yz_param.T
    colors = cm.coolwarm((P_YZ-P_YZ.min())/(P_YZ.max() - P_YZ.min()))
    rcount, ccount, _ = colors.shape
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(YY, ZZ, P_YZ, rcount=rcount, ccount=ccount,facecolors=colors, shade=True)
    surf.set_facecolor((0,0,0,0))
    ax.set_xlabel('y',fontsize=10)
    ax.set_ylabel('z', fontsize=10)
    ax.zaxis.set_rotate_label(False) 
    ax.set_zlabel(r'$P_{yz}$', fontsize=10)
    plt.show()
    plt.savefig(data_dir + 'fp_p_yz.png')
    
plot_p_yz_surface(fp_quantities['y'],fp_quantities['z'],fp_quantities['p_yz'],chartDataDir)

def plot_lambda(lambda_range,p_lambda, data_dir):
    fig = plt.figure(3)
    plt.plot(lambda_range,p_lambda)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$P_{\Lambda}$')
    plt.show()
    plt.savefig(data_dir+'fp_p_lambda.png')

plot_lambda(fp_quantities['lambda_range'], fp_quantities['p_lambda'], chartDataDir)

def plot_p_x(x,p_x, data_dir):
    fig = plt.figure(4)
    plt.plot(x / np.pi,p_x)
    plt.xlabel(r'$x/\pi$')
    plt.ylabel(r'$P_{X}$')
    plt.show()
    plt.savefig(data_dir+'fp_p_x.png')
 
plot_p_x(fp_quantities['x'], fp_quantities['p_x'], chartDataDir)


fig = plt.figure(5)
plt.plot(fp.lambda_range,fp.p_lambda)
plt.plot(xa[:-1]/xa[:-1].max(),ya/ya.max())

fig = plt.figure(6)
plt.plot(fp.x / np.pi,fp.p_x)
plt.plot(Zxa[:-1],Zya/Zya.max())
