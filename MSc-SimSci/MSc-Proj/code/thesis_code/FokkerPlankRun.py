import sys
import imp
import numpy as np
code_dir ='/home/ian/Desktop/MSc-Proj/code/thesis_code/'
data_dir ='/home/ian/Desktop/MSc-Proj/saved-datasets/'
img_save_dir = '/home/ian/Desktop/MSc-Proj/thesistemplate/imgs/'
chartDataDir='/home/ian/Desktop/MSc-Proj/chart-data/tau_0.01-T_5/'

T=float(sys.argv[1])
opt_params={}

opt_params['T']=T

if len(sys.argv) > 2:
    params_tmp=sys.argv[2].split(';')

if len(sys.argv) > 2:
    params_tmp=sys.argv[2].split(';')
    for p in params_tmp:
        kv = p.split("=")
        opt_params[kv[0]] = float(eval(kv[1]))            
                
opt_params_str=str(opt_params).replace("'","").replace(":","_").replace(" ","").replace("{","").replace("}","").replace(",","-")

chartDataDir=chartDataDir + str(opt_params_str) + '/'
import os
os.mkdir(chartDataDir)       

print "Running Fokker-Planck Solver for T = " + str(T) + "\n"
print "Params\n" 
print opt_params
fp_solver = imp.load_source('FokkerPlank', code_dir +'FokkerPlankSolve.py')
fp = fp_solver.FokkerPlank(T,opt_params)
fp.solve()

#save to file

fp_quantities={'x':None,'y':None,'z':None,'p_x':None,'p_y':None,'p_z':None,'p_xy':None,'p_yz':None, 'p_yz_anal': None,'p_lambda':None,'lambda_range':None}
for q in fp_quantities.keys():
    file_name = 'fp_<QUANT>.csv'.replace('<QUANT>',q)
    np.savetxt(chartDataDir + file_name , getattr(fp,q), delimiter=",",fmt='%1.10e')
#load from file
for q in fp_quantities.keys():
    file_name = 'fp_<QUANT>.csv'.replace('<QUANT>',q + '_tau_eq_0.01')
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
    #plt.show()
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
    #plt.show()
    plt.savefig(data_dir + 'fp_p_yz.png')
    
plot_p_yz_surface(fp_quantities['y'],fp_quantities['z'],fp_quantities['p_yz'],chartDataDir)

def plot_lambda(lambda_range,p_lambda, data_dir):
    import matplotlib.pyplot as plt
    fig = plt.figure(3)
    plt.plot(lambda_range,p_lambda)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$P_{\Lambda}$')
    #plt.show()
    plt.savefig(data_dir+'fp_p_lambda.png')

plot_lambda(fp_quantities['lambda_range'], fp_quantities['p_lambda'], chartDataDir)

from scipy import stats
stats.moment(fp_quantities['p_lambda'], moment=4)
np.mean(fp_quantities['p_lambda'])

def plot_p_x(x,p_x, data_dir):
    import matplotlib.pyplot as plt
    fig = plt.figure(4)
    plt.plot(x / np.pi,p_x)
    plt.xlabel(r'$x/\pi$')
    plt.ylabel(r'$P_{X}$')
    #plt.show()
    plt.savefig(data_dir+'fp_p_x.png')
 
plot_p_x(fp_quantities['x'], fp_quantities['p_x'], chartDataDir)
