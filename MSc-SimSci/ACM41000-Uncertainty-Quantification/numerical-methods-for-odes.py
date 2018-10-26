### numerical methods for ODE;s

def func(x,y):
    return x + 2*y


##1. Eulers method

def euler(f,x_0, y_0, max_iters, delta_h):
    
    x_t_minus_1 = x_0
    y_t_minus_1 = y_0
    x_t = 0
    y_t = 0
    cnter = 0
    while(True):
        print('{:d} = estimate({:>10f},{:>10f}) true({:>10f},{:>10f})').format(cnter, x_t_minus_1, y_t_minus_1, x_t_minus_1, 0.25*exp(2*x_t_minus_1)-0.5*x_t_minus_1 - 0.25 )
        cnter = cnter + 1
        x_t = x_t_minus_1 + delta_h 
        y_t = y_t_minus_1 + f(x_t_minus_1, y_t_minus_1)*delta_h
        if cnter > max_iters:
            break
        else:
            x_t_minus_1 = x_t
            y_t_minus_1 = y_t
    

euler(func, 0, 0, 10, 0.25)

##2. Improved Eulers method

def improved_euler(f,x_0, y_0, max_iters, delta_h):
    x_t_minus_1 = x_0
    y_t_minus_1 = y_0
    x_t = 0
    y_t = 0
    cnter = 0
    while(True):
        print('{:d} = estimate({:>10f},{:>10f}) true({:>10f},{:>10f})').format(cnter, x_t_minus_1, y_t_minus_1, x_t_minus_1, 0.25*exp(2*x_t_minus_1)-0.5*x_t_minus_1 - 0.25 )
        cnter = cnter + 1
        y_t_trial = y_t_minus_1 + f(x_t_minus_1, y_t_minus_1)*delta_h
        y_t = y_t_minus_1 + 0.5*(f(x_t_minus_1, y_t_minus_1) + f(x_t_minus_1, y_t_trial))*delta_h
        x_t = x_t_minus_1 + delta_h 
        if cnter > max_iters:
            break
        else:
            x_t_minus_1 = x_t
            y_t_minus_1 = y_t
    

improved_euler(func, 0, 0, 10, 0.25)

##3. Runga Kutta method

def runga_kutta(f, x_0, y_0, max_iters, delta_t):
    x_t_minus_1 = x_0
    y_t_minus_1 = y_0
    x_t = 0
    y_t = 0
    cnter = 0;
    while(True):
        cnter = cnter + 1
        print('{:d} = estimate({:>10f},{:>10f}) true({:>10f},{:>10f})').format(cnter, x_t_minus_1, y_t_minus_1, x_t_minus_1, 0.25*exp(2*x_t_minus_1)-0.5*x_t_minus_1 - 0.25 )
        k1 = f(x_t_minus_1, y_t_minus_1)*delta_t
        k2 = f(x_t_minus_1, y_t_minus_1 + 0.5*k1)*delta_t
        k3 = f(x_t_minus_1, y_t_minus_1 + 0.5*k2)*delta_t
        k4 = f(x_t_minus_1, y_t_minus_1 + k3)*delta_t
        y_t = y_t_minus_1 + (k1 + 2*k2 + 2*k3 + k4) / 6        
        x_t = x_t_minus_1 + delta_t 
        if cnter > max_iters:
            break
        else:
            x_t_minus_1 = x_t
            y_t_minus_1 = y_t
    
runga_kutta(func, 0, 0, 10, 0.25)
