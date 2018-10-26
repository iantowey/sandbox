from scipy.integrate import simps
import numpy as np

def f1(x):
    return x**2

x = np.array(range(0,10))

y1 = f1(x)

I1 = simps(y1,x)

print(I1)

