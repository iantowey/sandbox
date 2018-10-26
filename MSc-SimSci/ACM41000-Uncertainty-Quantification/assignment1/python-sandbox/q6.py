#

def f(u, a, b):
    return a - b*u[0] - exp(-u[0])

f(0)

U = np.linspace(start = -.5, stop = .5, num = 1000, endpoint=True)
F_U = [f([u],2,-2.5) for u in U]
plt.plot(U, F_U)
plt.grid(True)
plt.text(0.65, 0, r'$f(u) = a - bu - e^{-u}$')
plt.title('f(u) vs u, (a=2,b=2.5)')
plt.xlabel('u')
plt.ylabel('f(u)')
plt.show()
