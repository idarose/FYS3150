"""
Created on Wed Sep  7 12:49:37 2022

@author: danie
"""
import numpy as np
import matplotlib.pyplot as plt

def g(n):
    h = 1/n
    b = np.ones(n-1)*2
    g = np.zeros(n-1)
    x = np.linspace(0,1,n+1)
    a = np.ones(n-1)*(-1)
    c = a
    
    for i in range(1,n-1):
        g[i] = h**2 *100*np.exp(-10*x[i])
    
    for i in range(n-2):
        b[i+1] -= c[i] * a[i]/b[i]
        g[i+1] -= g[i]*a[i]/b[i]
        
    for i in reversed(range(1,n-1)):
        g[i] = g[i]/b[i]
        g[i-1] += g[i]
    g = np.pad(g, (1,1))
    return x, g

n = [10, 100, 1000, 10000]
fig, ax = plt.subplots(len(n),1)
fig.suptitle('Numeriske beregninger av poissons ligning med source term f(x)')
for i in range(len(n)):
    ni = n[i]
    x, gi = g(ni)
    ax[i].plot(x,gi)
    ax[i].set_title(f'n={ni}:')
plt.show()
