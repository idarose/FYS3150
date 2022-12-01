# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 22:14:34 2022

@author: danie
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from seaborn import heatmap
from matplotlib.animation import FuncAnimation

M=198
N = 322

data_real = pd.read_csv("real_new_U.csv", header=None)
data_imag = pd.read_csv("imag_new_U.csv", header = None)
T = []
for t in range(N):
    T.append(f'{t}')
data_real.columns = T
data_imag.columns = T

    

mat_N = np.zeros((M,M))

fig2=plt.figure()

def func(l):
    t = T[l]
    column_real = data_real[t]
    column_imag = data_imag[t]
    for i in range(M):
        for j in range(M):
            k = (M)*j +i
            a = np.sqrt(column_real[k]**2 + column_imag[k]**2)
            mat_N[i,j] = a
    return heatmap(np.transpose(mat_N), cbar=False)

animation = FuncAnimation(fig2, func, repeat=False, frames = np.arange(0,N,1), interval = 0.001)
#plt.show()
animation.save('aniwave.gif', writer="pillow", fps=40)
#%%

fig, ax = plt.subplots(1,1, figsize=(10,10))
P_data = pd.read_csv("prob.csv")

P_data.columns = ['T','P']

t = P_data['T']
y = P_data['P']

ax.grid()
ax.set_ylim(0.9999,1.0001)
ax.plot(t,y)


    
    
