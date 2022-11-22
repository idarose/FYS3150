# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 16:16:27 2022

@author: danie
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 20})
strin = 'L_40.csv'

data_hist = pd.read_csv('energy_T1_ordered.csv')
data_hist.columns = ['E']


data = pd.read_csv(strin)
data2 = pd.read_csv('L_60.csv')
data.columns = ['T', 'C', 'X']
data2.columns = ['T', 'C', 'X']

data3 = pd.read_csv('L_80.csv')
data3.columns = ['T', 'C', 'X']
data4 = pd.read_csv('L_100.csv')
data4.columns = ['T', 'C', 'X']


y = data['X']
x = data['T']
y2 = data2['X']
x2 = data2['T']
y3 = data3['X']
x3 = data3['T']
y4 = data4['X']
x4 = data4['T']
i1 = y.idxmax()
i2 = y2.idxmax()
i3 = y3.idxmax()
i4 = y4.idxmax()

y_hist = data_hist['E']
var = y_hist.var()
b = y_hist.max()
a = y_hist.min()
delt_E = 0.01
nr = int((b-a)/delt_E)

fig2 = plt.figure(figsize = (10,10))
fig2 = y_hist.hist(density=True, bins=nr+2)
fig2.set_title('Probability density for T=2.4 $J/k_b$ starting with \n a random system')
fig2.set_xlabel('Energy per spin [$J$]')
fig2.set_ylabel('Probability density')

#fig.set_xlim(-2,-1.90) #Optional
#plt.savefig('T2_unordered_hist.pdf')



fig1, ax = plt.subplots(1,1, figsize=(10,10))
ax.set_title('Magnetic susceptibility as a function of temperature')
ax.set_xlabel('Temperature [$J/k_b$]')
ax.set_ylabel('Susceptibility [$J^{-1}$]')
ax.grid()
ax.plot(x,y)
ax.plot(x2,y2)
ax.plot(x3,y3)
ax.plot(x4,y4)
ax.legend(['L=40', 'L=60', 'L=80', 'L=100'])
ax.plot(x[i1],y[i1], 'bo')
ax.plot(x2[i2],y2[i2], 'o', color = 'orange')
ax.plot(x3[i3],y3[i3], 'o', color = 'green')
ax.plot(x4[i4],y4[i4], 'o', color = 'red')
fig1.savefig('Xi_T.pdf')

a, b = np.polyfit([1/40, 1/60, 1/80, 1/100], [x[i1], x2[i2], x3[i3], x4[i4]],1)
xi = np.linspace(0,1/30,100)
yi = a*xi + b
fig = plt.figure(figsize=(10,10))
plt.plot(0,2.269, '*', color='red', label = 'exact')
plt.plot(1/40,x[i1], 'bo')
plt.plot(xi,yi)
plt.legend(['Exact value for $T_c(\infty)$', 'Critical temperatures','Fitted line'])
plt.plot(1/60, x2[i2], 'bo')
plt.plot(1/80, x3[i3], 'bo')
plt.plot(1/100, x4[i4], 'bo')
plt.plot(0,b, 'bo')
plt.ylim(2.25,2.32)
plt.title('Approximation of critical temperature at $L=\infty$')
plt.xlabel('$L^{-1}$')
plt.ylabel('Temperature [$J/k_b$]')


plt.grid()

#fig.savefig('T_C_inf.pdf')




