from array import array
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#N=9
data_a = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/data_prob_6a.csv", sep=',', names = ['v1', 'v2', 'v3'])
#N=99
data_b = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/data_prob_6b.csv", sep=',', names = ['v1', 'v2', 'v3'])

#vector = data_a['v1'].values
vector = data_b['v1'].values

y1 = np.zeros(len(vector[0].split()))
y2 = np.zeros(len(vector[2].split()))
y3 = np.zeros(len(vector[4].split()))
for i in range(len(vector[0].split())):
    y1[i] = float(vector[0].split()[i])
    y2[i] = float(vector[2].split()[i])
    y3[i] = float(vector[4].split()[i])

x = np.linspace(0,1,len(vector[0].split()))

plt.plot(x,y1, 'r-', label = 'Smallest eig')
plt.plot(x,y2, 'b-', label = '2. smallest eig')
plt.plot(x,y3, 'g-', label = '3. smallest eig')

plt.legend(bbox_to_anchor =(0.5,1.16), loc='upper center',fancybox=True, shadow=True, ncol=3)
plt.xlabel('x/L')
plt.grid()
plt.ylabel('v')
plt.title('v/x for n=100')
plt.savefig('C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/p6b.pdf')