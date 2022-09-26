import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data_a = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/data_prob_5.csv", sep=',', names = ['n', 'iter'], dtype=np.float64)
data_b = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/dense_boii_data.csv", sep=',', names = ['n', 'iter'], dtype=np.float64)

data_a.replace(',',''); data_b.replace(',','')

fig, axs = plt.subplots(2,1, figsize = (10,10), tight_layout=True)

axs[0].plot((data_a['n']), data_a['iter'])
axs[0].set_xlabel('n')
axs[0].set_ylabel('iterations')
axs[0].grid()
axs[0].set_title('N/iterations for tridiagonal A')
axs[1].plot((data_b['n']), data_b['iter'])
axs[1].set_xlabel('(n)')
axs[1].set_ylabel('iterations')
axs[1].grid()
axs[1].set_title('N/iterations for dense A')
fig.savefig('C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt2/5_ab.pdf')
fig.show()