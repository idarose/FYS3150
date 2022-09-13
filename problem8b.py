import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt1/relative_error.csv", sep=',', names = ['x', 'err'], dtype=np.float64)


data.replace(',','')

fig, ax = plt.subplots(4,1, figsize= (10,10), sharex=True)

ax[0].plot(data['x'].head(9), data['err'].head(9))
ax[1].plot(data['x'].iloc[13:107], data['err'].iloc[13:107])
ax[2].plot(data['x'].iloc[115:1106], data['err'].iloc[115:1106])
ax[3].plot(data['x'].iloc[1118:11100], data['err'].iloc[1118:11100])

ax[3].set_xlabel('x')
ax[0].set_ylabel('log10 of relative error')
ax[1].set_ylabel('log10 of relative error')
ax[2].set_ylabel('log10 of relative error')
ax[3].set_ylabel('log10 of relative error')

fig.savefig('problem8b.png')


