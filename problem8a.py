import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt1/absolute_error.csv", sep=',', names = ['x', 'err'], dtype=np.float64)

data.replace(',','')

plt.plot(data['x'].head(9), data['err'].head(9),'r-', label='n=10')
plt.plot(data['x'].iloc[13:107], data['err'].iloc[13:107], 'b-', label = 'n=100')
plt.plot(data['x'].iloc[115:1106], data['err'].iloc[115:1106], 'g-', label = 'n=1000')
plt.plot(data['x'].iloc[1118:11100].astype(float), data['err'].iloc[1118:11100], label = 'n=10000')
#plt.plot(data['x'].iloc[11120:101120], data['err'].iloc[11120:101120])
plt.xlabel('x')
plt.legend()
plt.ylabel('log10(absolute error)')

#plt.savefig
plt.show()

#ax[0].plot(data['x'].head(9), data['err'].head(9))
#ax[1].plot(data['x'].iloc[13:107], data['err'].iloc[13:107])
#ax[2].plot(data['x'].iloc[115:1106], data['err'].iloc[115:1106])
#ax[3].plot(data['x'].iloc[1118:11108].astype(float), data['err'].iloc[1118:11108])
#print(data['x'].iloc[1118:11108])
#fig.savefig('problem8a.png')


