import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt1/max_rel_err.csv", sep=',', names = ['n', 'err'], dtype=np.float64)

data.replace(',','')
plt.plot(data['n'], data['err'])

plt.show()

