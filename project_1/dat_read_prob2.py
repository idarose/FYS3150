import matplotlib.pyplot as plt
import pandas as pd
import csv
#headers = ['count','x','u']
data = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt1/data.csv", names=['x','y'], sep=',')

plt.plot(data['x'], data['y'])
plt.title('n=10')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.show()
