import numpy as np
import matplotlib.pyplot as plt

#Loading time data
data = np.loadtxt("time_algorithms_tocalc_std.csv", delimiter=",")

#Defining data arrays for special and general algoritm
general = data[:,0] # [s]
special = data[:,1] # [s]


#Function to calculate mean time and
#standard error of the mean
def std(t):
    """
    t: time measurements array

    returns:
    [0] mean time
    [1] standard error
    [2] standard error of the mean
    """

    N = len(t)

    #Finding mean time
    m = np.mean(t)

    #Calculating the standard error
    s = 0
    for i in range(N):
        s += (t[i] - m)**2

    std = 1/(N-1) * np.sqrt(s)

    #Calculating standard error of the mean
    std_m = std/(np.sqrt(N))

    return m, std, std_m 

#Creating list to index the data
indices = [-1,20,40,60,80,100,120]

G_ray = np.zeros(6)
S_ray = np.zeros(6)
x_ray = np.zeros(6)


#Printing the results
for i in range(6):
    i1 = indices[i] + 1
    i2 = indices[i+1]
    tg, ts = general[i1:i2], special[i1:i2]

    Gm, Gstd, Gstd_m = std(tg)
    Sm, Sstd, Sstd_m = std(ts)

    G_ray[i] = Gm 
    S_ray[i] = Sm 
    x_ray[i] = i+1 
    
    print(f"For N = N^{i+1}:")
    print(f"t_general = {Gm:.3e} +- {Gstd_m:.3e} s.")
    print(f"t_special = {Sm:.3e} +- {Sstd_m:.3e} s.")
    
plt.title("Runtime measurements")
plt.scatter(x_ray, G_ray, color="red", label="General")
plt.scatter(x_ray, S_ray, color="blue", label="Special")
plt.legend()
plt.xlabel("log_10(N_steps)")
plt.ylabel("Runtime [s]")
plt.savefig("Runtime.pdf")



















