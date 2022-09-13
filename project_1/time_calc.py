import numpy as np

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

#Printing the results
for i in range(6):
    i1 = indices[i] + 1
    i2 = indices[i+1]
    tg, ts = general[i1:i2], special[i1:i2]

    Gm, Gstd, Gstd_m = std(tg)
    Sm, Sstd, Sstd_m = std(ts)
    
    print(f"For N = N^{i+1}:")
    print(f"t_general = {Gm:.3e} +- {Gstd_m:.3e} s.")
    print(f"t_special = {Sm:.3e} +- {Sstd_m:.3e} s.")
    









