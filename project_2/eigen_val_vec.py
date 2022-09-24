import numpy as np 
import matplotlib.pyplot as plt 


N = 6
i = np.linspace(1,N,N)
h = 1/N 

a = -1/(h**2)
d = 2/(h**2)


def eig_val(a, d, i):
    N = len(i)
    return d + 2*a*np.cos (i*np.pi/(N+1))

print(eig_val(a,d,i))

v_mat = np.zeros((N,N), dtype=float)
for i in range(N):
    i += 1
    vi = np.zeros(N)
    for k in range(N):
        k += 1
        vi_el = np.sin(k*i*np.pi/(N+1))
        k -= 1
        vi[k] = vi_el 

    i -=1
    v_norm = vi/(np.linalg.norm(vi))
    v_mat[:,i] = vi

print(v_mat)


"""
The igenvalues are correct but not the eigenvectors:(
"""

"""
>  python3 eigen_val_vec.py
[  7.13024151  27.10873427  55.97849276  88.02150724 116.89126573
 136.86975849]
[[ 0.43388374  0.78183148  0.97492791  0.97492791  0.78183148  0.43388374]
 [ 0.78183148  0.97492791  0.43388374 -0.43388374 -0.97492791 -0.78183148]
 [ 0.97492791  0.43388374 -0.78183148 -0.78183148  0.43388374  0.97492791]
 [ 0.97492791 -0.43388374 -0.78183148  0.78183148  0.43388374 -0.97492791]
 [ 0.78183148 -0.97492791  0.43388374  0.43388374 -0.97492791  0.78183148]
 [ 0.43388374 -0.78183148  0.97492791 -0.97492791  0.78183148 -0.43388374]]
 """