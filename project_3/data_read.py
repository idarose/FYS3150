import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#
#Check the directory of your files. 
#

#Read data from z-t 
data = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data1_zt.csv", sep=',', names = ['z', 't'], dtype=np.float64)
data.replace(',','')
#Plot zt
plt.plot(data['t'],data['z'], 'xkcd:pumpkin')
plt.xlabel('t[$\mu$s]', fontsize= 20)
plt.ylabel('z[$\mu$m]', fontsize = 20)
plt.title('Oscillations of a particle in z-direction', fontsize = 20)
plt.savefig('C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/z_t_plot.pdf')

#Read data from files
data1_w_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_xy_wint.csv", sep=',', names = ['x1', 'y1', 'x2', 'y2'], dtype=np.float64)
data1_w_int.replace(',','')
data1_wo_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_xy_woint.csv", sep=',', names = ['x1', 'y1', 'x2', 'y2'], dtype=np.float64)
data1_wo_int.replace(',','')
data1_px = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_xvx_wint.csv", sep=',', names = ['x1', 'vx1', 'x2', 'vx2'], dtype=np.float64)
data1_px.replace(',','')
data1_px_wo_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_xvx_woint.csv", sep=',', names = ['x1', 'vx1', 'x2', 'vx2'], dtype=np.float64)
data1_px_wo_int.replace(',','')
data1_pz = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_zvz_wint.csv", sep=',', names = ['z1', 'vz1', 'z2', 'vz2'], dtype=np.float64)
data1_pz.replace(',','')
data1_pz_wo_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_zvz_woint.csv", sep=',', names = ['z1', 'vz1', 'z2', 'vz2'], dtype=np.float64)
data1_pz_wo_int.replace(',','')
data1_3d = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data1_3d_w_int.csv", sep=',', names = ['x', 'y', 'z'], dtype=np.float64)
data1_3d.replace(',','')
data1_3d_wo_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data1_3d_wo_int.csv", sep=',', names = ['x', 'y', 'z'], dtype=np.float64)
data1_3d_wo_int.replace(',','')

data2_3d = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data2_3d_w_int.csv", sep=',', names = ['x', 'y', 'z'], dtype=np.float64)
data2_3d.replace(',','')
data2_3d_wo_int = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data2_3d_wo_int.csv", sep=',', names = ['x', 'y', 'z'], dtype=np.float64)
data2_3d_wo_int.replace(',','')

data_trap0 = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_trap_f0.csv", sep=',', names = ['num0', 'num1', 'num2', 'wv'], dtype=np.float64)
data_trap0.replace(',','')

data_trap_wint = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_trap_fine_w_int.csv", sep=',', names = ['num', 'wv'], dtype=np.float64)
data_trap_wint.replace(',','')
data_trap_woint = pd.read_csv("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/data_trap_fine_wo_int.csv", sep=',', names = ['num', 'wv'], dtype=np.float64)
data_trap_woint.replace(',','')

#Plot and save phaseplots in figure
fig = plt.figure(figsize = (10,20), tight_layout = True)
ax1 = fig.add_subplot(4, 1, 1)
ax2 = fig.add_subplot(4, 1, 2,sharex = ax1, sharey = ax1)
ax3 = fig.add_subplot(4, 1, 3)
ax4 = fig.add_subplot(4, 1, 4, sharex = ax3, sharey = ax3)

ax1.plot(data1_px['x1'], data1_px['vx1'], 'xkcd:minty green',label='P1 w/ int')
ax1.plot(data1_px_wo_int['x1'], data1_px_wo_int['vx1'], 'xkcd:barbie pink',label='P1 w/o int')
ax2.plot(data1_px['x2'], data1_px['vx2'], 'xkcd:minty green',label='P2 w/ int')
ax2.plot(data1_px_wo_int['x2'], data1_px_wo_int['vx2'], 'xkcd:barbie pink',label='P2 w/o int')
ax1.legend(fontsize = 20)
ax1.legend(bbox_to_anchor =(0.6,1.06), loc='lower right',fancybox=True, shadow=True, ncol=3, fontsize = 20)
ax1.set_xlabel('$v_x$[$\mu$m/s]', fontsize = 20)
ax1.set_ylabel('x [$\mu$m]', fontsize = 20)
ax2.legend(bbox_to_anchor =(0.6,1.06), loc='lower right',fancybox=True, shadow=True, ncol=3, fontsize = 20)
ax2.set_xlabel('$v_x$ [$\mu$m/s]', fontsize = 20)
ax2.set_ylabel('x [$\mu$m]', fontsize = 20)

ax3.plot(data1_pz['z1'], data1_pz['vz1'], 'xkcd:minty green',label='P1 w/ int')
ax4.plot(data1_pz['z2'], data1_pz['vz2'], 'xkcd:minty green',label='P2 w/ int')
ax3.plot(data1_pz_wo_int['z1'], data1_pz_wo_int['vz1'], 'xkcd:barbie pink',label='P1 w/o int')
ax4.plot(data1_pz_wo_int['z2'], data1_pz_wo_int['vz2'], 'xkcd:barbie pink',label='P2 w/o int')
ax3.legend(bbox_to_anchor =(0.6,1.06), loc='lower right',fancybox=True, shadow=True, ncol=3, fontsize = 20)
ax3.set_xlabel('$v_z$ [$\mu$m/s]', fontsize = 20)
ax3.set_ylabel('z [$\mu$m]', fontsize = 20)
ax4.legend(bbox_to_anchor =(0.6,1.06), loc='lower right',fancybox=True, shadow=True, ncol=3, fontsize = 20)
ax4.set_xlabel('$v_z$ [$\mu$m/s]', fontsize = 20)
ax4.set_ylabel('z [$\mu$m]', fontsize = 20)

fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/Phaseplot.pdf")

#Plot and save trajectories
fig, axs = plt.subplots(2,1,figsize = (10,15), tight_layout=True)
axs[0].plot(data1_w_int['x1'], data1_w_int['y1'], 'xkcd:minty green',label='P1 w/ int')
axs[0].plot(data1_wo_int['x1'], data1_wo_int['y1'], 'xkcd:barbie pink',label='P1 w/o int')
axs[1].plot(data1_w_int['x2'], data1_w_int['y2'], 'xkcd:minty green',label='P2 w/ int')
axs[1].plot(data1_wo_int['x2'], data1_wo_int['y2'], 'xkcd:barbie pink',label='P2 w/o int')
axs[0].legend(bbox_to_anchor =(0.6,1.02), loc='lower right',fancybox=True, shadow=True, ncol=3, fontsize = 20)
axs[0].set_xlabel('x [$\mu$m]', fontsize = 20)
axs[0].set_ylabel('y [$\mu$m]', fontsize = 20)
axs[0].axis('equal')
axs[1].legend(bbox_to_anchor =(0.6,-0.2), loc='lower right',fancybox=True, shadow=True, ncol=2, fontsize = 20)
axs[1].set_xlabel('x [$\mu$m]', fontsize = 20)
axs[1].axis('equal')
axs[1].set_ylabel('y [$\mu$m]', fontsize = 20)
fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/Posplot.pdf")

#Plot and save 3D figures for particles 1 and 2
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection = '3d')
ax.plot3D(data1_3d['x'], data1_3d['y'], data1_3d['z'], 'xkcd:minty green', label='W/ Interaction')
ax.plot3D(data1_3d_wo_int['x'], data1_3d_wo_int['y'], data1_3d_wo_int['z'], 'xkcd:barbie pink', label='W/o Interaction')
ax.set_xlabel('x[$\mu$m]', fontsize = 20)
ax.set_ylabel('y[$\mu$m]', fontsize = 20)
ax.set_zlabel('z[$\mu$m]', fontsize = 20)
ax.legend(fontsize = 20)
ax.set_title('Particle 1 in Penning Trap', fontsize = 20)
fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/3d_1.pdf")
fig = plt.figure(figsize = (10,10))
ax = plt.axes(projection = '3d')
ax.plot3D(data2_3d['x'], data2_3d['y'], data2_3d['z'])
ax.plot3D(data2_3d_wo_int['x'], data2_3d_wo_int['y'], data2_3d_wo_int['z'])
ax.plot3D(data2_3d['x'], data2_3d['y'], data2_3d['z'], 'xkcd:minty green', label='W/ Interaction')
ax.plot3D(data2_3d_wo_int['x'], data2_3d_wo_int['y'], data2_3d_wo_int['z'], 'xkcd:barbie pink', label='W/o Interaction')
ax.set_xlabel('x[$\mu$m]', fontsize = 20)
ax.set_ylabel('y[$\mu$m]', fontsize = 20)
ax.set_zlabel('z[$\mu$m]', fontsize = 20)
ax.legend(fontsize = 20)
ax.set_title('Particle 2 in Penning Trap', fontsize = 20)
fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/3d_2.pdf")

#Plot and save frequency plots 
fig, axs = plt.subplots(3,1, figsize=(10,10), sharey=True, sharex = True)
axs[0].plot(data_trap0['wv'], data_trap0['num0'],'xkcd:pumpkin', label='f0 = 0.1')
axs[0].set_ylabel('%particles in trap',fontsize = 20)
axs[1].plot(data_trap0['wv'], data_trap0['num1'],'xkcd:pumpkin', label='f1 = 0.4')
axs[1].set_ylabel('%particles in trap', fontsize = 20)
axs[2].plot(data_trap0['wv'], data_trap0['num2'],'xkcd:pumpkin', label='f0 = 0.7')
axs[2].set_xlabel(r'$\omega_v$', fontsize = 20)
axs[2].set_ylabel('%particles in trap', fontsize = 20)
axs[0].set_title('f = 0.1', fontsize = 20); axs[1].set_title('f = 0.4', fontsize = 20); axs[2].set_title('f = 0.7', fontsize = 20)
fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/particles_in_trap_f0.pdf")

#Plot and save fine frequncy search
fig= plt.figure(figsize = (10,10))
ax = plt.axes()
ax.plot(data_trap_wint['wv'], data_trap_wint['num'], 'xkcd:red pink',label='w/o int')
ax.set_xlabel('$\omega_v$', fontsize = 20)
ax.set_ylabel('%particles in trap', fontsize = 20)
ax.plot(data_trap_woint['wv'], data_trap_woint['num'], 'xkcd:shamrock',label='w/ int')
ax.grid()
ax.legend(bbox_to_anchor = (0.5,1.11), loc='upper center',fancybox=True, shadow=True, ncol=2, fontsize = 20)
fig.savefig("C:/Users/Lenovo/Documents/Høst22/FYS4150/Prosjekt3/particles_in_trap_fine.pdf")
