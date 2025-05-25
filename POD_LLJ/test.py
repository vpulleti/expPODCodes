import os
import sys
import numpy as np
import scipy.sparse as scysparse
from pdb import set_trace as keyboard
from time import sleep
import scipy.sparse as scysparse
import scipy.sparse.linalg as spysparselinalg  # sparse linear algebra
import scipy.linalg as scylinalg               # non-sparse linear algebra
import matplotlib as mpl
mpl.use('Agg')
import pylab as plt
import pickle
from matplotlib import rc as matplotlibrc
from matplotlib import cm
import matplotlib.mlab as mlab
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path
import POD_functions
import code_functions

# please note that this program is written to read the inlet case. The name can be changed accordingly for the other cases and also the timesteps have to be changed to 2000 

# This program reads the data and stores in a 3D numpy array for U,V and 2D array for X,Y as those do not change at each time step (these can be a 2D arrays). 

grid_y = 314
grid_x = 207
num_tstp = 1000    
prefix = 'Inlet'

bin_name ="/scratch/brown/vpulleti/LLJ_PIV/binary/inlet/Inlet_" # path name can be changed to directory where you want to store the binary files

if Path(bin_name+"U").exists() and Path(bin_name+"V").exists() and Path(bin_name+"X").exists() and Path(bin_name+"Y").exists():
  
  print("The binary files exists and loading")

  with open(bin_name+"U",'rb') as file:
    U = pickle.load(file)
   
  with open(bin_name+"V",'rb') as file:
    V = pickle.load(file)
    
  with open(bin_name+"X",'rb') as file:
    X = pickle.load(file)
    
  with open(bin_name+"Y",'rb') as file:
    Y = pickle.load(file)
 
 # print("X\n",X[:,:])
 # print("U\n",U[:,:,999])
  print("The files are loaded and moving onto calculation")
else:

  print("The binary files do not exist and creating them")
  X = np.zeros((grid_y,grid_x))
  Y = np.zeros((grid_y,grid_x))
  U = np.zeros((grid_y,grid_x,num_tstp))
  V = np.zeros((grid_y,grid_x,num_tstp))
 
  for i in range(num_tstp):
    name = "Inlet000"+str("%03d"%i)+".T000.D000.P000.H000.L.vec"
    fold = "/scratch/brown/vpulleti/LLJ_PIV/inlet/" #pathname where the data files are there
    path = fold+name
    with open(path,'r') as file:
      lines = file.readlines()
      lines = lines[1:]
      x = np.zeros(len(lines))
      y = np.zeros(len(lines))
      u = np.zeros(len(lines))
      v = np.zeros(len(lines))
      for count in range(len(lines)):
        string = lines[count].split(',')
        y[count] = np.float(string[0])
        x[count] = np.float(string[1])
        u[count] = np.float(string[3])
        v[count] = np.float(string[2]) 
    X[:,:] = np.transpose(np.reshape(x,(grid_x,grid_y)))   
    Y[:,:] = np.transpose(np.reshape(y,(grid_x,grid_y)))
    U[:,:,i] = np.transpose(np.reshape(u,(grid_x,grid_y)))
    V[:,:,i] = np.transpose(np.reshape(v,(grid_x,grid_y)))
#    keyboard()
    print("Done timestep %d"%i)
  #print(U[:,:,999])
  with open(bin_name+"X",'wb') as file:
    pickle.dump(X,file)

  with open(bin_name+"Y",'wb') as file:
    pickle.dump(Y,file)

  with open(bin_name+"U",'wb') as file:
    pickle.dump(U,file)

  with open(bin_name+"V",'wb') as file:
    pickle.dump(V,file)

  print("The binary files are created and moving further")

X_new = -1.*X[:,39:165]
Y_new = Y[:,39:165]
U_new = -1.*U[:,39:165,:]
V_new = V[:,39:165,:]

U_sum = U_new.sum(2)
V_sum = V_new.sum(2)

U_mean = U_sum/float(num_tstp)
V_mean = V_sum/float(num_tstp)

U_fluc = np.zeros(np.shape(U_new))
V_fluc = np.zeros(np.shape(V_new))

for i in range(num_tstp):
  U_fluc[:,:,i] = U_new[:,:,i]-U_mean
  V_fluc[:,:,i] = V_new[:,:,i]-U_mean

#print("U",U_new)
#print("V",V_new)
#print("X",X_new)
#print("Y",Y_new)



if Path(bin_name+"eigvals").exists() and Path(bin_name+"eigvecs").exists():
  with open(bin_name+'eigvals','rb') as file:
    eigvals = pickle.load(file)
  with open(bin_name+"eigvecs",'rb') as file:
    eigvecs = pickle.load(file)

  print("The eigen valuess and eigen vectors exist and loaded")
else:
  U_sing = np.zeros((num_tstp,(U_fluc[:,:,0]).size))
  V_sing = np.zeros((num_tstp,(V_fluc[:,:,0]).size))

#keyboard()
  for i in range(num_tstp):
    U_sing[i,:] = (U_fluc[:,:,i]).flatten(order='C')
    V_sing[i,:] = (V_fluc[:,:,i]).flatten(order='C')

  A = U_sing+V_sing

  area_xy = (np.max(np.max(X_new))-np.min(np.min(X_new)))*(np.max(np.max(Y_new))-np.min(np.min(Y_new)))

  A = A*area_xy*1E-06
  coeff_mat = np.matmul(A,A.transpose())

  eigvals,eigvecs = scylinalg.eig(coeff_mat)

  with open(bin_name+"eigvals",'wb') as file:
    pickle.dump(eigvals,file)

  with open(bin_name+"eigvecs",'wb') as file:
    pickle.dump(eigvecs,file)
  print("The eigen values and eigen vectors do not exist and created")

#keyboard()
sum_eigvals = abs(eigvals).sum()

eigvals = abs(eigvals)/sum_eigvals

#print(100*eigvals)

phi_u = POD_functions.eigfunc(U_fluc,eigvecs)
phi_v = POD_functions.eigfunc(V_fluc,eigvecs)

out = POD_functions.normlze(X_new,Y_new,phi_u,phi_v)
#keyboard()
phi_u = out['phiu']
phi_v = out['phiv']

#keyboard()
b_u = POD_functions.timcoeff(X_new,Y_new,U_fluc,phi_u)
b_v = POD_functions.timcoeff(X_new,Y_new,V_fluc,phi_v)

mod_reqd = 1
u_recons_1 = POD_functions.reconst(b_u,phi_u,mod_reqd)
v_recons_1 = POD_functions.reconst(b_v,phi_v,mod_reqd)
print("Reconstructed using first mode")

mod_reqd = 10
u_recons_10 = POD_functions.reconst(b_u,phi_u,mod_reqd)
v_recons_10 = POD_functions.reconst(b_v,phi_v,mod_reqd)
print("Reconstructed using 10 modes")

mod_reqd = 100
u_recons_100 = POD_functions.reconst(b_u,phi_u,mod_reqd)
v_recons_100 = POD_functions.reconst(b_v,phi_v,mod_reqd)
print("Reconstructed using 100 modes")


out = code_functions.normalize(prefix)

vel = out['vel']
D  = out['D']
length = out['length']


U_mean = U_mean/vel
V_mean = V_mean/vel
X_new = 0.001*X_new/D
Y_new = 0.001*(Y_new-200)/D
u_recons_1 = u_recons_1/vel
v_recons_1 = u_recons_1/vel
u_recons_10 = u_recons_10/vel
v_recons_10 = v_recons_10/vel
u_recons_100 = u_recons_100/vel
v_recons_100 = v_recons_100/vel


out = code_functions.contrwritedata(X_new,Y_new,U_mean+U_fluc[:,:,1],bin_name+"u_inst_m0.dat")
out = code_functions.contrwritedata(X_new,Y_new,V_mean+V_fluc[:,:,1],bin_name+"v_inst_m0.dat")

out = code_functions.contrwritedata(X_new,Y_new,U_mean+u_recons_1[:,:,1],bin_name+"u_inst_m1.dat")
out = code_functions.contrwritedata(X_new,Y_new,V_mean+v_recons_1[:,:,1],bin_name+"v_inst_m1.dat")


out = code_functions.contrwritedata(X_new,Y_new,U_mean+u_recons_10[:,:,1],bin_name+"u_inst_m10.dat")
out = code_functions.contrwritedata(X_new,Y_new,V_mean+v_recons_10[:,:,1],bin_name+"v_inst_m10.dat")


out = code_functions.contrwritedata(X_new,Y_new,U_mean+u_recons_100[:,:,1],bin_name+"u_inst_m100.dat")
out = code_functions.contrwritedata(X_new,Y_new,V_mean+v_recons_100[:,:,1],bin_name+"v_inst_m100.dat")



figure_name = bin_name+"mode_energy.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.bar(np.arange(len(eigvals)),eigvals,align='center',alpha=0.5)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('mode number',fontsize=textFontSize)
ax.set_ylabel('residual',fontsize=textFontSize,rotation=90)
ax.set_title('Modal energy',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()
#keyboard()
cmin = np.min(U_fluc[:,:,1])
cmax = np.max(U_fluc[:,:,1])

figure_name = bin_name+"u_fluctuation.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0, figsize=(figwidth,figheight))
#fig = plt.figure()
ax_lefttop = fig.add_subplot(2,2,1)
ax_righttop = fig.add_subplot(2,2,2)
ax_leftbot = fig.add_subplot(2,2,3)
ax_rightbot = fig.add_subplot(2,2,4)
vmin = np.min(U_fluc[:,:,1])
vmax = np.max(U_fluc[:,:,1])
  

ax = ax_lefttop
plt.axes(ax)
plt.contourf(X_new,Y_new,U_mean+U_fluc[:,:,1],cmap='jet')
plt.colorbar()
#plt.clim(vmin=cmin,vmax=cmax)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('u_fluc_inst',fontsize = textFontSize)


ax = ax_righttop
plt.axes(ax)
plt.contourf(X_new,Y_new,U_mean+u_recons_1[:,:,1],cmap='jet')
plt.colorbar()
#plt.clim(vmin=cmin,vmax=cmax)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('u_recons_1',fontsize = textFontSize)

ax = ax_leftbot
plt.axes(ax)
plt.contourf(X_new,Y_new,U_mean+u_recons_10[:,:,1],cmap='jet')
plt.colorbar()
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('u_recons_10',fontsize = textFontSize)

ax = ax_rightbot
plt.axes(ax)
plt.contourf(X_new,Y_new,U_mean+u_recons_100[:,:,1],cmap='jet')
plt.colorbar()
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('u_recons_1000',fontsize = textFontSize)  

print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


figure_name = bin_name+"v_contour.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0, figsize=(figwidth,figheight))
#fig = plt.figure()
ax_lefttop = fig.add_subplot(2,2,1)
ax_righttop = fig.add_subplot(2,2,2)
ax_leftbot = fig.add_subplot(2,2,3)
ax_rightbot = fig.add_subplot(2,2,4)
cmin = np.min(V_fluc[:,:,1])
cmax = np.max(V_fluc[:,:,1])

ax = ax_lefttop
plt.axes(ax)
plt.contourf(X_new,Y_new,V_mean+V_fluc[:,:,1],cmap='jet')
#plt.clim(cmin,cmax)
plt.colorbar()
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('v_fluc_inst',fontsize = textFontSize)


ax = ax_righttop
plt.axes(ax)
plt.contourf(X_new,Y_new,V_mean+v_recons_1[:,:,1],cmap='jet')
#plt.clim(cmin,cmax)
plt.colorbar()
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('v_recons_1',fontsize = textFontSize)

ax = ax_leftbot
plt.axes(ax)
plt.contourf(X_new,Y_new,V_mean+v_recons_10[:,:,1],cmap='jet')
plt.colorbar()
#plt.clim(cmin,cmax)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('v_recons_10',fontsize = textFontSize)

ax = ax_rightbot
plt.axes(ax)
plt.contourf(X_new,Y_new,V_mean+v_recons_100[:,:,1],cmap='jet')
plt.colorbar()
#plt.clim(vmin,vmax)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel('x',fontsize=textFontSize)
ax.set_ylabel('y',fontsize=textFontSize,rotation=90)
ax.set_title('v_recons_1000',fontsize = textFontSize)

print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

