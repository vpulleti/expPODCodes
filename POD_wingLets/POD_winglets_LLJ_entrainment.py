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
import pylab as plt
import pickle
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path
import matplotlib.pyplot as plt
import code_functions
import POD_functions
from pdb import set_trace

mpl.rc('text', usetex = True)
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


nx = 159
ny = 154
num_tstp = 1700
loc_u = np.array([197,10])
loc_v = np.array([197,10])

# base case
'''
bin_path = "/scratch/bell/vpulleti/winglet_data/binary/BL_w0_2D/"
ascii_path = "/scratch/bell/vpulleti/winglet_data/BL_w0_2D/"
plot_path = "/scratch/bell/vpulleti/winglet_data/plots/BL_w0_2D/"
prefix = "BL_w0_2D"
bin_name ="/scratch/bell/vpulleti/winglet_data/binary/BL_w0_2D/BL_w0_2D_" # path name can be changed to directory where you want to store the binary files
'''


bin_path = sys.argv[1]
ascii_path = sys.argv[2]
plot_path = sys.argv[3]
prefix = sys.argv[4]
bin_name = sys.argv[5]

if not os.path.isdir(bin_path):
  print(bin_path+' does not exists and creating it\n')
  os.makedirs(bin_path)

if not os.path.isdir(plot_path):
  print(plot_path+' does not exists and creating it\n')
  os.makedirs(plot_path)

# reading the data
out = code_functions.bin_creation(ny,nx,num_tstp,prefix,ascii_path,bin_path,bin_name)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out
# reading the data
inlet_Ascii = "/scratch/bell/vpulleti/winglet_data/LLJ_inlet/"
inlet_path = "/scratch/bell/vpulleti/winglet_data/binary/LLJ_inlet/"
inlet_bin = "/scratch/bell/vpulleti/winglet_data/binary/LLJ_inlet/LLJ_inlet_"
out = code_functions.bin_creation(ny,nx,1600,'try',inlet_Ascii,inlet_path,inlet_bin)

u_inlet = out['u']
del out

U_inlet = u_inlet[20:130,5:148,:]
Umean_inlet = np.mean(np.mean(U_inlet,axis=2),axis=0)


# discarding the corrupted data

X_new = x[20:130,5:148]
Y_new = y[20:130,5:148]

Y_new = Y_new+372

lmts = np.array([500,520])
avgChord_len = 0.1687 #Leo's 2015 paper on wingtip vortices

out = code_functions.normalize(prefix)
D = out['D']
delta = out['delta']
U_hub = out['vel']

if prefix=="LLJ_POS_w0_2D" or prefix=="LLJ_POS_w1_2D" or prefix=="LLJ_POS_w2_2D" or prefix=="LLJ_POS_w3_2D":
  X_norm = X_new/D

elif prefix=="LLJ_POS_w0_7D" or prefix=="LLJ_POS_w1_7D" or prefix=="LLJ_POS_w2_7D" or prefix=="LLJ_POS_w3_7D":
  X_norm = X_new/D+5

if prefix=="LLJ_HUB_w0_2D" or prefix=="LLJ_HUB_w1_2D" or prefix=="LLJ_HUB_w2_2D" or prefix=="LLJ_HUB_w3_2D":
  X_norm = X_new/D

elif prefix=="LLJ_HUB_w0_7D" or prefix=="LLJ_HUB_w1_7D" or prefix=="LLJ_HUB_w2_7D" or prefix=="LLJ_HUB_w3_7D":
  X_norm = X_new/D+5


Y_norm = (Y_new-delta)/D

U_new = u[20:130,5:148,:]
V_new = v[20:130,5:148,:]

Nx = np.shape(U_new)[0]
Ny = np.shape(U_new)[1]
#keyboard()
# finding the time mean of the data
U_mean = np.mean(U_new,axis=2)
V_mean = np.mean(V_new,axis=2)


# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U_new))
V_fluc = np.zeros(np.shape(V_new))

for i in range(num_tstp):
  U_fluc[:,:,i] = U_new[:,:,i]-U_mean
  V_fluc[:,:,i] = V_new[:,:,i]-V_mean


uv_real = np.zeros(np.shape(U_new))
vv_real = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_real[:,:,i] = -U_fluc[:,:,i]*V_fluc[:,:,i]
  vv_real[:,:,i] = V_fluc[:,:,i]**2

uv_realMean = np.mean(uv_real,axis=2)/U_hub**2
vv_realMean = np.mean(vv_real,axis=2)/U_hub**2

Uuv_realMean = U_mean*uv_realMean/U_hub


UMean_2D = U_mean
Vmean_2D = V_mean


vort_real3D = code_functions.get_vorticity(X_new,Y_new,U_new,V_new)

if Path(bin_name+"eigvals").exists() and Path(bin_name+"eigvecs").exists() and Path(bin_path+"phi_u_normlzd").exists() and Path(bin_path+"time_coeff_u").exists() and Path(bin_path+"phi_v_normlzd").exists() and Path(bin_path+"time_coeff_v").exists():
  with open(bin_name+'eigvals','rb') as file:
    eigvals = pickle.load(file)
  with open(bin_name+"eigvecs",'rb') as file:
    eigvecs = pickle.load(file)
  with open(bin_name+'phi_u_normlzd','rb') as file:
    phi_u = pickle.load(file)
  with open(bin_path+"time_coeff_u",'rb') as file:
    timcoeff_u = pickle.load(file)
  with open(bin_name+'phi_v_normlzd','rb') as file:
    phi_v = pickle.load(file)
  with open(bin_path+"time_coeff_v",'rb') as file:
    timcoeff_v = pickle.load(file)

  print("The eigen values, eigen vectors, phi_vel and time_coeff exist and loaded")
else:
  U_sing = np.zeros((num_tstp,(U_fluc[:,:,0]).size))
  V_sing = np.zeros((num_tstp,(V_fluc[:,:,0]).size))

#keyboard()
  for i in range(num_tstp):
    U_sing[i,:] = (U_fluc[:,:,i]).flatten()
    V_sing[i,:] = (V_fluc[:,:,i]).flatten()

  A = np.concatenate((U_sing,V_sing),axis=1)

  coeff_mat = (1./(num_tstp-1))*np.matmul(A,A.transpose())

  eigvals,eigvecs = scylinalg.eig(coeff_mat)

  phi_vel = POD_functions.eigfunc(A,eigvecs)

  phi_vel = POD_functions.normlze(phi_vel)

  phi_u = phi_vel[:int(Ny*Nx),:]
  phi_v = phi_vel[int(Ny*Nx):,:]

  with open(bin_path+"eigvals",'wb') as file:
    pickle.dump(eigvals,file)

  with open(bin_path+"eigvecs",'wb') as file:
    pickle.dump(eigvecs,file)

  with open(bin_path+"phi_u_normlzd",'wb') as file:
    pickle.dump(phi_u,file)

  with open(bin_path+"phi_v_normlzd",'wb') as file:
    pickle.dump(phi_v,file)  

  timcoeff_u = POD_functions.timcoeff(U_sing,phi_u)
  timcoeff_v = POD_functions.timcoeff(V_sing,phi_v)

  with open(bin_path+"time_coeff_u",'wb') as file:
    pickle.dump(timcoeff_u,file)

  with open(bin_path+"time_coeff_v",'wb') as file:
    pickle.dump(timcoeff_v,file)

  print("The eigen values, eigen vectors, phi_vel_normlzd and time_coeff did not exist and created")

sum_eigvals = abs(eigvals).sum()

eigvals = abs(eigvals)/sum_eigvals

figwidth  = 12
figheight = 10
lineWidth = 3
textFontSize = 32
gcafontSize = 28

fig = plt.figure(1, figsize=(figwidth,figheight))
ax = plt.axes()
plt.loglog(np.arange(np.size(eigvals)-1)+1,eigvals[:-1],'r',linewidth=2)
ax.set_xlabel(r'n',fontsize=textFontSize)
ax.set_ylabel(r'$\lambda_n/\beta$',fontsize=textFontSize)
#plt.xlim([0,800])
plt.ylim([0.0001,0.1])
plt.setp(ax.get_xticklabels(), fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.savefig(plot_path+prefix+'_energySpecEvals.pdf',dpi=500)
print("saved file "+prefix+'_energySpecEvals.pdf')
plt.close()

cum_eigvals = np.cumsum(eigvals)

Q_12 = np.zeros((Nx,Ny,num_tstp))

for mode in range(1,num_tstp+1):
  u_recons = POD_functions.reconst(timcoeff_u,phi_u,mode-1,mode,Nx,Ny)
  v_recons = POD_functions.reconst(timcoeff_v,phi_v,mode-1,mode,Nx,Ny)
  
  uv_recons = np.zeros(np.shape(U_new))
  for i in range(num_tstp):
    uv_recons[:,:,i] = -u_recons[:,:,i]*v_recons[:,:,i]

  Q_12[:,:,mode-1] = np.mean(uv_recons,axis=2)/U_hub**2

u_reconsTot = POD_functions.reconst(timcoeff_u,phi_u,0,num_tstp,Nx,Ny)
v_reconsTot = POD_functions.reconst(timcoeff_v,phi_v,0,num_tstp,Nx,Ny)

uv_reconsTot = np.zeros(np.shape(U_new))

for i in range(num_tstp):
  uv_reconsTot[:,:,i] = -u_reconsTot[:,:,i]*v_reconsTot[:,:,i]

Q_12T = np.mean(uv_reconsTot,axis=2)/U_hub**2

data = {}
data['Xnorm'] = X_norm
data['Ynorm'] = Y_norm
data['uvNorm_Modal'] = Q_12
data['uvNorm_Total'] = Q_12T
data['Readme'] = ["The data contains normalized X and Y with D, and uv for each mode (uvNorm_Modal) and uv for full signal (uvNorm_Total)"]

with open(bin_path+prefix+"dataForQ12.pickle",'wb') as file:
  pickle.dump(data,file)

print("Saved data for Q12 "+prefix)

'''
figwidth  = 12
figheight = 10
lineWidth = 3
textFontSize = 32
gcafontSize = 28

fig = plt.figure(1, figsize=(figwidth,figheight))
ax = plt.axes()
plt.loglog(np.arange(np.size(eigvals)),Qn_12,'r',linewidth=2)
ax.set_xlabel(r'n',fontsize=textFontSize)
ax.set_ylabel(r'$Q^n_{12}/Q^T_{12}$',fontsize=textFontSize)
#plt.xlim([0,800])
plt.ylim([1E-06,0.02])
plt.setp(ax.get_xticklabels(), fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.savefig(plot_path+prefix+'_Qn12.pdf',dpi=500)
print("save file "+prefix+'_Qn12.pdf')
plt.close()
'''




