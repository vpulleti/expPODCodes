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


ny = 314
nx = 207
num_tstp = 2000
loc_u = np.array([197,10])
loc_v = np.array([197,10])

# base case
bin_path = "/scratch/bell/vpulleti/LLJ_PIV/binary/single_7D_L/single_7D_L_"
ascii_path = "/scratch/bell/vpulleti/LLJ_PIV/single_7D_L/"
plot_path = "/scratch/bell/vpulleti/LLJ_PIV/binary/single_7D_L/"
prefix = "single_7D_L"
bin_name ="/scratch/bell/vpulleti/LLJ_PIV/binary/single_7D_L/single_7D_L_" # path name can be changed to directory where you want to store the binary files


if not os.path.exists('bin_path'):
    os.makedirs('bin_path')

if not os.path.exists('plot_path'):
    os.makedirs('plot_path')

# reading the data
out = code_functions.bin_creation(ny,nx,num_tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# reading the data
inlet_Ascii = "/scratch/bell/vpulleti/LLJ_PIV/inlet/"
inlet_bin = "/scratch/bell/vpulleti/LLJ_PIV/binary/inlet/Inlet_"
out = code_functions.bin_creation(ny,nx,1000,'Inlet',inlet_Ascii,inlet_bin)

u_inlet = out['u']
del out

U_inlet = -1.*u_inlet[32:,39:165,:]
U_inlet_sum = U_inlet.sum(2)

Umean_inlet = U_inlet_sum/float(1000)

Umean_inlet = np.mean(Umean_inlet,axis=1)

# discarding the corrupted data

X_new = -1.*x[32:,39:165]
Y_new = y[32:,39:165]

out = code_functions.normalize(prefix)
D = out['D']
delta = out['delta']
U_hub = out['vel']
vel_norm = code_functions.normalize('Inlet')['vel']


if prefix =="Single_2D_L" or prefix=='single_2D_H' or prefix=='Single_2D_M':
  X_norm = (X_new+D)/D
elif prefix=="single_7D_L" or prefix=='Single_7D_H' or prefix=='single_7D_M' or prefix=='array_7D_M' or prefix=='array_7D_L' or prefix=='array_7D_H':
  X_norm = (X_new+7*D-X_new[10,0])/D
else:
  X_norm = (X_new-X_new[10,0])/D

Y_norm = (Y_new-delta)/D


U_new = -1.*u[32:,39:165,:]
V_new = v[32:,39:165,:]


Ny = np.shape(U_new)[0]
Nx = np.shape(U_new)[1]
#keyboard()
# finding the time mean of the data
U_sum = U_new.sum(2)
V_sum = V_new.sum(2)

U_mean = U_sum/float(num_tstp)
V_mean = V_sum/float(num_tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U_new))
V_fluc = np.zeros(np.shape(V_new))

for i in range(num_tstp):
  U_fluc[:,:,i] = U_new[:,:,i]-U_mean
  V_fluc[:,:,i] = V_new[:,:,i]-V_mean


uv_real = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_real[:,:,i] = -U_fluc[:,:,i]*V_fluc[:,:,i]

uv_realMean = np.mean(uv_real,axis=2)

Uuv_realMean = U_mean*uv_realMean/U_hub**3

UMean_2D = U_mean
Vmean_2D = V_mean

U_mean = np.mean(U_mean,axis=1)
V_mean = np.mean(V_mean,axis=1)


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
  #A = U_sing+V_sing

  #area_xy = 1E-06*(np.max(np.max(X_new))-np.min(np.min(X_new)))*(np.max(np.max(Y_new))-np.min(np.min(Y_new)))

  #A = A*area_xy
  coeff_mat = (1./(num_tstp-1))*np.matmul(A,A.transpose())

  eigvals,eigvecs = scylinalg.eig(coeff_mat)

  phi_vel = POD_functions.eigfunc(A,eigvecs)

  phi_vel = POD_functions.normlze(phi_vel)

  phi_u = phi_vel[:np.int(Ny*Nx),:]
  phi_v = phi_vel[np.int(Ny*Nx):,:]

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


#keyboard()
sum_eigvals = abs(eigvals).sum()

eigvals = 100* abs(eigvals)/sum_eigvals


cum_eigvals = np.cumsum(eigvals)

figwidth  = 12
figheight = 10
lineWidth = 3
textFontSize = 32
gcafontSize = 28

fig = plt.figure(1, figsize=(figwidth,figheight))
ax1 = plt.axes()
ax2 = ax1.twinx()
ax1.bar(np.arange(520)+1,eigvals[:520],2)
ax2.plot(np.arange(520)+1,cum_eigvals[:520],'r',linewidth=2)
ax1.set_xlabel(r'\textrm{Modes}',fontsize=textFontSize)
ax1.set_ylabel(r'\textrm{Modal energy}',fontsize=textFontSize)
ax2.set_ylabel(r'\textrm{Cumulative energy} \%',fontsize=textFontSize)
plt.xlim([1,520])
plt.setp(ax1.get_xticklabels(), fontsize=gcafontSize)
plt.setp(ax1.get_yticklabels(),fontsize=gcafontSize)
plt.setp(ax2.get_yticklabels(),fontsize=gcafontSize)
plt.savefig(plot_path+prefix+'_eigvals.pdf',dpi=500)
print("save file "+prefix+'_eigvals.pdf')
plt.close()

'''
mod_reqd = 10
u_recons_10 = POD_functions.reconst(timcoeff_u,phi_u,0,mod_reqd,Ny,Nx)
v_recons_10 = POD_functions.reconst(timcoeff_v,phi_v,0,mod_reqd,Ny,Nx)

print("Reconstructed using 10 mode")

uv_10 = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_10[:,:,i] = -u_recons_10[:,:,i]*v_recons_10[:,:,i]

uv_10Mean = np.mean(uv_10,axis=2)

turb_prod_10 = uv_10Mean*np.gradient(UMean_2D,Y_new[2,10]-Y_new[1,10],axis=0)

Uuv_10Mean = UMean_2D*uv_10Mean/U_hub**3


mod_reqd = 50
u_recons_50 = POD_functions.reconst(timcoeff_u,phi_u,0,mod_reqd,Ny,Nx)
v_recons_50 = POD_functions.reconst(timcoeff_v,phi_v,0,,mod_reqd,Ny,Nx)
print("Reconstructed using 50 mode")

uv_50 = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_50[:,:,i] = -u_recons_50[:,:,i]*v_recons_50[:,:,i]

uv_50Mean = np.mean(uv_50,axis=2)

Uuv_50Mean = UMean_2D*uv_50Mean/U_hub**3


mod_reqd = 100
u_recons_100 = POD_functions.reconst(timcoeff_u,phi_u,0,mod_reqd,Ny,Nx)
v_recons_100 = POD_functions.reconst(timcoeff_v,phi_v,0,mod_reqd,Ny,Nx)

print("Reconstructed using 100 mode")

uv_100 = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_100[:,:,i] = -u_recons_100[:,:,i]*v_recons_100[:,:,i]

uv_100Mean = np.mean(uv_100,axis=2)

Uuv_100Mean = UMean_2D*uv_100Mean/U_hub**3
'''
mod_reqd = 350
u_recons_210 = POD_functions.reconst(timcoeff_u,phi_u,0,mod_reqd,Ny,Nx)
v_recons_210 = POD_functions.reconst(timcoeff_v,phi_v,0,mod_reqd,Ny,Nx)
print("Reconstructed using 210 mode")

uv_210 = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_210[:,:,i] = -u_recons_210[:,:,i]*v_recons_210[:,:,i]

uv_210Mean = np.mean(uv_210,axis=2)

Uuv_210Mean = UMean_2D*uv_210Mean/U_hub**3


#finding SSM
u_recons_ssm = POD_functions.reconst(timcoeff_u,phi_u,mod_reqd,num_tstp,Ny,Nx)
v_recons_ssm = POD_functions.reconst(timcoeff_v,phi_v,mod_reqd,num_tstp,Ny,Nx)
print("Reconstructed SSM")

uv_ssm = np.zeros(np.shape(U_new))
for i in range(num_tstp):
  uv_ssm[:,:,i] = -u_recons_ssm[:,:,i]*v_recons_ssm[:,:,i]

uv_ssmMean = np.mean(uv_ssm,axis=2)

Uuv_ssmMean = UMean_2D*uv_ssmMean/U_hub**3


data = {}

data['Uuv_realMean'] = Uuv_realMean
data['Uuv_210Mean'] = Uuv_210Mean
data['Uuvu_ssmMean'] = Uuv_ssmMean
data['X_norm'] = X_norm
data['Y_norm'] = Y_norm

with open(bin_path+prefix+"data.pickle",'wb') as file:
  pickle.dump(data,file)

figwidth  = 19
figheight = 15
lineWidth = 3
textFontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


#fig = plt.figure(2, figsize=(figwidth,figheight))
fig,axs = plt.subplots(2, 2, figsize=(figwidth,figheight))

ax_lefttop = axs[0,0]
ax_righttop = axs[0,1]
ax_leftbot = axs[1,0]
ax_rightbot = axs[1,1]

min_100 = np.min(np.min(Uuv_100Mean))
max_100 = np.max(np.max(Uuv_100Mean))
min_act = np.min(np.min(Uuv_realMean))
max_act = np.max(np.max(Uuv_realMean))
min_10 = np.min(np.min(Uuv_10Mean))
max_10 = np.max(np.max(Uuv_10Mean))
min_210 = np.min(np.min(Uuv_210Mean))
max_210 = np.max(np.max(Uuv_210Mean))


min_avg = np.mean([min_act,min_10,min_100,min_210])
max_avg = np.mean([max_act,max_10,max_100,max_210])

im1 = ax_lefttop.contourf(X_norm,Y_norm,Uuv_realMean,np.linspace(min_act,max_act,100),cmap='jet')
for c in im1.collections:
    c.set_edgecolor("face")
ax_lefttop.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_lefttop.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_lefttop.set_title(r'\textrm{Multi-scale}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im1, ax=axs[0, 0], orientation='vertical',shrink=0.8)




im2 = ax_righttop.contourf(X_norm,Y_norm,Uuv_10Mean,np.linspace(min_avg,max_avg,100),cmap='jet')
for c in im2.collections:
    c.set_edgecolor("face")
ax_righttop.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_righttop.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_righttop.set_title(r'\textrm{10 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im2, ax=axs[0, 1], orientation='vertical',shrink=0.8)



im3 = ax_leftbot.contourf(X_norm,Y_norm,Uuv_100Mean,np.linspace(min_act,max_act,100),cmap='jet')
for c in im3.collections:
    c.set_edgecolor("face")
ax_leftbot.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_leftbot.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_leftbot.set_title(r'\textrm{100 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im3, ax=axs[1, 0], orientation='vertical',shrink=0.8)



im4 = ax_rightbot.contourf(X_norm,Y_norm,Uuv_210Mean,np.linspace(min_act,max_act,100),cmap='jet')
for c in im4.collections:
    c.set_edgecolor("face")
ax_rightbot.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_rightbot.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_rightbot.set_title(r'\textrm{350 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im4, ax=axs[1, 1], orientation='vertical',shrink=0.8)

cbar = fig.colorbar(im2, ax=axs[1, :], orientation='horizontal')

cbar.set_label(r'-$\overline{U}\overline{uv}/U^3_{hub}$',fontsize=gcafontSize)

plt.savefig(plot_path+prefix+'entrainmentFlux_contour.pdf',dpi=500)

plt.close()





# velocity flux

figwidth  = 19
figheight = 15
lineWidth = 3
textFontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


#fig = plt.figure(2, figsize=(figwidth,figheight))
fig,axs = plt.subplots(2, 2, figsize=(figwidth,figheight))

ax_lefttop = axs[0,0]
ax_righttop = axs[0,1]
ax_leftbot = axs[1,0]
ax_rightbot = axs[1,1]


min_act = np.min(np.min(U_fluc[:,:,500]/np.max(np.max(U_fluc[:,:,500]))))
max_act = np.max(np.max(U_fluc[:,:,500]/np.max(np.max(U_fluc[:,:,500]))))
min_100 = np.min(np.min(u_recons_100[:,:,500]/np.max(np.max(u_recons_100[:,:,500]))))
max_100 = np.max(np.max(u_recons_100[:,:,500]/np.max(np.max(u_recons_100[:,:,500]))))
min_10 = np.min(np.min(u_recons_10[:,:,500]/np.max(np.max(u_recons_10[:,:,500]))))
max_10 = np.max(np.max(u_recons_10[:,:,500]/np.max(np.max(u_recons_10[:,:,500]))))
min_210 = np.min(np.min(u_recons_210[:,:,500]/np.max(np.max(u_recons_210[:,:,500]))))
max_210 = np.max(np.max(u_recons_210[:,:,500]/np.max(np.max(u_recons_210[:,:,500]))))

min_avg = np.mean([min_act,min_10,min_100,min_210])
max_avg = np.mean([max_act,max_10,max_100,max_210])


im1 = ax_lefttop.contourf(X_norm,Y_norm,U_fluc[:,:,500]/np.max(np.max(U_fluc[:,:,500])),np.linspace(min_act,max_act,100),cmap='jet')
for c in im1.collections:
    c.set_edgecolor("face")
ax_lefttop.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_lefttop.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_lefttop.set_title(r'\textrm{Instantaneous}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im1, ax=axs[0, 0], orientation='vertical',shrink=0.8)


im2 = ax_righttop.contourf(X_norm,Y_norm,u_recons_10[:,:,500]/np.max(np.max(u_recons_10[:,:,500])),np.linspace(min_avg,max_avg,100),cmap='jet')
for c in im2.collections:
    c.set_edgecolor("face")
ax_righttop.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_righttop.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_righttop.set_title(r'\textrm{10 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im2, ax=axs[0, 1], orientation='vertical',shrink=0.8)



im3 = ax_leftbot.contourf(X_norm,Y_norm,u_recons_100[:,:,500]/np.max(np.max(u_recons_100[:,:,500])),np.linspace(min_act,max_act,100),cmap='jet')
for c in im3.collections:
    c.set_edgecolor("face")
ax_leftbot.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_leftbot.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_leftbot.set_title(r'\textrm{100 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)
#cbar = fig.colorbar(im3, ax=axs[1, 0], orientation='vertical',shrink=0.8)


im4 = ax_rightbot.contourf(X_norm,Y_norm,u_recons_210[:,:,500]/np.max(np.max(u_recons_210[:,:,500])),np.linspace(min_avg,max_act,100),cmap='jet')
for c in im4.collections:
    c.set_edgecolor("face")
ax_rightbot.set_xlabel(r'$x/D$',fontsize=gcafontSize)
ax_rightbot.set_ylabel(r'$(y-H)/D$',fontsize=gcafontSize)
ax_rightbot.set_title(r'\textrm{210 modes}',fontsize=gcafontSize)
plt.xticks(fontsize=gcafontSize)
plt.yticks(fontsize=gcafontSize)


cbar = fig.colorbar(im2, ax=axs[1, :], orientation='horizontal')

cbar.set_label(r'$u/u_{max}$',fontsize=gcafontSize)

plt.savefig(plot_path+prefix+'_contour.pdf',dpi=500)

plt.close


print("Done for "+prefix)

