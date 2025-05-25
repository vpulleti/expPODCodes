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
import scipy as spy
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path
import matplotlib.pyplot as plt
import code_functions
import POD_functions
from pdb import set_trace

mpl.rc('text', usetex = True)
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})


nx = 500
ny = 334
num_tstp = 4000
loc_u = np.array([197,10])
loc_v = np.array([197,10])
'''
# base case
bin_path = "/scratch/bell/vpulleti/flatPlate_exp/UIUC_BL/sm_HD_25Hz/binary/"
ascii_path = "/scratch/bell/vpulleti/flatPlate_exp/UIUC_BL/sm_HD_25Hz/img/"
plot_path = "/scratch/bell/vpulleti/flatPlate_exp/UIUC_BL/sm_HD_25Hz/plots/"
prefix = "B0"
bin_name ="/scratch/bell/vpulleti/flatPlate_exp/UIUC_BL/sm_HD_25Hz/binary/check/" # path name can be changed to directory where you want to store the binary files
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
out = code_functions.flatPlate_creation(ny,nx,num_tstp,prefix,ascii_path,bin_path,bin_name)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out
#masking

x = x[20:-21,20:-21]
y = y[20:-21,20:-21]
u = u[20:-21,20:-21]
v = v[20:-21,20:-21]

uinst = np.flip(np.flip(u,axis=0),axis=1)
vinst = np.flip(np.flip(v,axis=0),axis=1)

U = np.mean(uinst,axis=2)
V = np.mean(vinst,axis=2)

num_tstp = np.shape(uinst)[2]

ufluc = np.zeros(np.shape(uinst))
vfluc = np.zeros(np.shape(vinst))

for time in range(np.shape(ufluc)[2]):
  ufluc[:,:,time] = uinst[:,:,time]-U
  vfluc[:,:,time] = vinst[:,:,time]-V

uvMean = np.mean(ufluc*vfluc,axis=2)
TKEMean  = np.mean(0.5*(ufluc**2+vfluc**2),axis=2)

uhat_fluc = np.zeros([1024,np.shape(ufluc)[1],np.shape(ufluc)[2]],dtype=complex)
vhat_fluc = np.zeros([1024,np.shape(vfluc)[1],np.shape(vfluc)[2]],dtype=complex)

for tstp in range(num_tstp):
  uhat_fluc[:,:,tstp] = np.fft.fft(ufluc[:,:,tstp],n=1024,axis=0,norm="forward")
  vhat_fluc[:,:,tstp] = np.fft.fft(vfluc[:,:,tstp],n=1024,axis=0,norm="forward")


kx = 2*np.pi*np.linspace(0,np.size(uhat_fluc[:,10,10]),np.size(uhat_fluc[:,10,10]))/(x[-1,10]-x[0,10])

lambx = 2*np.pi/kx

E_uu = np.mean(uhat_fluc*np.conjugate(uhat_fluc),axis=2)/np.max(U[10,:])**2
E_vv = np.mean(vhat_fluc*np.conjugate(vhat_fluc),axis=2)/np.max(U[10,:])**2

phi_uu = np.zeros(np.shape(E_uu))
phi_vv = np.zeros(np.shape(E_vv))


for j in range(np.shape(E_uu)[1]):
  phi_uu[:,j] = kx*E_uu[:,j]
  phi_vv[:,j] = kx*E_vv[:,j]

with open(plot_path+'premultspectra_uu.pickle','wb') as file:
  pickle.dump({'phi_uu':phi_uu,'kx':kx,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)

with open(plot_path+'premultspectra_vv.pickle','wb') as file:
  pickle.dump({'phi_vv':phi_vv,'kx':kx,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)

with open(plot_path+'Umean.pickle','wb') as file:
  pickle.dump({'U':U,'x':x,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)

with open(plot_path+'Vmean.pickle','wb') as file:
  pickle.dump({'V':V,'x':x,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)
 
with open(plot_path+'uvmean.pickle','wb') as file:
  pickle.dump({'uv':uvMean,'x':x,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)

with open(plot_path+'TKEmean.pickle','wb') as file:
  pickle.dump({'TKE':TKEMean,'x':x,'y':y},file,protocol=pickle.HIGHEST_PROTOCOL)


'''

set_trace()

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

if prefix=="BL_w0_2D" or prefix=="BL_w1_2D" or prefix=="BL_w2_2D" or prefix=="BL_w3_2D":
  X_norm = X_new/D

elif prefix=="BL_w0_7D" or prefix=="BL_w1_7D" or prefix=="BL_w2_7D" or prefix=="BL_w3_7D":
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

mod_reqd = (np.abs(cum_eigvals-80)).argmin()

figwidth  = 12
figheight = 10
lineWidth = 3
textFontSize = 32
gcafontSize = 28

fig = plt.figure(1, figsize=(figwidth,figheight))
ax1 = plt.axes()
ax2 = ax1.twinx()
ax1.bar(np.arange(620)+1,eigvals[:620],width=10)
ax2.plot(np.arange(620)+1,cum_eigvals[:620],'r',linewidth=2)
ax1.set_xlabel(r'\textrm{Modes}',fontsize=textFontSize)
ax1.set_ylabel(r'\textrm{Modal energy}',fontsize=textFontSize)
ax2.set_ylabel(r'\textrm{Cumulative energy} \%',fontsize=textFontSize)
plt.xlim([1,620])
plt.setp(ax1.get_xticklabels(), fontsize=gcafontSize)
plt.setp(ax1.get_yticklabels(),fontsize=gcafontSize)
plt.setp(ax2.get_yticklabels(),fontsize=gcafontSize)
plt.savefig(plot_path+prefix+'_eigvals.pdf',dpi=500)
print("save file "+prefix+'_eigvals.pdf')
plt.close()

u_recons_LSM = POD_functions.reconst(timcoeff_u,phi_u,0,mod_reqd,Nx,Ny)
v_recons_LSM = POD_functions.reconst(timcoeff_v,phi_v,0,mod_reqd,Nx,Ny)

print("Reconstructed using "+str(mod_reqd)+" modes")

u_recons_SSM = POD_functions.reconst(timcoeff_u,phi_u,mod_reqd,num_tstp-1,Nx,Ny)
v_recons_SSM = POD_functions.reconst(timcoeff_v,phi_v,mod_reqd,num_tstp-1,Nx,Ny)

print("Reconstructed using "+str(num_tstp-mod_reqd)+" modes")

TKE_LSM = 0.5*(np.mean(u_recons_LSM**2+v_recons_LSM**2,axis=2))/U_hub**2

TKE_SSM = 0.5*(np.mean(u_recons_SSM**2+v_recons_SSM**2,axis=2))/U_hub**2

TKE_real = 0.5*(np.mean(U_fluc**2+V_fluc**2,axis=2))/U_hub**2

uv_LSM = np.zeros(np.shape(U_new))
uv_SSM = np.zeros(np.shape(U_new))
vv_LSM = np.zeros(np.shape(U_new))
vv_SSM = np.zeros(np.shape(U_new))

for i in range(num_tstp):
  uv_LSM[:,:,i] = -u_recons_LSM[:,:,i]*v_recons_LSM[:,:,i]
  uv_SSM[:,:,i] = -u_recons_SSM[:,:,i]*v_recons_SSM[:,:,i]  
  vv_LSM[:,:,i] = v_recons_LSM[:,:,i]**2
  vv_SSM[:,:,i] = v_recons_SSM[:,:,i]**2

uv_LSMMean = np.mean(uv_LSM,axis=2)
uv_SSMMean = np.mean(uv_SSM,axis=2)

vv_LSMMean = np.mean(vv_LSM,axis=2)
vv_SSMMean = np.mean(vv_SSM,axis=2)

Uuv_LSMMean = UMean_2D*uv_LSMMean/U_hub**3
Uuv_SSMMean = UMean_2D*uv_SSMMean/U_hub**3

uv_LSMMean = uv_LSMMean/U_hub**2
uv_SSMMean = uv_SSMMean/U_hub**2

vv_LSMMean = vv_LSMMean/U_hub**2
vv_SSMMean = vv_SSMMean/U_hub**2

uinst_LSM = np.zeros(np.shape(u_recons_LSM))
vinst_LSM = np.zeros(np.shape(v_recons_LSM))

uinst_SSM = np.zeros(np.shape(u_recons_SSM))
vinst_SSM = np.zeros(np.shape(v_recons_SSM))

for snap in range(num_tstp):
  uinst_LSM[:,:,snap] = u_recons_LSM[:,:,snap]+UMean_2D
  vinst_LSM[:,:,snap] = v_recons_LSM[:,:,snap]+Vmean_2D

  uinst_SSM[:,:,snap] = u_recons_SSM[:,:,snap]+UMean_2D
  vinst_SSM[:,:,snap] = v_recons_SSM[:,:,snap]+Vmean_2D

print("Calculating the lambda2 for LSM")
lambdaci_LSM = code_functions.get_lambda2(u_recons_LSM,v_recons_LSM,X_new,Y_new,lmts)

print("Calculating the lambda2 for SSM")
lambdaci_SSM = code_functions.get_lambda2(u_recons_SSM,v_recons_SSM,X_new,Y_new,lmts)

print("Calculating the vorticity for real")
lambdaci_real = code_functions.get_lambda2(U_new,V_new,X_new,Y_new,lmts)


data = {}

data['Uuv_real'] = Uuv_realMean
data['uv_real'] = uv_realMean
data['vv_real'] = vv_realMean
data['Uuv_LSM'] = Uuv_LSMMean
data['uv_LSM'] = uv_LSMMean
data['vv_LSM'] = vv_LSMMean
data['Uuv_SSM'] = Uuv_SSMMean
data['uv_SSM'] = uv_SSMMean
data['vv_SSM'] = vv_SSMMean
data['lambda_LSM'] = lambdaci_LSM*avgChord_len/U_hub
data['lambda_SSM'] = lambdaci_SSM*avgChord_len/U_hub
data['lambda_real'] = lambdaci_real*avgChord_len/U_hub
#data['Uuvu_ssmMean'] = Uuv_ssmMean
data['TKE_LSM'] = TKE_LSM
data['TKE_SSM'] = TKE_SSM
data['TKE_real'] = TKE_real
data['X_norm'] = X_norm
data['Y_norm'] = Y_norm
data['uinst_LSM'] = uinst_LSM
data['uinst_SSM'] = uinst_SSM
data['vinst_LSM'] = vinst_LSM
data['vinst_SSM'] = vinst_SSM
data['uinst_real'] = U_new
data['vinst_real'] = V_new

with open(bin_path+prefix+"data.pickle",'wb') as file:
  pickle.dump(data,file)

print("Dumped the data for "+prefix)


'''