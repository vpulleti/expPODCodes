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
import code_functions

ny = 314
nx = 207
tstp = 1000
loc_u = np.array([197,10])
loc_v = np.array([197,10])

# base case
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/inlet/Inlet_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/inlet/"
plot_path = "/scratch/brown/vpulleti/LLJ_PIV/"
prefix = "Inlet"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_b = urms.mean(1)
vrms_b = vrms.mean(1)
uv_b = uv.mean(1)
U_mean_b = U_mean.mean(1)
V_mean_b = V_mean.mean(1)

kflux_b = uv*U_mean

out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_b = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_b = prod.mean(1)

# two point correlations

#keyboard()

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_b = rho_uu
rho_vv_b = rho_vv
'''
rho_uu_b = rho_uu.mean(1)
rho_vv_b = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_b = urms_b/vel
vrms_b = vrms_b/vel
uv_b = uv_b/vel**2
tke_spaavg_b = tke_spaavg_b/(vel**2)
U_mean_b = U_mean_b/vel
V_mean_b = V_mean_b/vel
prod_b = 1000*D*prod_b/vel**3
kflux_b = kflux_b/vel**3
# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_b = out['spectra']
ku_b  = out['k']
#keyboard()



# change in number of timesteps
tstp = 2000

# array_L

bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/array_L/array_L_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/array_L/"

prefix = "array_7D_L"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]
#keyboard()
# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_arrL = urms.mean(1)
vrms_arrL = vrms.mean(1)
uv_arrL = uv.mean(1)
U_mean_arrL = U_mean.mean(1)
V_mean_arrL = V_mean.mean(1)

kflux_arrL = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrL = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_arrL = prod.mean(1)

# two point correlations
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_arrL = rho_uu
rho_vv_arrL = rho_vv
'''
rho_uu_arrL = rho_uu.mean(1)
rho_vv_arrL = rho_vv.mean(1)

'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_arrL = urms_arrL/vel
vrms_arrL = vrms_arrL/vel
uv_arrL = uv_arrL/vel**2
tke_spaavg_arrL = tke_spaavg_arrL/(vel**2)
U_mean_arrL = U_mean_arrL/vel
V_mean_arrL = V_mean_arrL/vel
prod_arrL = 1000*D*prod_arrL/vel**3
kflux_arrL = kflux_arrL/vel**3

# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_arrL = out['spectra']
ku_arrL  = out['k']


# array_M

bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/array_M/array_M_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/array_M/"

prefix = "array_7D_M"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_arrM = urms.mean(1)
vrms_arrM = vrms.mean(1)
uv_arrM = uv.mean(1)
U_mean_arrM = U_mean.mean(1)
V_mean_arrM = V_mean.mean(1)

kflux_arrM = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrM = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_arrM = prod.mean(1)


# two point correlations

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out
rho_uu_arrM = rho_uu
rho_vv_arrM = rho_vv

'''
rho_uu_arrM = rho_uu.mean(1)
rho_vv_arrM = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_arrM = urms_arrM/vel
vrms_arrM = vrms_arrM/vel
uv_arrM = uv_arrM/vel**2
tke_spaavg_arrM = tke_spaavg_arrM/(vel**2)
U_mean_arrM = U_mean_arrM/vel
V_mean_arrM = V_mean_arrM/vel
prod_arrM = 1000*D*prod_arrM/vel**3
kflux_arrM = kflux_arrM/vel**3

# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_arrM = out['spectra']
ku_arrM  = out['k']


# array_H
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/array_H/array_H_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/array_H/"

prefix = "array_7D_H"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_arrH = urms.mean(1)
vrms_arrH = vrms.mean(1)
uv_arrH = uv.mean(1)
U_mean_arrH = U_mean.mean(1)
V_mean_arrH = V_mean.mean(1)

kflux_arrH = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrH = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_arrH = prod.mean(1)

# two point correlations
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_arrH = rho_uu
rho_vv_arrH = rho_vv
'''
rho_uu_arrH = rho_uu.mean(1)
rho_vv_arrH = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_arrH = urms_arrH/vel
vrms_arrH = vrms_arrH/vel
uv_arrH = uv_arrH/vel**2
tke_spaavg_arrH = tke_spaavg_arrH/(vel**2)
U_mean_arrH = U_mean_arrH/vel
V_mean_arrH = V_mean_arrH/vel
prod_arrH = 1000*D*prod_arrH/vel**3
kflux_arrH = kflux_arrH/vel**3

# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_arrH = out['spectra']
ku_arrH  = out['k']


# Single_2D_L
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_L/single_2D_L_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_2D_L/"

prefix = "Single_2D_L"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s2DL = urms.mean(1)
vrms_s2DL = vrms.mean(1)
uv_s2DL = uv.mean(1)
U_mean_s2DL = U_mean.mean(1)
V_mean_s2DL = V_mean.mean(1)

kflux_s2DL = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DL = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s2DL = prod.mean(1)


# two point correlations

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DL = rho_uu
rho_vv_s2DL = rho_vv

#Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DL = urms_s2DL/vel
vrms_s2DL = vrms_s2DL/vel
uv_s2DL = uv_s2DL/vel**2
tke_spaavg_s2DL = tke_spaavg_s2DL/(vel**2)
U_mean_s2DL = U_mean_s2DL/vel
V_mean_s2DL = V_mean_s2DL/vel
prod_s2DL = 1000*D*prod_s2DL/vel**3
kflux_s2DL = kflux_s2DL/vel**3
# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s2DL = out['spectra']
ku_s2DL  = out['k']

#single_2D_M

bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_M/single_2D_M_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_2D_M/"

prefix = "Single_2D_M"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s2DM = urms.mean(1)
vrms_s2DM = vrms.mean(1)
uv_s2DM = uv.mean(1)
U_mean_s2DM = U_mean.mean(1)
V_mean_s2DM = V_mean.mean(1)

kflux_s2DM = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DM = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s2DM = prod.mean(1)

# two point correlations

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DM = rho_uu
rho_vv_s2DM = rho_vv

# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DM = urms_s2DM/vel
vrms_s2DM = vrms_s2DM/vel
uv_s2DM = uv_s2DM/vel**2
tke_spaavg_s2DM = tke_spaavg_s2DM/(vel**2)
U_mean_s2DM = U_mean_s2DM/vel
V_mean_s2DM = V_mean_s2DM/vel
prod_s2DM = 1000*D*prod_s2DM/vel**3
kflux_s2DM =kflux_s2DM/vel**3
# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s2DM = out['spectra']
ku_s2DM  = out['k']

#single_2D_H
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_H/single_2D_H_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_2D_H/"

prefix = "single_2D_H"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s2DH = urms.mean(1)
vrms_s2DH = vrms.mean(1)
uv_s2DH = uv.mean(1)
U_mean_s2DH = U_mean.mean(1)
V_mean_s2DH = V_mean.mean(1)
kflux_s2DH = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DH = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s2DH = prod.mean(1)

# two point correlations
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DH = rho_uu
rho_vv_s2DH = rho_vv

# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DH = urms_s2DH/vel
vrms_s2DH = vrms_s2DH/vel
uv_s2DH = uv_s2DH/vel**2
tke_spaavg_s2DH = tke_spaavg_s2DH/(vel**2)
U_mean_s2DH = U_mean_s2DH/vel
V_mean_s2DH = V_mean_s2DH/vel
prod_s2DH = 1000*D*prod_s2DH/vel**3
kflux_s2DH = kflux_s2DH/vel**3

# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s2DH = out['spectra']
ku_s2DH  = out['k']


# Single_7D_L
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_7D_L/single_7D_L_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_7D_L/"

prefix = "single_7D_L"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s7DL = urms.mean(1)
vrms_s7DL = vrms.mean(1)
uv_s7DL = uv.mean(1)
U_mean_s7DL = U_mean.mean(1)
V_mean_s7DL = V_mean.mean(1)
kflux_s7DL = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s7DL = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s7DL = prod.mean(1)

# two point correlations
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s7DL = rho_uu
rho_vv_s7DL = rho_vv

# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s7DL = urms_s7DL/vel
vrms_s7DL = vrms_s7DL/vel
uv_s7DL = uv_s7DL/vel**2
tke_spaavg_s7DL = tke_spaavg_s7DL/(vel**2)
U_mean_s7DL = U_mean_s7DL/vel
V_mean_s7DL = V_mean_s7DL/vel
prod_s7DL = 1000*D*prod_s7DL/vel**3
kflux_s7DL = kflux_s7DL/vel**3
# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s7DL = out['spectra']
ku_s7DL  = out['k']


#single_7D_M

bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_7D_M/single_7D_M_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_7D_M/"

prefix = "single_7D_M"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s7DM = urms.mean(1)
vrms_s7DM = vrms.mean(1)
uv_s7DM = uv.mean(1)
U_mean_s7DM = U_mean.mean(1)
V_mean_s7DM = V_mean.mean(1)

kflux_s7DM = uv*U_mean

out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s7DM = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s7DM = prod.mean(1)


# two point correlations

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s7DM = rho_uu
rho_vv_s7DM = rho_vv

# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s7DM = urms_s7DM/vel
vrms_s7DM = vrms_s7DM/vel
uv_s7DM = uv_s7DM/vel**2
tke_spaavg_s7DM = tke_spaavg_s7DM/(vel**2)
U_mean_s7DM = U_mean_s7DM/vel
V_mean_s7DM = V_mean_s7DM/vel
prod_s7DM = 1000*D*prod_s7DM/vel**3
kflux_s7DM = kflux_s7DM/vel**3
# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s7DM = out['spectra']
ku_s7DM  = out['k']


#single_7D_H
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_7D_H/single_7D_H_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_7D_H/"

prefix = "Single_7D_H"

# reading the data
out = code_functions.bin_creation(ny,nx,tstp,prefix,ascii_path,bin_path)

x = out['x']
y = out['y']
u = out['u']
v = out['v']

del out

# discarding the corrupted data

X = -1.*x[:,39:165]
Y = y[:,39:165]
U = -1.*u[:,39:165,:]
V = v[:,39:165,:]

# finding the time mean of the data
U_sum = U.sum(2)
V_sum = V.sum(2)

U_mean = U_sum/float(tstp)
V_mean = V_sum/float(tstp)

# finding the fluctuation of the u and v
U_fluc = np.zeros(np.shape(U))
V_fluc = np.zeros(np.shape(V))

for i in range(tstp):
  U_fluc[:,:,i] = U[:,:,i]-U_mean
  V_fluc[:,:,i] = V[:,:,i]-V_mean

# single point statistics
out = code_functions.single_point_stats(U_fluc,V_fluc)

urms = out['urms']
vrms = out['vrms']
uv = out['uv']

del out

urms_s7DH = urms.mean(1)
vrms_s7DH = vrms.mean(1)
uv_s7DH = uv.mean(1)
U_mean_s7DH = U_mean.mean(1)
V_mean_s7DH = V_mean.mean(1)

kflux_s7DH = uv*U_mean


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s7DH = tke.sum(1)/tstp

# production of turbulence
dUdy = code_functions.firstderv(Y,U_mean,0) # 0 for dy and 1 for dx

prod = -uv*dUdy

prod_s7DH = prod.mean(1)


# two point correlations

out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s7DH = rho_uu
rho_vv_s7DH = rho_vv

# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s7DH = urms_s7DH/vel
vrms_s7DH = vrms_s7DH/vel
uv_s7DH = uv_s7DH/vel**2
tke_spaavg_s7DH = tke_spaavg_s7DH/(vel**2)
U_mean_s7DH = U_mean_s7DH/vel
V_mean_s7DH = V_mean_s7DH/vel
prod_s7DH = 1000*D*prod_s7DH/vel**3
kflux_s7DH = kflux_s7DH/vel**3

# spectra
out = code_functions.spectra_1D(U_fluc/vel,X/(D*1000))

Euu_s7DH = out['spectra']
ku_s7DH  = out['k']


Y_rhouu = 0.001*(Y-Y[loc_u[0],loc_u[1]])/D
Y_rhovv = 0.001*(Y-Y[loc_v[0],loc_v[1]])/D
Y = 0.001*(Y-200)/D
X_rhouu = 0.001*(X-X[loc_u[0],loc_u[1]])/D
X_rhovv = 0.001*(X-X[loc_v[0],loc_v[1]])/D
X = 0.001*X/D


# Writing data to a dat files

# coordinations
#out = code_functions.writedata(X[1,:].transpose(),plot_path+"x.dat")
#out = code_functions.writedata(Y[:,1],plot_path+"y.dat")

#mean

#U
out = code_functions.writedata(Y.mean(1),U_mean_b,plot_path+"umean_b.dat")
out = code_functions.writedata(Y.mean(1),U_mean_arrL,plot_path+"umean_arrL.dat")
out = code_functions.writedata(Y.mean(1),U_mean_arrM,plot_path+"umean_arrM.dat")
out = code_functions.writedata(Y.mean(1),U_mean_arrH,plot_path+"umean_arrH.dat")
out = code_functions.writedata(Y.mean(1),U_mean_s2DL,plot_path+"umean_s2DL.dat")
out = code_functions.writedata(Y.mean(1),U_mean_s2DM,plot_path+"umean_s2DM.dat")
out = code_functions.writedata(Y.mean(1),U_mean_s2DH,plot_path+"umean_s2DH.dat")

out = code_functions.writedata(Y.mean(1),U_mean_s7DL,plot_path+"umean_s7DL.dat")
out = code_functions.writedata(Y.mean(1),U_mean_s7DM,plot_path+"umean_s7DM.dat")
out = code_functions.writedata(Y.mean(1),U_mean_s7DH,plot_path+"umean_s7DH.dat")

#V
out = code_functions.writedata(Y.mean(1),V_mean_b,plot_path+"vmean_b.dat")
out = code_functions.writedata(Y.mean(1),V_mean_arrL,plot_path+"vmean_arrL.dat")
out = code_functions.writedata(Y.mean(1),V_mean_arrM,plot_path+"vmean_arrM.dat")
out = code_functions.writedata(Y.mean(1),V_mean_arrH,plot_path+"vmean_arrH.dat")
out = code_functions.writedata(Y.mean(1),V_mean_s2DL,plot_path+"vmean_s2DL.dat")
out = code_functions.writedata(Y.mean(1),V_mean_s2DM,plot_path+"vmean_s2DM.dat")
out = code_functions.writedata(Y.mean(1),V_mean_s2DH,plot_path+"vmean_s2DH.dat")

out = code_functions.writedata(Y.mean(1),V_mean_s7DL,plot_path+"vmean_s7DL.dat")
out = code_functions.writedata(Y.mean(1),V_mean_s7DM,plot_path+"vmean_s7DM.dat")
out = code_functions.writedata(Y.mean(1),V_mean_s7DH,plot_path+"vmean_s7DH.dat")

# RMS
#U
out = code_functions.writedata(Y.mean(1),urms_b,plot_path+"urms_b.dat")
out = code_functions.writedata(Y.mean(1),urms_arrL,plot_path+"urms_arrL.dat")
out = code_functions.writedata(Y.mean(1),urms_arrM,plot_path+"urms_arrM.dat")
out = code_functions.writedata(Y.mean(1),urms_arrH,plot_path+"urms_arrH.dat")
out = code_functions.writedata(Y.mean(1),urms_s2DL,plot_path+"urms_s2DL.dat")
out = code_functions.writedata(Y.mean(1),urms_s2DM,plot_path+"urms_s2DM.dat")
out = code_functions.writedata(Y.mean(1),urms_s2DH,plot_path+"urms_s2DH.dat")

out = code_functions.writedata(Y.mean(1),urms_s7DL,plot_path+"urms_s7DL.dat")
out = code_functions.writedata(Y.mean(1),urms_s7DM,plot_path+"urms_s7DM.dat")
out = code_functions.writedata(Y.mean(1),urms_s7DH,plot_path+"urms_s7DH.dat")

#V
out = code_functions.writedata(Y.mean(1),vrms_b,plot_path+"vrms_b.dat")
out = code_functions.writedata(Y.mean(1),vrms_arrL,plot_path+"vrms_arrL.dat")
out = code_functions.writedata(Y.mean(1),vrms_arrM,plot_path+"vrms_arrM.dat")
out = code_functions.writedata(Y.mean(1),vrms_arrH,plot_path+"vrms_arrH.dat")
out = code_functions.writedata(Y.mean(1),vrms_s2DL,plot_path+"vrms_s2DL.dat")
out = code_functions.writedata(Y.mean(1),vrms_s2DM,plot_path+"vrms_s2DM.dat")
out = code_functions.writedata(Y.mean(1),vrms_s2DH,plot_path+"vrms_s2DH.dat")

out = code_functions.writedata(Y.mean(1),vrms_s7DL,plot_path+"vrms_s7DL.dat")
out = code_functions.writedata(Y.mean(1),vrms_s7DM,plot_path+"vrms_s7DM.dat")
out = code_functions.writedata(Y.mean(1),vrms_s7DH,plot_path+"vrms_s7DH.dat")


#uv

out = code_functions.writedata(Y.mean(1),uv_b,plot_path+"uv_b.dat")
out = code_functions.writedata(Y.mean(1),uv_arrL,plot_path+"uv_arrL.dat")
out = code_functions.writedata(Y.mean(1),uv_arrM,plot_path+"uv_arrM.dat")
out = code_functions.writedata(Y.mean(1),uv_arrH,plot_path+"uv_arrH.dat")
out = code_functions.writedata(Y.mean(1),uv_s2DL,plot_path+"uv_s2DL.dat")
out = code_functions.writedata(Y.mean(1),uv_s2DM,plot_path+"uv_s2DM.dat")
out = code_functions.writedata(Y.mean(1),uv_s2DH,plot_path+"uv_s2DH.dat")

out = code_functions.writedata(Y.mean(1),uv_s7DL,plot_path+"uv_s7DL.dat")
out = code_functions.writedata(Y.mean(1),uv_s7DM,plot_path+"uv_s7DM.dat")
out = code_functions.writedata(Y.mean(1),uv_s7DH,plot_path+"uv_s7DH.dat")


# tke
out = code_functions.writedata(Y.mean(1),tke_spaavg_b,plot_path+"tke_spaavg_b.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_arrL,plot_path+"tke_spaavg_arrL.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_arrM,plot_path+"tke_spaavg_arrM.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_arrH,plot_path+"tke_spaavg_arrH.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_s2DL,plot_path+"tke_spaavg_s2DL.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_s2DM,plot_path+"tke_spaavg_s2DM.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_s2DH,plot_path+"tke_spaavg_s2DH.dat")

out = code_functions.writedata(Y.mean(1),tke_spaavg_s7DL,plot_path+"tke_spaavg_s7DL.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_s7DM,plot_path+"tke_spaavg_s7DM.dat")
out = code_functions.writedata(Y.mean(1),tke_spaavg_s7DH,plot_path+"tke_spaavg_s7DH.dat")


# production
out = code_functions.writedata(Y.mean(1),prod_b,plot_path+"prod_b.dat")
out = code_functions.writedata(Y.mean(1),prod_arrL,plot_path+"prod_arrL.dat")
out = code_functions.writedata(Y.mean(1),prod_arrM,plot_path+"prod_arrM.dat")
out = code_functions.writedata(Y.mean(1),prod_arrH,plot_path+"prod_arrH.dat")
out = code_functions.writedata(Y.mean(1),prod_s2DL,plot_path+"prod_s2DL.dat")
out = code_functions.writedata(Y.mean(1),prod_s2DM,plot_path+"prod_s2DM.dat")
out = code_functions.writedata(Y.mean(1),prod_s2DH,plot_path+"prod_s2DH.dat")

out = code_functions.writedata(Y.mean(1),prod_s7DL,plot_path+"prod_s7DL.dat")
out = code_functions.writedata(Y.mean(1),prod_s7DM,plot_path+"prod_s7DM.dat")
out = code_functions.writedata(Y.mean(1),prod_s7DH,plot_path+"prod_s7DH.dat")


# kflux

out = code_functions.writedata(Y.mean(1),kflux_b.mean(1),plot_path+"kfluxmean_b.dat")
out = code_functions.writedata(Y.mean(1),kflux_arrL.mean(1),plot_path+"kfluxmean_arrL.dat")
out = code_functions.writedata(Y.mean(1),kflux_arrM.mean(1),plot_path+"kfluxmean_arrM.dat")
out = code_functions.writedata(Y.mean(1),kflux_arrH.mean(1),plot_path+"kfluxmean_arrH.dat")
out = code_functions.writedata(Y.mean(1),kflux_s2DL.mean(1),plot_path+"kfluxmean_s2DL.dat")
out = code_functions.writedata(Y.mean(1),kflux_s2DM.mean(1),plot_path+"kfluxmean_s2DM.dat")
out = code_functions.writedata(Y.mean(1),kflux_s2DH.mean(1),plot_path+"kfluxmean_s2DH.dat")

out = code_functions.writedata(Y.mean(1),kflux_s7DL.mean(1),plot_path+"kfluxmean_s7DL.dat")
out = code_functions.writedata(Y.mean(1),kflux_s7DM.mean(1),plot_path+"kfluxmean_s7DM.dat")
out = code_functions.writedata(Y.mean(1),kflux_s7DH.mean(1),plot_path+"kfluxmean_s7DH.dat")



# plotting

U_mean = U_mean.mean(1)
V_mean = V_mean.mean(1)


umin = np.min([np.min(U_mean_b),np.min(U_mean_arrL),np.min(U_mean_arrM),np.min(U_mean_arrH)])-0.1
umax = np.max([np.max(U_mean_b),np.max(U_mean_arrL),np.max(U_mean_arrM),np.max(U_mean_arrH)])+0.3

figure_name = plot_path+"u_mean_array.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(U_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(U_mean_arrL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(U_mean_arrM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(U_mean_arrH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(umin,umax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{U}}{U_{hub}}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$\frac{U_{mean}}{U_{hub}}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


figure_name = plot_path+"u_mean_single_2D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(U_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(U_mean_s2DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(U_mean_s2DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(U_mean_s2DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(umin,umax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{U}}{U_{hub}}$",fontsize=textFontSize)
ax.legend(['base','Single_2D_L','Single_2D_M','Single_2D_H'])
ax.set_title(r"$U_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"u_mean_single_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(U_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(U_mean_s7DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(U_mean_s7DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(U_mean_s7DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(umin,umax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{U}}{U_{hub}}$",fontsize=textFontSize)
ax.legend(['base','Single_7D_L','Single_7D_M','Single_7D_H'])
ax.set_title(r"$U_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()



# vmean
vmin = np.min([np.min(V_mean_b),np.min(V_mean_arrL),np.min(V_mean_arrM),np.min(V_mean_arrH)])-0.1
vmax = np.max([np.max(V_mean_b),np.max(V_mean_arrL),np.max(V_mean_arrM),np.max(V_mean_arrH)])+0.1

figure_name = plot_path+"v_mean_array.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(V_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(V_mean_arrL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(V_mean_arrM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(V_mean_arrH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vmin,vmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{V}}{U_{hub}}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$V_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"v_mean_single_2D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(V_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(V_mean_s2DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(V_mean_s2DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(V_mean_s2DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vmin,vmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{V}}{U_{hub}}$",fontsize=textFontSize)
ax.legend(['base','Single_2D_L','Single_2D_M','Single_2D_H'])
ax.set_title(r"$V_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"v_mean_single_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(V_mean_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(V_mean_s7DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(V_mean_s7DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(V_mean_s7DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vmin,vmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{\overline{V}}{U_hub}$",fontsize=textFontSize)
ax.legend(['base','Single_7D_L','Single_7D_M','Single_7D_H'])
ax.set_title(r"$V_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

uvmin = np.min([np.min(uv_b),np.min(uv_s2DL),np.min(uv_s2DM),np.min(uv_s2DH)])
uvmax = np.min([np.max(uv_b),np.max(uv_s2DL),np.max(uv_s2DM),np.max(uv_s2DH)])

figure_name = plot_path+"uv_single_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(uv_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(uv_s7DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(uv_s7DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(uv_s7DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(uvmin,uvmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_7D_L','Single_7D_M','Single_7D_H'])
ax.set_title(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"uv_single_2D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(uv_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(uv_s2DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(uv_s2DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(uv_s2DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(uvmin,uvmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_2D_L','Single_2D_M','Single_2D_H'])
ax.set_title(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"uv_array_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(uv_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(uv_arrL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(uv_arrM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(uv_arrH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vmin,vmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$\frac{-\overline{uv}}{U_{hub}^2}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

#rms 

urmin = np.min([np.min(urms_b),np.min(urms_s2DL),np.min(urms_s2DM),np.min(urms_s2DH)])-0.1
urmax = np.max([np.max(urms_b),np.max(urms_s2DL),np.max(urms_s2DM),np.max(urms_s2DH)])+0.1

figure_name = plot_path+"urms_single_2D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(urms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(urms_s2DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(urms_s2DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(urms_s2DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(urmin,urmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{u_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_2D_L','Single_2D_M','Single_2D_H'])
ax.set_title(r"$u_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


figure_name = plot_path+"urms_single_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(urms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(urms_s7DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(urms_s7DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(urms_s7DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(urmin,urmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{u_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_7D_L','Single_7D_M','Single_7D_H'])
ax.set_title(r"$u_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


figure_name = plot_path+"urms_array_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(urms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(urms_arrL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(urms_arrM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(urms_arrH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(urmin,urmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{u_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$u_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

#vrms
vrmin = np.min([np.min(vrms_s2DL),np.min(vrms_s2DM),np.min(vrms_s2DH)])
vrmax = np.max([np.max(vrms_s2DL),np.max(vrms_s2DM),np.max(vrms_s2DH)])

figure_name = plot_path+"vrms_single_2D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(vrms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(vrms_s2DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(vrms_s2DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(vrms_s2DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vrmin,vrmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{v_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_2D_L','Single_2D_M','Single_2D_H'])
ax.set_title(r"$v_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"vrms_single_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(vrms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(vrms_s7DL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(vrms_s7DM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(vrms_s7DH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vrmin,vrmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{v_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','Single_7D_L','Single_7D_M','Single_7D_H'])
ax.set_title(r"$v_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

figure_name = plot_path+"vrms_array_7D.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(vrms_b,Y.mean(1),'k',linewidth=lineWidth)
plt.plot(vrms_arrL,Y.mean(1),'r',linewidth=lineWidth)
plt.plot(vrms_arrM,Y.mean(1),'b',linewidth=lineWidth)
plt.plot(vrms_arrH,Y.mean(1),'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(vrmin,vrmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{v_{rms}}{U_{hub}^2}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$v_{rms}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

# rho_uu

rho_uu_b[np.isnan(rho_uu_b)]=0.0
rho_vv_b[np.isnan(rho_vv_b)]=0.0
rho_uu_arrL[np.isnan(rho_uu_arrL)]=0.0
rho_vv_arrL[np.isnan(rho_vv_arrL)]=0.0
rho_uu_arrM[np.isnan(rho_uu_arrM)]=0.0
rho_vv_arrM[np.isnan(rho_vv_arrM)]=0.0
rho_uu_arrH[np.isnan(rho_uu_arrH)]=0.0
rho_vv_arrH[np.isnan(rho_vv_arrH)]=0.0
rho_uu_s2DL[np.isnan(rho_uu_s2DL)]=0.0
rho_uu_s2DM[np.isnan(rho_uu_s2DM)]=0.0
rho_uu_s2DH[np.isnan(rho_uu_s2DH)]=0.0
rho_uu_s7DL[np.isnan(rho_uu_s7DL)]=0.0
rho_uu_s7DM[np.isnan(rho_uu_s7DM)]=0.0
rho_uu_s7DH[np.isnan(rho_uu_s7DH)]=0.0

rho_vv_s2DL[np.isnan(rho_vv_s2DL)]=0.0
rho_vv_s2DM[np.isnan(rho_vv_s2DM)]=0.0
rho_vv_s2DH[np.isnan(rho_vv_s2DH)]=0.0
rho_vv_s7DL[np.isnan(rho_vv_s7DL)]=0.0
rho_vv_s7DM[np.isnan(rho_vv_s7DM)]=0.0
rho_vv_s7DH[np.isnan(rho_vv_s7DH)]=0.0


'''
rho_uu_b = rho_uu_b.mean(1)
rho_vv_b = rho_vv_b.mean(1)
rho_uu_arrL=rho_uu_arrL.mean(1)
rho_vv_arrL=rho_vv_arrL.mean(1)
rho_uu_arrM=rho_vv_arrM.mean(1)
rho_vv_arrM=rho_vv_arrM.mean(1)
rho_uu_arrH=rho_vv_arrH.mean(1)
rho_vv_arrH=rho_vv_arrH.mean(1)

rho_uu_b = rho_uu_b[:,1]
rho_vv_b = rho_vv_b[:,1]
rho_uu_arrL=rho_uu_arrL[:,1]
rho_vv_arrL=rho_vv_arrL[:,1]
rho_uu_arrM=rho_vv_arrM[:,1]
rho_vv_arrM=rho_vv_arrM[:,1]
rho_uu_arrH=rho_vv_arrH[:,1]
rho_vv_arrH=rho_vv_arrH[:,1]
'''
# writing correlations
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_b,plot_path+"Hrho_uu_b.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_b,plot_path+"Hrho_vv_b.dat")

out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_arrL,plot_path+"Hrho_uu_arrL.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_arrM,plot_path+"Hrho_uu_arrM.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_arrH,plot_path+"Hrho_uu_arrH.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_arrL,plot_path+"Hrho_vv_arrL.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_arrM,plot_path+"Hrho_vv_arrM.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_arrH,plot_path+"Hrho_vv_arrH.dat")

out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s2DL,plot_path+"Hrho_uu_s2DL.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s2DM,plot_path+"Hrho_uu_s2DM.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s2DH,plot_path+"Hrho_uu_s2DH.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s2DL,plot_path+"Hrho_vv_s2DL.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s2DM,plot_path+"Hrho_vv_s2DM.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s2DH,plot_path+"Hrho_vv_s2DH.dat")


out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s7DL,plot_path+"Hrho_uu_s7DL.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s7DM,plot_path+"Hrho_uu_s7DM.dat")
out = code_functions.contrwritedata(X_rhouu,Y_rhouu,rho_uu_s7DH,plot_path+"Hrho_uu_s7DH.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s7DL,plot_path+"Hrho_vv_s7DL.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s7DM,plot_path+"Hrho_vv_s7DM.dat")
out = code_functions.contrwritedata(X_rhovv,Y_rhovv,rho_vv_s7DH,plot_path+"Hrho_vv_s7DH.dat")


wave_num = np.zeros(np.shape(Euu_b))

sh = np.shape(Euu_b)
for count in range(sh[0]):
  wave_num[count,:] = ku_b

#keyboard()
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_b[:,1:65],plot_path+"Euu_b.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_arrL[:,1:65],plot_path+"Euu_arrL.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_arrM[:,1:65],plot_path+"Euu_arrM.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_arrH[:,1:65],plot_path+"Euu_arrH.dat")

out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s2DL[:,1:65],plot_path+"Euu_s2DL.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s2DM[:,1:65],plot_path+"Euu_s2DM.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s2DH[:,1:65],plot_path+"Euu_s2DH.dat")

out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s7DL[:,1:65],plot_path+"Euu_s7DL.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s7DM[:,1:65],plot_path+"Euu_s7DM.dat")
out = code_functions.contrwritedata(wave_num[:,:64],Y[:,:64],wave_num[:,:64]*Euu_s7DH[:,1:65],plot_path+"Euu_s7DH.dat")

# Energy flux
#keyboard()

out = code_functions.contrwritedata(X,Y,kflux_b,plot_path+"kflux_b.dat")
out = code_functions.contrwritedata(X,Y,kflux_arrL,plot_path+"kflux_arrL.dat")
out = code_functions.contrwritedata(X,Y,kflux_arrM,plot_path+"kflux_arrM.dat")
out = code_functions.contrwritedata(X,Y,kflux_arrH,plot_path+"kflux_arrH.dat")

out = code_functions.contrwritedata(X,Y,kflux_s2DL,plot_path+"kflux_s2DL.dat")
out = code_functions.contrwritedata(X,Y,kflux_s2DM,plot_path+"kflux_s2DM.dat")
out = code_functions.contrwritedata(X,Y,kflux_s2DH,plot_path+"kflux_s2DH.dat")

out = code_functions.contrwritedata(X,Y,kflux_s7DL,plot_path+"kflux_s7DL.dat")
out = code_functions.contrwritedata(X,Y,kflux_s7DM,plot_path+"kflux_s7DM.dat")
out = code_functions.contrwritedata(X,Y,kflux_s7DH,plot_path+"kflux_s7DH.dat")


levels = np.linspace(-0.1,1,8)

figure_name = plot_path+"rho_uu.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
CS=plt.contourf(X_rhouu,Y_rhouu,rho_uu_b,levels,cmap='jet',linewidths=1.5)
plt.colorbar()
#ax.clabel(CS,levels,inline=1,fmt='%1.1f',fontsize=14)
#plt.plot(Y_rhouu.mean(1),rho_uu_b,'k',linewidth=lineWidth)
#plt.plot(Y_rhouu.mean(1),rho_uu_arrL,'r',linewidth=lineWidth)
#plt.plot(Y_rhouu.mean(1),rho_uu_arrM,'b',linewidth=lineWidth)
#plt.plot(Y_rhouu.mean(1),rho_uu_arrH,'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
#ax.set_ylim(-0.1,1)
#ax.set_xlim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{x}{D}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$\rho_{uu}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


figure_name = plot_path+"rho_vv.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
#plt.contourf(X,Y,rho_vv_b,cmap='jet',vmin=np.min(rho_vv_b),vmax=np.max(rho_vv_b))
#plt.colorbar()
plt.plot(Y_rhovv.mean(1),rho_vv_b,'k',linewidth=lineWidth)
plt.plot(Y_rhovv.mean(1),rho_vv_arrL,'r',linewidth=lineWidth)
plt.plot(Y_rhovv.mean(1),rho_vv_arrM,'b',linewidth=lineWidth)
plt.plot(Y_rhovv.mean(1),rho_vv_arrH,'g',linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
#ax.set_ylim(-0.1,1)
#ax.set_xlim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$\rho_{vv}$",fontsize=textFontSize)
ax.set_xlabel(r"$\frac{y}{D}$",fontsize=textFontSize)
ax.legend(['base','array_L','array_M','array_H'])
ax.set_title(r"$\rho_{vv}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


#keyboard()
'''
Y_spectra = 

figure_name = plot_path+"spectra_inlet.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.contourf(
'''

