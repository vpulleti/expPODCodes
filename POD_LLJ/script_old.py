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


# base case
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/inlet/Inlet_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/inlet/"

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_b = tke.sum(1)/tstp
# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

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
uv_b = uv_b/vel
tke_spaavg_b = tke_spaavg_b/(vel**2)

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrL = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

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
uv_arrL = uv_arrL/vel
tke_spaavg_arrL = tke_spaavg_arrL/(vel**2)

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrM = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

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
uv_arrM = uv_arrM/vel
tke_spaavg_arrM = tke_spaavg_arrM/(vel**2)

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_arrH = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

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
uv_arrH = uv_arrH/vel
tke_spaavg_arrH = tke_spaavg_arrH/(vel**2)

# Single_2D_L
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_L/single_2D_L_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_2D_L/"

prefix = "single_2D_L"

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DL = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DL = rho_uu.mean(1)
rho_vv_s2DL = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DL = urms_s2DL/vel
vrms_s2DL = vrms_s2DL/vel
uv_s2DL = uv_s2DL/vel
tke_spaavg_s2DL = tke_spaavg_s2DL/(vel**2)

#single_2D_M

bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_M/single_2D_M_"
ascii_path = "/scratch/brown/vpulleti/LLJ_PIV/single_2D_M/"

prefix = "single_2D_M"

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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DM = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DM = rho_uu.mean(1)
rho_vv_s2DM = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DM = urms_s2DM/vel
vrms_s2DM = vrms_s2DM/vel
uv_s2DM = uv_s2DM/vel
tke_spaavg_s2DM = tke_spaavg_s2DM/(vel**2)

#single_@D_H
bin_path = "/scratch/brown/vpulleti/LLJ_PIV/binary/single_2D_M/single_2D_H_"
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


out = code_functions.tke(U_fluc,V_fluc)

tke = out['tke']

tke_spaavg_s2DH = tke.sum(1)/tstp

# two point correlations
'''
out = code_functions.two_point_corrl(loc_u,loc_v,urms,vrms,,U_fluc,V_fluc)

rho_uu = out['rhouu']
rho_vv = out['rhovv']

del out

rho_uu_s2DH = rho_uu.mean(1)
rho_vv_s2DH = rho_vv.mean(1)
'''
# Normalizing
out = code_functions.normalize(prefix)

length = out['length']
vel = out['vel']
D = out['D']

urms_s2DH = urms_s2DH/vel
vrms_s2DH = vrms_s2DH/vel
uv_s2DH = uv_s2DH/vel
tke_spaavg_s2DH = tke_spaavg_s2DH/(vel**2)








Y = 0.001*(Y-200)/D






# plotting

U_mean = U_mean.mean(1)
V_mean = V_mean.mean(1)

cmin = np.min(U_mean)
cmax = np.max(U_mean)

figure_name = bin_path+"u_mean.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(U_mean,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$\overline{U}$",fontsize=textFontSize)
ax.set_title(r"$U_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


cmin = np.min(V_mean)
cmax = np.max(V_mean)

figure_name = bin_path+"v_mean.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(V_mean,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$\overline{V}$",fontsize=textFontSize)
ax.set_title(r"$V_{mean}$",fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()



cmin = np.min(urms)
cmax = np.max(urms)

figure_name = bin_path+"u_rms.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(urms,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$u_{rms}$",fontsize=textFontSize)
ax.set_title('urms',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


cmin = np.min(vrms)
cmax = np.max(vrms)

figure_name = bin_path+"v_rms.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(vrms,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$v_{rms}$",fontsize=textFontSize)
ax.set_title('vrms',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()


cmin = np.min(uv)
cmax = np.max(uv)

figure_name = bin_path+"uv.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(uv,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$-\overline{uv}$",fontsize=textFontSize)
ax.set_title('uv',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()

cmin = np.min(tke_spaavg)
cmax = np.max(tke_spaavg)

figure_name = bin_path+"tke.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(tke_spaavg,Y.mean(1),linewidth=lineWidth)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlim(cmin,cmax)
ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$TKE$",fontsize=textFontSize)
ax.set_title('TKE',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()



figure_name = bin_path+"periodicity_verification.pdf"
figwidth  = 18
figheight = 10
lineWidth = 3
textFontSize = 28
gcafontSize = 20

fig = plt.figure(0,figsize=(figwidth,figheight))
ax  = plt.gca()
plt.plot(U[100,:,100])
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
#ax.set_xlim(cmin,cmax)
#ax.set_ylim(np.min(Y),np.max(Y))
ax.set_ylabel(r"$y$",fontsize=textFontSize)
ax.set_xlabel(r"$-\overline{uv}$",fontsize=textFontSize)
ax.set_title('uv',fontsize=textFontSize)
print("Saving figure:"+figure_name)
plt.tight_layout()
plt.savefig(figure_name)
plt.close()



