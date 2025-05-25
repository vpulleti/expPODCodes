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


ascii_path ='/scratch/brown/vpulleti/LLJ_PIV/base_case/'
bin_path = '/scratch/brown/vpulleti/LLJ_PIV/binary/base_case/'
plot_path = '/scratch/brown/vpulleti/LLJ_PIV/base_case/data/'

ny = 251
nx = 133
with open(ascii_path+'notree.dat','r') as file:
  lines = file.readlines()
  lines_up = lines[4:33387]
  lines_2D = lines[33389:66772]
  lines_4D = lines[66774:100157]
  lines_6D = lines[100159:]
#  keyboard()
  x_up = np.zeros(len(lines_up))
  x_2D = np.zeros(len(lines_2D))
  x_4D = np.zeros(len(lines_4D))
  x_6D = np.zeros(len(lines_6D))
  
  y_up = np.zeros(len(lines_up))
  y_2D = np.zeros(len(lines_2D))
  y_4D = np.zeros(len(lines_4D))
  y_6D = np.zeros(len(lines_6D))

  y_up = np.zeros(len(lines_up))
  y_2D = np.zeros(len(lines_2D))
  y_4D = np.zeros(len(lines_4D))
  y_6D = np.zeros(len(lines_6D))

  y_up = np.zeros(len(lines_up))
  y_2D = np.zeros(len(lines_2D))
  y_4D = np.zeros(len(lines_4D))
  y_6D = np.zeros(len(lines_6D))

  u_up = np.zeros(len(lines_up))
  u_2D = np.zeros(len(lines_2D))
  u_4D = np.zeros(len(lines_4D))
  u_6D = np.zeros(len(lines_6D))
  
  v_up = np.zeros(len(lines_up))
  v_2D = np.zeros(len(lines_2D))
  v_4D = np.zeros(len(lines_4D))
  v_6D = np.zeros(len(lines_6D))

  tke_up = np.zeros(len(lines_up))
  tke_2D = np.zeros(len(lines_2D))
  tke_4D = np.zeros(len(lines_4D))
  tke_6D = np.zeros(len(lines_6D))

  uv_up = np.zeros(len(lines_up))
  uv_2D = np.zeros(len(lines_2D))
  uv_4D = np.zeros(len(lines_4D))
  uv_6D = np.zeros(len(lines_6D))

  urms_up = np.zeros(len(lines_up))
  urms_2D = np.zeros(len(lines_2D))
  urms_4D = np.zeros(len(lines_4D))
  urms_6D = np.zeros(len(lines_6D))

  vrms_up = np.zeros(len(lines_up))
  vrms_2D = np.zeros(len(lines_2D))
  vrms_4D = np.zeros(len(lines_4D))
  vrms_6D = np.zeros(len(lines_6D))

  for count in range(len(lines_up)):
    string_up = lines_up[count].split(' ')
    x_up[count] = np.float(string_up[0])
    y_up[count] = np.float(string_up[1])
    u_up[count] = np.float(string_up[2])
    v_up[count] = np.float(string_up[3])
    tke_up[count] = np.float(string_up[4])
    uv_up[count] = np.float(string_up[5])
    vrms_up[count] =np.float(string_up[6])
    urms_up[count] =np.float(string_up[7])
    
    string_2D = lines_2D[count].split(' ')
    x_2D[count] = np.float(string_2D[0])
    y_2D[count] = np.float(string_2D[1])
    u_2D[count] = np.float(string_2D[2])
    v_2D[count] = np.float(string_2D[3])
    tke_2D[count] = np.float(string_2D[4])
    uv_2D[count] = np.float(string_2D[5])
    vrms_2D[count] =np.float(string_2D[6])
    urms_2D[count] =np.float(string_2D[7])

    string_4D = lines_4D[count].split(' ')
    x_4D[count] = np.float(string_4D[0])
    y_4D[count] = np.float(string_4D[1])
    u_4D[count] = np.float(string_4D[2])
    v_4D[count] = np.float(string_4D[3])
    tke_4D[count] = np.float(string_4D[4])
    uv_4D[count] = np.float(string_4D[5])
    vrms_4D[count] =np.float(string_4D[6])
    urms_4D[count] =np.float(string_4D[7])

    string_6D = lines_6D[count].split(' ')
    x_6D[count] = np.float(string_6D[0])
    y_6D[count] = np.float(string_6D[1])
    u_6D[count] = np.float(string_6D[2])
    v_6D[count] = np.float(string_6D[3])
    tke_6D[count] = np.float(string_6D[4])
    uv_6D[count] = np.float(string_6D[5])
    vrms_6D[count] =np.float(string_6D[6])
    urms_6D[count] =np.float(string_6D[7])

X_up = x_up.reshape(nx,ny).transpose()
Y_up = y_up.reshape(nx,ny).transpose()
U_up = u_up.reshape(nx,ny).transpose()
V_up = v_up.reshape(nx,ny).transpose()
UV_up = uv_up.reshape(nx,ny).transpose()
TKE_up = tke_up.reshape(nx,ny).transpose()
Urms_up = urms_up.reshape(nx,ny).transpose()
Vrms_up = vrms_up.reshape(nx,ny).transpose()
kflux_up = UV_up*U_up


X_2D = x_2D.reshape(nx,ny).transpose()
Y_2D = y_2D.reshape(nx,ny).transpose()
U_2D = u_2D.reshape(nx,ny).transpose()
V_2D = v_2D.reshape(nx,ny).transpose()
UV_2D = uv_2D.reshape(nx,ny).transpose()
TKE_2D = tke_2D.reshape(nx,ny).transpose()
Urms_2D = urms_2D.reshape(nx,ny).transpose()
Vrms_2D = vrms_2D.reshape(nx,ny).transpose()
kflux_2D = UV_2D*U_2D

X_4D = x_4D.reshape(nx,ny).transpose()
Y_4D = y_4D.reshape(nx,ny).transpose()
U_4D = u_4D.reshape(nx,ny).transpose()
V_4D = v_4D.reshape(nx,ny).transpose()
UV_4D = uv_4D.reshape(nx,ny).transpose()
TKE_4D = tke_4D.reshape(nx,ny).transpose()
Urms_4D = urms_4D.reshape(nx,ny).transpose()
Vrms_4D = vrms_4D.reshape(nx,ny).transpose()
kflux_4D = UV_4D*U_4D

X_6D = x_6D.reshape(nx,ny).transpose()
Y_6D = y_6D.reshape(nx,ny).transpose()
U_6D = u_6D.reshape(nx,ny).transpose()
V_6D = v_6D.reshape(nx,ny).transpose()
UV_6D = uv_6D.reshape(nx,ny).transpose()
TKE_6D = tke_6D.reshape(nx,ny).transpose()
Urms_6D = urms_6D.reshape(nx,ny).transpose()
Vrms_6D = vrms_6D.reshape(nx,ny).transpose()
kflux_6D = UV_6D*U_6D

umean_up = U_up.mean(1)
vmean_up = V_up.mean(1)
uvmean_up = UV_up.mean(1)
tkemean_up = TKE_up.mean(1)
urmsmean_up = Urms_up.mean(1)
vrmsmean_up = Vrms_up.mean(1)
kfluxmean_up = kflux_up.mean(1)

umean_2D = U_2D.mean(1)
vmean_2D = V_2D.mean(1)
uvmean_2D = UV_2D.mean(1)
tkemean_2D = TKE_2D.mean(1)
urmsmean_2D = Urms_2D.mean(1)
vrmsmean_2D = Vrms_2D.mean(1)
kfluxmean_2D = kflux_2D.mean(1)

umean_4D = U_4D.mean(1)
vmean_4D = V_4D.mean(1)
uvmean_4D = UV_4D.mean(1)
tkemean_4D = TKE_4D.mean(1)
urmsmean_4D = Urms_4D.mean(1)
vrmsmean_4D = Vrms_4D.mean(1)
kfluxmean_4D = kflux_4D.mean(1)

umean_6D = U_6D.mean(1)
vmean_6D = V_6D.mean(1)
uvmean_6D = UV_6D.mean(1)
tkemean_6D = TKE_6D.mean(1)
urmsmean_6D = Urms_6D.mean(1)
vrmsmean_6D = Vrms_6D.mean(1)
kfluxmean_6D = kflux_6D.mean(1)

keyboard()

out = code_functions.writedata(Y_up.mean(1),umean_4D,plot_path+"umean_4D.dat")
out = code_functions.writedata(Y_up.mean(1),uvmean_4D,plot_path+"uv_4D.dat")
out = code_functions.writedata(Y_up.mean(1),vmean_4D,plot_path+"vmean_4D.dat")
out = code_functions.writedata(Y_up.mean(1),urmsmean_4D,plot_path+"urms_4D.dat")
out = code_functions.writedata(Y_up.mean(1),vrmsmean_4D,plot_path+"vrms_4D.dat")
out = code_functions.writedata(Y_up.mean(1),kfluxmean_4D,plot_path+"kflux_4D.dat")


out = code_functions.writedata(Y_up.mean(1),umean_2D,plot_path+"umean_2D.dat")
out = code_functions.writedata(Y_up.mean(1),uvmean_2D,plot_path+"uv_2D.dat")
out = code_functions.writedata(Y_up.mean(1),vmean_2D,plot_path+"vmean_2D.dat")
out = code_functions.writedata(Y_up.mean(1),urmsmean_2D,plot_path+"urms_2D.dat")
out = code_functions.writedata(Y_up.mean(1),vrmsmean_2D,plot_path+"vrms_2D.dat")
out = code_functions.writedata(Y_up.mean(1),kfluxmean_2D,plot_path+"kflux_2D.dat")





  
#  keyboard()


