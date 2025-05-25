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
import pandas as pd
from scipy import signal
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path
import matplotlib.pyplot as plt
import code_functions
from pdb import set_trace
import matplotlib.colors as mcolors

mpl.rc('text', usetex = True)
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

colormap = pd.DataFrame({"name" : ["cus_red","cus_blue","cus_grey"],"color":["#cb0000","#000071","#959595"]})

colormap["name"] = colormap["name"].apply(lambda x: x.lower())   
c = dict(zip(*colormap.values.T))
mcolors.get_named_colors_mapping().update(c)

DPI = 400


plot_folder = '/scratch/bell/vpulleti/winglet_data/binary/plots/'
'''
files_path = ['/scratch/bell/vpulleti/winglet_data/binary/BL_w0_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w0_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w1_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w1_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w2_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w2_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w3_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/BL_w3_7D/']
'''
files_path = ['/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w0_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w0_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w1_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w1_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w2_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w2_7D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w3_2D/',
         '/scratch/bell/vpulleti/winglet_data/binary/LLJ_POS_w3_7D/']


'''
with open(files_path[0]+'BL_w0_2Ddata.pickle','rb') as file:
  data_w02D = pickle.load(file)

with open(files_path[1]+'BL_w0_7Ddata.pickle','rb') as file:
  data_w07D = pickle.load(file)

with open(files_path[2]+'BL_w1_2Ddata.pickle','rb') as file:
  data_w12D = pickle.load(file)

with open(files_path[3]+'BL_w1_7Ddata.pickle','rb') as file:
  data_w17D = pickle.load(file)

with open(files_path[4]+'BL_w2_2Ddata.pickle','rb') as file:
  data_w22D = pickle.load(file)

with open(files_path[5]+'BL_w2_7Ddata.pickle','rb') as file:
  data_w27D = pickle.load(file)

with open(files_path[6]+'BL_w3_2Ddata.pickle','rb') as file:
  data_w32D = pickle.load(file)

with open(files_path[7]+'BL_w3_7Ddata.pickle','rb') as file:
  data_w37D = pickle.load(file)

'''
with open(files_path[0]+'LLJ_POS_w0_2Ddata.pickle','rb') as file:
  data_w02D = pickle.load(file)

with open(files_path[1]+'LLJ_POS_w0_7Ddata.pickle','rb') as file:
  data_w07D = pickle.load(file)

with open(files_path[2]+'LLJ_POS_w1_2Ddata.pickle','rb') as file:
  data_w12D = pickle.load(file)

with open(files_path[3]+'LLJ_POS_w1_7Ddata.pickle','rb') as file:
  data_w17D = pickle.load(file)

with open(files_path[4]+'LLJ_POS_w2_2Ddata.pickle','rb') as file:
  data_w22D = pickle.load(file)

#with open(files_path[5]+'LLJ_POS_w2_7Ddata.pickle','rb') as file:
#  data_w27D = pickle.load(file)

with open(files_path[6]+'LLJ_POS_w3_2Ddata.pickle','rb') as file:
  data_w32D = pickle.load(file)

with open(files_path[7]+'LLJ_POS_w3_7Ddata.pickle','rb') as file:
  data_w37D = pickle.load(file)


avgChord_len = 0.1687 #Leo's 2015 paper on wingtip vortices


y_norm = data_w02D['Y_norm'][10,:]
Y_norm = data_w02D['Y_norm']
X_norm_2D = data_w02D['X_norm']
X_norm_7D = data_w07D['X_norm']


X_2D,Y_2D = np.meshgrid(data_w02D['X_norm'][:,10],y_norm)


num_snaps = np.shape(data_w02D['lambda_real'])[2]

loc1T_dot5D = np.abs(X_norm_2D[:,10]-0.5).argmin()
loc1T_1D = np.abs(X_norm_2D[:,10]-1.0).argmin()
loc1T_1dot5D = np.abs(X_norm_2D[:,10]-1.5).argmin()


loc1T_2D = np.abs(X_norm_2D[:,10]-2.0).argmin()
loc1T_2dot5D = np.abs(X_norm_2D[:,10]-2.5).argmin()


loc1T_6D = np.abs(X_norm_7D[:,10]-6.0).argmin()
loc1T_6dot5D = np.abs(X_norm_7D[:,10]-6.5).argmin()


loc1T_7D = np.abs(X_norm_7D[:,10]-7.0).argmin()
loc1T_7dot5D = np.abs(X_norm_7D[:,10]-7.5).argmin()

code_functions.tecplot_writer(plot_folder+'lambda2_w17D_510.dat',{'lambda2_LSM':data_w17D['lambda_LSM'][:,:,510],'lambda2_real':data_w17D['lambda_real'][:,:,510]},X_norm_7D,data_w07D['Y_norm'])

plot_w01D = {}
plot_w01D['uv_realMean'] = signal.savgol_filter(np.mean(data_w02D['uv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['uv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w02D['uv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w01D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w02D = {}
plot_w02D['uv_realMean'] = signal.savgol_filter(np.mean(data_w02D['uv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w02D['Uuv_realMean'] = np.mean(data_w02D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0)
plot_w02D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w02D['Uuv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['uv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['uv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w02D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w02D['TKE_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w06D = {}
plot_w06D['uv_realMean'] = signal.savgol_filter(np.mean(data_w07D['uv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['uv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['uv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w06D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w07D = {}
plot_w07D['uv_realMean'] = signal.savgol_filter(np.mean(data_w07D['uv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w07D['Uuv_realMean'] = np.mean(data_w07D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0)
plot_w07D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['Uuv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['uv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['uv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w07D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w07D['TKE_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

# w1

plot_w11D = {}
plot_w11D['uv_realMean'] = signal.savgol_filter(np.mean(data_w12D['uv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['uv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['uv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w11D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')


plot_w12D = {}
plot_w12D['uv_realMean'] = signal.savgol_filter(np.mean(data_w12D['uv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w12D['Uuv_realMean'] = np.mean(data_w12D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0)
plot_w12D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['uv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['Uuv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['uv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w12D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w12D['TKE_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w16D = {}
plot_w16D['uv_realMean'] = signal.savgol_filter(np.mean(data_w17D['uv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['uv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['uv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w16D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w17D = {}
plot_w17D['uv_realMean'] = signal.savgol_filter(np.mean(data_w17D['uv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w17D['Uuv_realMean'] = np.mean(data_w17D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0)
plot_w17D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['uv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['Uuv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['uv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w17D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w17D['TKE_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

# w2

plot_w21D = {}
plot_w21D['uv_realMean'] = signal.savgol_filter(np.mean(data_w22D['uv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['uv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['uv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w21D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w22D = {}
plot_w22D['uv_realMean'] = signal.savgol_filter(np.mean(data_w22D['uv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w22D['Uuv_realMean'] = np.mean(data_w22D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0)
plot_w22D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['uv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['Uuv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['uv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w22D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w22D['TKE_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

'''
plot_w26D = {}
plot_w26D['uv_realMean'] = signal.savgol_filter(np.mean(data_w27D['uv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['uv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['uv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w26D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w27D = {}
plot_w27D['uv_realMean'] = signal.savgol_filter(np.mean(data_w27D['uv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w27D['Uuv_realMean'] = np.mean(data_w27D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0)
plot_w27D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['uv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['Uuv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['uv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w27D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w27D['TKE_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
'''

# w3

plot_w31D = {}
plot_w31D['uv_realMean'] = signal.savgol_filter(np.mean(data_w32D['uv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['uv_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['uv_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_real'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_LSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w31D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_SSM'][loc1T_dot5D:loc1T_1dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w32D = {}
plot_w32D['uv_realMean'] = signal.savgol_filter(np.mean(data_w32D['uv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w32D['Uuv_realMean'] = np.mean(data_w32D['Uuv_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0)
plot_w32D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['uv_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['Uuv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['uv_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_real'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_LSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w32D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w32D['TKE_SSM'][loc1T_1dot5D:loc1T_2dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w36D = {}
plot_w36D['uv_realMean'] = signal.savgol_filter(np.mean(data_w37D['uv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['uv_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['uv_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_real'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_LSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w36D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_SSM'][:loc1T_6dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')

plot_w37D = {}
plot_w37D['uv_realMean'] = signal.savgol_filter(np.mean(data_w37D['uv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['Uuv_realMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
#plot_w37D['Uuv_realMean'] = np.mean(data_w37D['Uuv_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0)
plot_w37D['Uuv_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['uv_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['uv_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['Uuv_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['Uuv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['uv_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['uv_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['TKE_realMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_real'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['TKE_LSMMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_LSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')
plot_w37D['TKE_SSMMean'] = signal.savgol_filter(np.mean(data_w37D['TKE_SSM'][loc1T_6dot5D:loc1T_7dot5D,:],axis=0),window_length=15,polyorder=2,mode='nearest')


figwidth  = 10
figheight = 9
lineWidth = 2
textfontSize = 32
gcafontSize = 28

mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize

indx_need= np.linspace(2,2+20*2,2).astype(int)
phas_data=data_w02D['lambda_real'][:,:,indx_need]

phas_lamb=np.mean(phas_data,axis=2)
fig = plt.figure(0,figsize=(figwidth,figheight))
ax=plt.axes()
'''
for snap in range(10):
  plt.streamplot(X_2D,Y_2D,np.transpose(data_w12D['uinst_real'][:,:,snap]),np.transpose(data_w12D['vinst_real'][:,:,snap]),density=2)
  plt.xlim([0.5,2.6])
  plt.ylim([-1.0,1.5])
  plt.xlabel(r'$x/D$',fontsize=textfontSize)
  plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
  plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
  plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
  plt.tight_layout()
  plt.show()
'''
cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['lambda_real'][:,:,8],np.linspace(0,0.007,100),cmap='Reds')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
 c.set_edgecolor("face")
plt.streamplot(X_2D,Y_2D,np.transpose(data_w12D['uinst_real'][:,:,8]),np.transpose(data_w12D['vinst_real'][:,:,8]),density=2)
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.show()

set_trace()
'''
for snap in range(3,27,5):
  cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['lambda_real'][:,:,snap],np.linspace(0,0.01,100),cmap='jet')
  cb = plt.colorbar()
  cb.ax.tick_params(labelsize=gcafontSize)
  for c in cnt.collections:
   c.set_edgecolor("face")
  plt.xlim([0.5,2.6])
  plt.ylim([-1.0,1.5])
  plt.xlabel(r'$x/D$',fontsize=textfontSize)
  plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
  plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
  plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
  plt.tight_layout()
  plt.show()
set_trace()
'''


Names  = {'FileName':plot_folder+'TKE_LLJ_w0_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w02D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w02D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w02D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'TKE_LLJ_w0_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w07D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w07D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w07D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

#w1

Names  = {'FileName':plot_folder+'TKE_LLJ_w1_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w12D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w12D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w12D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'TKE_LLJ_w1_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w17D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w17D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w17D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

#w2

Names  = {'FileName':plot_folder+'TKE_LLJ_w2_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w22D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w22D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w22D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'TKE_w2_2D.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

'''
Names  = {'FileName':plot_folder+'TKE_LLJ_w2_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w27D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w27D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w27D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(plot_folder+'TKE_w2_7D.png',dpi=DPI)
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()
'''

#w3
Names  = {'FileName':plot_folder+'TKE_LLJ_w3_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w32D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w32D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w32D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'TKE_LLJ_w3_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w37D['TKE_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w37D['TKE_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w37D['TKE_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([0,0.025])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$\textrm{TKE}/U^2_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


#line EF
Names  = {'FileName':plot_folder+'line_EF_LLJ_w0_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w02D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w02D['Uuv_LSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.plot(plot_w02D['Uuv_SSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'line_EF_LLJ_w0_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w07D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w07D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w07D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

#w1

Names  = {'FileName':plot_folder+'line_EF_LLJ_w1_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w12D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w12D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w12D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'line_EF_LLJ_w1_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w17D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w17D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w17D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

#w2

Names  = {'FileName':plot_folder+'line_EF_LLJ_w2_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w22D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w22D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w22D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


'''
Names  = {'FileName':plot_folder+'line_EF_LLJ_w2_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w27D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w27D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w27D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()
'''

#w3
Names  = {'FileName':plot_folder+'line_EF_LLJ_w3_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w32D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w32D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w32D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'line_EF_LLJ_w3_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w37D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=3)
plt.plot(plot_w37D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w37D['Uuv_SSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{All modes}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'line_EF_LLJ_LSMcomparison_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w02D['Uuv_SSMMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w12D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w22D['Uuv_LSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.plot(plot_w32D['Uuv_LSMMean'],y_norm,color='green',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
ax.legend([r'$\textrm{w}_0$',r'$\textrm{w}_1$',r'$\textrm{w}_2$',r'$\textrm{w}_3$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'line_EF_LSMcomparison_2D.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'line_EF_LLJ_LSMcomparison_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w07D['Uuv_LSMMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w17D['Uuv_LSMMean'],y_norm,color='cus_red',linewidth=lineWidth)
#plt.plot(plot_w27D['Uuv_LSMMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.plot(plot_w37D['Uuv_LSMMean'],y_norm,color='green',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
#ax.legend([r'$\textrm{w}_0$',r'$\textrm{w}_1$',r'$\textrm{w}_2$',r'$\textrm{w}_3$'],fontsize=gcafontSize)
ax.legend([r'$\textrm{w}_0$',r'$\textrm{w}_1$',r'$\textrm{w}_3$'],fontsize=gcafontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'line_EF_LSMcomparison_7D.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'line_EF_LLJ_comparison_2D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w02D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w12D['Uuv_realMean'],y_norm,color='cus_red',linewidth=lineWidth)
plt.plot(plot_w22D['Uuv_realMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.plot(plot_w32D['Uuv_realMean'],y_norm,color='green',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.show()
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'line_EF_LLJ_comparison_7D.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
plt.plot(plot_w07D['Uuv_realMean'],y_norm,'ok',markevery=2,markersize=5)
plt.plot(plot_w17D['Uuv_realMean'],y_norm,color='cus_red',linewidth=lineWidth)
#plt.plot(plot_w27D['Uuv_realMean'],y_norm,color='cus_blue',linewidth=lineWidth)
plt.plot(plot_w37D['Uuv_realMean'],y_norm,color='green',linewidth=lineWidth)
plt.axhline(y=0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.axhline(y=0.0,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.8)
plt.axhline(y=-0.5,linestyle='dashed',color='cus_grey',linewidth=lineWidth-1.5)
plt.xlim([-0.0075,0.0075])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_\textrm{hub}$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.show()
print("Saving file name :" + Names['FileName'])
plt.close()
#uv

min_val = np.min(np.min(-data_w02D['uv_real']))
max_val = np.max(np.max(-data_w02D['uv_real']))

figwidth  = 14
figheight = 10
lineWidth = 1
textfontSize = 32
gcafontSize = 28

Names  = {'FileName':plot_folder+'w02D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['uv_LSM'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w02D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['uv_real'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'w07D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w07D['uv_LSM'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w07D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w07D['uv_real'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w12D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w12D['uv_LSM'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w12D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w12D['uv_real'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w17D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w17D['uv_LSM'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w17D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w17D['uv_real'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w22D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w22D['uv_LSM'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w22D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w12D['uv_real'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

'''
Names  = {'FileName':plot_folder+'w27D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w27D['uv_LSM'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w27D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w27D['uv_real'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()
'''

#w3
Names  = {'FileName':plot_folder+'w32D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w32D['uv_LSM'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w32D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w32D['uv_real'],np.linspace(-0.015,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w37D_RS_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w37D['uv_LSM'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w37D_RS_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w37D['uv_real'],np.linspace(-0.002,0.006,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


#Uuv

Names  = {'FileName':plot_folder+'w02D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['Uuv_LSM'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w02D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_2D,Y_norm,data_w02D['Uuv_real'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


Names  = {'FileName':plot_folder+'w07D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w07D['Uuv_LSM'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w07D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w07D['Uuv_real'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w12D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w12D['Uuv_LSM'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w12D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w12D['Uuv_real'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w17D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w17D['Uuv_LSM'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w17D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w17D['Uuv_real'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w22D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w22D['Uuv_LSM'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'w22D_EF_LSM.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w22D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w22D['Uuv_real'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'w22D_EF_real.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

'''
Names  = {'FileName':plot_folder+'w27D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w27D['Uuv_LSM'],np.linspace(-0.002,0.0065,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'w27D_EF_LSM.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w27D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w27D['Uuv_real'],np.linspace(-0.002,0.0065,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
plt.savefig(plot_folder+'w27D_EF_real.png',dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()
'''

#w3
Names  = {'FileName':plot_folder+'w32D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w32D['Uuv_LSM'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w32D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_2D,Y_norm,data_w32D['Uuv_real'],np.linspace(-0.01,0.01,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([0.5,2.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w37D_EF_LLJ_LSM.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=cnt=plt.contourf(X_norm_7D,Y_norm,data_w37D['Uuv_LSM'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

Names  = {'FileName':plot_folder+'w37D_EF_LLJ_real.pdf'}
fig = plt.figure(0,figsize=(figwidth,figheight))
ax= plt.axes()
cnt=plt.contourf(X_norm_7D,Y_norm,data_w37D['Uuv_real'],np.linspace(-0.002,0.005,100),cmap='jet')
cb = plt.colorbar()
cb.ax.tick_params(labelsize=gcafontSize)
for c in cnt.collections:
    c.set_edgecolor("face")
plt.xlim([5.5,7.6])
plt.ylim([-1.0,1.5])
plt.xlabel(r'$x/D$',fontsize=textfontSize)
plt.ylabel(r'$(y-h)/D$',fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()



