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
from pdb import set_trace

mpl.rc('text', usetex = True)
mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

DPI = 400


plot_folder = '/scratch/bell/vpulleti/LLJ_PIV/binary/'

files_path = ['/scratch/bell/vpulleti/LLJ_PIV/binary/array_L/array_L_',
         '/scratch/bell/vpulleti/LLJ_PIV/binary/array_H/array_H_',
         '/scratch/bell/vpulleti/LLJ_PIV/binary/single_2D_L/single_2D_L_',
         '/scratch/bell/vpulleti/LLJ_PIV/binary/single_2D_H/single_2D_H_',
         '/scratch/bell/vpulleti/LLJ_PIV/binary/single_7D_L/single_7D_L_',
         '/scratch/bell/vpulleti/LLJ_PIV/binary/single_7D_H/single_7D_H_']

with open(files_path[0]+'array_7D_Ldata.pickle','rb') as file:
  data_arr7DL = pickle.load(file)

with open(files_path[1]+'array_7D_Hdata.pickle','rb') as file:
  data_arr7DH = pickle.load(file)

with open(files_path[2]+'Single_2D_Ldata.pickle','rb') as file:
  data_sing2DL = pickle.load(file)

with open(files_path[3]+'single_2D_Hdata.pickle','rb') as file:
  data_sing2DH = pickle.load(file)

with open(files_path[4]+'single_7D_Ldata.pickle','rb') as file:
  data_sing7DL = pickle.load(file)

with open(files_path[5]+'Single_7D_Hdata.pickle','rb') as file:
  data_sing7DH = pickle.load(file)


data_arr7DL['Uuv_realMean'] = np.mean(data_arr7DL['Uuv_realMean'],axis=1)
data_arr7DL['Uuv_210Mean'] = np.mean(data_arr7DL['Uuv_210Mean'],axis=1)
data_arr7DL['Uuvu_ssmMean'] = np.mean(data_arr7DL['Uuvu_ssmMean'],axis=1)
arr7DL_y = data_arr7DL['Y_norm'][:,10]

data_arr7DH['Uuv_realMean'] = np.mean(data_arr7DH['Uuv_realMean'],axis=1)
data_arr7DH['Uuv_210Mean'] = np.mean(data_arr7DH['Uuv_210Mean'],axis=1)
data_arr7DH['Uuvu_ssmMean'] = np.mean(data_arr7DH['Uuvu_ssmMean'],axis=1)

arr7DH_y = data_arr7DH['Y_norm'][:,10]

data_sing2DL['Uuv_realMean'] = np.mean(data_sing2DL['Uuv_realMean'],axis=1)
data_sing2DL['Uuv_210Mean'] = np.mean(data_sing2DL['Uuv_210Mean'],axis=1)
data_sing2DL['Uuvu_ssmMean'] = np.mean(data_sing2DL['Uuvu_ssmMean'],axis=1)

sing2DL_y = data_sing2DL['Y_norm'][:,10]

data_sing2DH['Uuv_realMean'] = np.mean(data_sing2DH['Uuv_realMean'],axis=1)
data_sing2DH['Uuv_210Mean'] = np.mean(data_sing2DH['Uuv_210Mean'],axis=1)
data_sing2DH['Uuvu_ssmMean'] = np.mean(data_sing2DH['Uuvu_ssmMean'],axis=1)

sing2DH_y = data_sing2DH['Y_norm'][:,10]

data_sing7DL['Uuv_realMean'] = np.mean(data_sing7DL['Uuv_realMean'],axis=1)
data_sing7DL['Uuv_210Mean'] = np.mean(data_sing7DL['Uuv_210Mean'],axis=1)
data_sing7DL['Uuv_ssmMean'] = np.mean(data_sing7DL['Uuvu_ssmMean'],axis=1)

sing7DL_y = data_sing7DL['Y_norm'][:,10]

data_sing7DH['Uuv_realMean'] = np.mean(data_sing7DH['Uuv_realMean'],axis=1)
data_sing7DH['Uuv_210Mean'] = np.mean(data_sing7DH['Uuv_210Mean'],axis=1)
data_sing7DH['Uuvu_ssmMean'] = np.mean(data_sing7DH['Uuvu_ssmMean'],axis=1)

sing7DH_y = data_sing7DH['Y_norm'][:,10]

figwidth  = 12
figheight = 10
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'positive_shear_array.pdf'}

plt.plot(data_arr7DL['Uuv_realMean'],arr7DL_y,'k',linewidth=lineWidth)
plt.plot(data_arr7DL['Uuv_210Mean'],arr7DL_y,'r',linewidth=lineWidth)
plt.plot(data_arr7DL['Uuvu_ssmMean'],arr7DL_y,'b',linewidth=lineWidth)

plt.grid(True,which='both',linewidth=1)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
plt.ylim([-0.5,2.0])
plt.xlim([-0.01,0.015])
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'negative_shear_array.pdf'}

plt.plot(data_arr7DH['Uuv_realMean'],arr7DH_y,'k',linewidth=lineWidth)
plt.plot(data_arr7DH['Uuv_210Mean'],arr7DH_y,'r',linewidth=lineWidth)
plt.plot(data_arr7DH['Uuvu_ssmMean'],arr7DH_y,'b',linewidth=lineWidth)
plt.ylim([-1.5,1.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'positive_shear_single7D.pdf'}

plt.plot(data_sing7DL['Uuv_realMean'],sing7DL_y,'k',linewidth=lineWidth)
plt.plot(data_sing7DL['Uuv_210Mean'],sing7DL_y,'r',linewidth=lineWidth)
plt.plot(data_sing7DL['Uuvu_ssmMean'],sing7DL_y,'b',linewidth=lineWidth)
plt.ylim([-0.5,2.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()

figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'negative_shear_single7D.pdf'}

plt.plot(data_sing7DH['Uuv_realMean'],sing7DH_y,'k',linewidth=lineWidth)
plt.plot(data_sing7DH['Uuv_210Mean'],sing7DH_y,'r',linewidth=lineWidth)
plt.plot(data_sing7DH['Uuvu_ssmMean'],sing7DH_y,'b',linewidth=lineWidth)
plt.ylim([-1.5,1.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'positive_shear_single2D.pdf'}

plt.plot(data_sing2DL['Uuv_realMean'],sing2DL_y,'k',linewidth=lineWidth)
plt.plot(data_sing2DL['Uuv_210Mean'],sing2DL_y,'r',linewidth=lineWidth)
plt.plot(data_sing2DL['Uuvu_ssmMean'],sing2DL_y,'b',linewidth=lineWidth)
plt.ylim([-0.5,2.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'negative_shear_single2D.pdf'}

plt.plot(data_sing2DH['Uuv_realMean'],sing2DH_y,'k',linewidth=lineWidth)
plt.plot(data_sing2DH['Uuv_210Mean'],sing2DH_y,'r',linewidth=lineWidth)
plt.plot(data_sing2DH['Uuvu_ssmMean'],sing2DH_y,'b',linewidth=lineWidth)
plt.ylim([-1.5,1.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.legend([r'$\textrm{Multi-scale}$',r'$\textrm{LSM}$',r'$\textrm{SSM}$'],fontsize=textfontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()


figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'positive_shear_arr_and_single7D.pdf'}

plt.plot(data_sing7DL['Uuv_realMean'],sing7DL_y,'--k',linewidth=lineWidth)
plt.plot(data_sing7DL['Uuv_210Mean'],sing7DL_y,'--r',linewidth=lineWidth)
plt.plot(data_sing7DL['Uuvu_ssmMean'],sing7DL_y,'--b',linewidth=lineWidth)

plt.plot(data_arr7DL['Uuv_realMean'],arr7DL_y,'k',linewidth=lineWidth)
plt.plot(data_arr7DL['Uuv_210Mean'],arr7DL_y,'r',linewidth=lineWidth)
plt.plot(data_arr7DL['Uuvu_ssmMean'],arr7DL_y,'b',linewidth=lineWidth)
plt.ylim([-0.5,2.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()



figwidth = 12
figheight = 11
lineWidth = 3
textfontSize = 28
gcafontSize = 24


mpl.rcParams['xtick.labelsize'] = gcafontSize
mpl.rcParams['ytick.labelsize'] = gcafontSize


fig = plt.figure(0,figsize=(figwidth,figheight))
ax = plt.axes()
Names  = {'FileName':plot_folder+'negative_shear_arr_and_single7D.pdf'}

plt.plot(data_sing7DH['Uuv_realMean'],sing7DH_y,'--k',linewidth=lineWidth)
plt.plot(data_sing7DH['Uuv_210Mean'],sing7DH_y,'--r',linewidth=lineWidth)
plt.plot(data_sing7DH['Uuvu_ssmMean'],sing7DH_y,'--b',linewidth=lineWidth)

plt.plot(data_arr7DH['Uuv_realMean'],arr7DH_y,'k',linewidth=lineWidth)
plt.plot(data_arr7DH['Uuv_210Mean'],arr7DH_y,'r',linewidth=lineWidth)
plt.plot(data_arr7DH['Uuvu_ssmMean'],arr7DH_y,'b',linewidth=lineWidth)
plt.ylim([-1.5,1.0])
plt.xlim([-0.01,0.015])
plt.grid(True,which='both',linewidth=1)
plt.setp(ax.get_xticklabels(),fontsize=gcafontSize)
plt.setp(ax.get_yticklabels(),fontsize=gcafontSize)
ax.set_xlabel(r'$-\overline{U}\overline{u\prime v\prime}/U^3_{hub}$',fontsize=textfontSize)
ax.set_ylabel(r'$(y-H)/D$',fontsize=textfontSize)
plt.tight_layout()
plt.savefig(Names['FileName'],dpi=DPI)
print("Saving file name :" + Names['FileName'])
plt.close()
