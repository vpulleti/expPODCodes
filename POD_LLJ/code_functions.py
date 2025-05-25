import os
import sys
import math
import numpy as np
import scipy.sparse as scysparse
from pdb import set_trace as keyboard
from time import sleep
import scipy.sparse as scysparse
import scipy.sparse.linalg as spysparselinalg  # sparse linear algebra
import scipy.linalg as scylinalg               # non-sparse linear algebra
import pylab as plt
import pickle
from scipy.fftpack import fft
from matplotlib import rc as matplotlibrc
from matplotlib import cm
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path
from pdb import set_trace
from scipy.signal import get_window

def get_Fileinfo(folder,ext):
  snaplist = list()
  snapiternumlist = list()
  #set_trace()
  for itemname in os.listdir(folder):
    extension = itemname.split(".")[-1]
    if (extension==ext):
      filenumber = itemname.split(".")[0].split('Inlet')[-1]
      filename = itemname.split(".")[0]
      snapiternumlist.append(int(float(filenumber)))
      snaplist.append(folder+"/"+itemname)
  totsnaps = len(snaplist)
  return {'snapnames_ref':snaplist,'snap_num_list':snapiternumlist}



def bin_creation(ny,nx,tstp,prefx,ascii_path,bin_path):

  print("Doing for "+prefx)
  if Path(bin_path+"U").exists() and Path(bin_path+"V").exists() and Path(bin_path+"X").exists() and Path(bin_path+"Y").exists():

    print("The binary files exist and loading")

    with open(bin_path+"U",'rb') as file:
      U = pickle.load(file)

    with open(bin_path+"V",'rb') as file:
      V = pickle.load(file)
   
    with open(bin_path+"X",'rb') as file:
      X = pickle.load(file)

    with open(bin_path+"Y",'rb') as file:
      Y = pickle.load(file)

    print("The files are loaded and moving onto calculation")
  
  else:
    
    print("The binary files do not exist and creating them")
   
    X = np.zeros((ny,nx))
    Y = np.zeros((ny,nx))
    U = np.zeros((ny,nx,tstp))
    V = np.zeros((ny,nx,tstp))

    for i in range(tstp):
      if i<1000:
        name = prefx+"000"+str("%03d"%i)+".T000.D000.P000.H000.L.vec"
      else:
        name = prefx+"00"+str("%03d"%i)+".T000.D000.P000.H000.L.vec"   
      with open(ascii_path+name,'r') as file:
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

      X[:,:] = np.transpose(np.reshape(x,(nx,ny)))
      Y[:,:] = np.transpose(np.reshape(y,(nx,ny)))
      U[:,:,i] = np.transpose(np.reshape(u,(nx,ny)))
      V[:,:,i] = np.transpose(np.reshape(v,(nx,ny)))
      print("Done timestep %d"%i)

    with open(bin_path+"X",'wb') as file:
      pickle.dump(X,file)

    with open(bin_path+"Y",'wb') as file:
      pickle.dump(Y,file)

    with open(bin_path+"U",'wb') as file:
      pickle.dump(U,file)
      
    with open(bin_path+"V",'wb') as file:
      pickle.dump(V,file)

  print("The binary files are created and moving further")

  return {'u':U,'v':V,'x':X,'y':Y}    

def single_point_stats(ufluc,vfluc):

  ufluc2 = ufluc**2
  vfluc2 = vfluc**2
  uv2    = ufluc*vfluc
  sh = np.shape(ufluc)
  tstp = sh[2]
  urms = np.sqrt(ufluc2.sum(2)/tstp) 
  vrms = np.sqrt(vfluc2.sum(2)/tstp)

  uv = uv2.sum(2)/tstp

  return {'urms':urms,'vrms':vrms,'uv':uv}

def two_point_corrl(loc_u,loc_v,urms,vrms,u,v):

  sh = np.shape(u)
  tstp = sh[2]
  rho_uu = np.zeros(sh)
  rho_vv = np.zeros(sh)
  ref_u = u[loc_u[0],loc_u[1],:]**2
  ref_v = v[loc_v[0],loc_u[1],:]**2
  #keyboard() 
  
  ref_u = np.sqrt(ref_u.mean(0))
  ref_v = np.sqrt(ref_v.mean(0))
   
  
  for i in range(tstp):
    rho_uu[:,:,i] = (u[loc_u[0],loc_u[1],i]*u[:,:,i])
    rho_vv[:,:,i] = (v[loc_v[0],loc_v[1],i]*v[:,:,i])
    
  rho_uu = rho_uu.sum(2)/tstp
  rho_vv = rho_vv.sum(2)/tstp

  fac_uu = 1./(ref_u*urms)
  fac_vv = 1./(ref_v*vrms)

  rho_uu = fac_uu*rho_uu
  rho_vv = fac_vv*rho_vv
  #keyboard() 
  return{'rhouu':rho_uu,'rhovv':rho_vv}

def tke(u,v):
  tstp = np.shape(u)[2]
  tke = 0.5*(u**2+v**2)
  tke = tke.sum(2)/tstp
  return {'tke':tke}
  
# calculating the spectra
def spectra_1D(fld,x):

  sh = np.shape(fld)
  fld_pad = np.zeros((sh[0],128,sh[2]))
  fld_pad[:,:sh[1],:] = fld
  
 
  # windowing to reduce spurious erros becasue of non-periodicity
  win = get_window('hanning',np.shape(fld_pad)[1])
  sca_fact = 1./np.sqrt(np.mean(win**2))
  
  fld_wind = np.zeros(np.shape(fld_pad))
 
  for i in range(np.shape(fld_pad)[1]):
    fld_wind[:,i,:] = win[i]*fld_pad[:,i,:]

  fld_wind_mean = fld_wind.mean(1) 
 
  fld_windod = np.zeros(np.shape(fld_wind))

  for xco in range(np.shape(fld_pad)[1]):
    for time in range(np.shape(fld_pad)[2]):
      fld_windod[:,xco,time] = sca_fact*(fld_wind[:,xco,time]-fld_wind_mean[xco,time])
 
  
  #keyboard()   
  fld_wind_hat=fft(fld_windod,axis=1)
  #keyboard()
  length = np.max(x)-np.min(x)
  wave_number = (np.pi/length)*np.linspace(1,np.shape(fld_pad)[1],np.shape(fld_pad)[1])
  #keyboard()
  #Euu = abs(fld_wind_hat*np.conj(fld_wind_hat)*np.shape(fld_pad)[1]).mean(2)
  Euu = abs(fld_wind_hat*np.conj(fld_wind_hat))
  #keyboard()
  Euu_mean = Euu.mean(2)
  return {'k':wave_number,'spectra':Euu_mean}

  
def normalize(case):
  if case=="array_7D_L" or case=="Single_2D_L" or case=='single_7D_L':
    D = 120
    vel = 6.510197676039683
    delta = 133.59584
  elif case=="array_7D_M" or case=="Single_2D_M" or case=="single_7D_M" or case=='Inlet':
    D = 120
    vel = 9.225284880269841
    delta = 193.59584
  elif case=="array_7D_H" or case=="single_2D_H" or case=="Single_7D_H":
    D = 120
    vel = 6.301402394301588
    delta = 253.59584
  return {'D':D,'vel':vel,'delta':delta}


def contrwritedata(var1,var2,data,name):
  array = (np.array([var1.flatten(order='C'),var2.flatten(order='C'),data.flatten(order='C')])).transpose()  
  np.savetxt(name,array,fmt="%2.6f")
  return 0


def writedata(var1,var2,name):
  array = (np.array([var1,var2])).transpose()
  np.savetxt(name,array,fmt="%2.6f")
  return 0


def firstderv(var1,var2,direct):
  sh = np.shape(var2)
  derv = np.zeros(sh)
  if direct==0:
    derv[0,:] = (var2[1,:]-var1[0,:])/(var1[1,:]-var1[0,:])
    for i in range(1,sh[0]-1):
      derv[i,:] = (var2[i+1,:]-var2[i-1,:])/(var1[i+1,:]-var1[i-1,:])
    derv[sh[0]-1,:] = (var2[sh[0]-1,:]-var2[sh[0]-2,:])/(var1[sh[0]-1,:]-var1[sh[0]-2,:])
  elif direct==1:
    derv[:,0] = (var2[:,1]-var1[:,0])/(var1[:,1]-var1[:,0]) 
    for i in range(1,sh[1]-1):
      derv[:,i] = (var2[:,i+1]-var2[:,i-1])/(var1[:,i+1]-var1[:,i-1])
    derv[:,sh[1]-1] =  (var2[:,sh[1]-1]-var2[:,sh[1]-2])/(var1[:,sh[1]-1]-var1[:,sh[1]-2])
  return derv

def tecplot_writer(filename, variables, X=[], Y=[], Z=[]):
    """
    X, Y, Z are the lists of xyz coordinates. If not provided, intergers
    from 0 will be used.
    `variables` is a dict of the variables to store with the variable names as
    the keys. Each variable should be 2 or 3 dimensional array using numpy's
    row-major order.
    Check the test function to see how to create input data structure.
    Notice that tecplot format use 'column-major order' as in Fortran, which is
    different from that of Numpy or C.
    """
    if filename[-4:] != '.dat':
        filename += '.dat'

    for var in variables:
      variables[var] = np.nan_to_num(variables[var])
    with open(filename, 'w') as f:
      if len(Z)==0:
        if len(Y)==0:
          ## 1D case
          f.write('Variables="X"')
          for key in variables.keys():
            f.write(',"' + key + '"')
          f.write('\n\nZone I='+str(np.size(X))+'.F=POINT\n')
          for i in range(np.shape(X)[0]):
            f.write(str(X[i]))
            for var in variables.values():
              f.write(' ' + str(var[i]))
            f.write('\n')
        elif len(X)==0:
          ## 1D case
          f.write('Variables="Y"')
          for key in variables.keys():
            f.write(',"' + key + '"')
          f.write('\n\nZone I='+str(np.size(Y))+'.F=POINT\n')
          for i in range(np.shape(Y)[0]):
            f.write(str(Y[i]))
            for var in variables.values():
              f.write(' ' + str(var[i]))
          f.write('\n')
        else:
          f.write('X, Y')
          for key in variables.keys():
            f.write(',"' + key +'"')
          f.write('\n')
          #f.write('\n\nZone I='+str(np.size(X))+', J='+str(np.size(Y))+', F=POINT\n')
          for j in range(np.shape(Y)[1]):
            for i in range(np.shape(X)[0]):
              f.write(str(X[i][j]) + ' ' +str(Y[i][j]))
              for var in variables.values():
                f.write(' ' + str(var[i][j]))
              f.write('\n')
      else:
        f.write('Variables = "X", "Y", "Z"')
        for key in variables.keys():
          f.write(', "' + key + '"')
        f.write('\n\nZone I=' + str(np.size(X)) + ', J=' + str(np.size(Y)) +', K=' + str(np.size(Z)) + ', F=POINT\n')
        for k in range(np.shape(Z)[2]):
          for j in range(np.shape(Y)[0]):
            for i in range(np.shape(X)[1]):
              f.write(str(X[j,i,k]) + ' ' + str(Y[j,i,k]) + ' ' + str(Z[j,i,k]))
              for var in variables.values():
                f.write(' ' + str(var[j,i,k]))
              f.write('\n')
      print("Saved :" +filename)
'''       
      ## 2D case
      if len(Z) == 0:
          f.write('Variables = "X", "Y"')
          for key in variables.keys():
            f.write(', "' + key + '"')
          f.write('\n\nZone I='+str(np.size(X))+', J='+str(np.size(Y))+', F=POINT\n')

          for j in range(np.shape(Y)[1]):
            for i in range(np.shape(X)[0]):
              f.write(str(X[i][j]) + ' ' + str(Y[i][j]))
              for var in variables.values():
                f.write(' ' + str(var[i][j]))
              f.write('\n')
'''
        ## 3D case



