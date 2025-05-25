import os
import sys
import math
import numpy as np
import scipy.sparse as scysparse
from pdb import set_trace
from time import sleep
import scipy.sparse as scysparse
import scipy.sparse.linalg as spysparselinalg  # sparse linear algebra
import scipy.linalg as scylinalg               # non-sparse linear algebra
import pylab as plt
import pickle
from pdb import set_trace
from matplotlib import rc as matplotlibrc
from matplotlib import cm
import time # has the equivalent of tic/toc
from numpy import linalg as LA
from pathlib import Path

###################################################################
######## function for calculating eigen(basis) function ##################
# U is velocity of interest (Y x X x time)
# vecs is eigenvectors (modes x modes)
# eigen functions (phi) (Y x X x modes)
##################################################################
def eigfunc(U,vecs):

  phi = np.matmul(U.transpose(),vecs)
  #sh = np.shape(U)
  #funcs = np.matmul(U.reshape((sh[0]*sh[1],sh[2])),vecs)
  #return funcs.reshape(sh)
  return phi
##################################################################
####### function for calculating time coefficient ################
# x is X (x meshgrid)
# y is Y (y meshgrid)
# U is streamwise velocity (Y x X x time)
# V is wallnormal velocity (Y x X x time)
# phi_u is basis function of U (Y x X x modes)
# phi_v is basis funciton of V (Y x X x modes)
# coe_mat is coefficient matrix (time x modes)
#################################################################
def timcoeff_old(x,y,U,V,phi_u,phi_v):
  sh = np.shape(U)
  area = 1E-06*(np.max(x)-np.min(x))*(np.max(y)-np.min(y))
  U_tvsxy = (U.reshape((sh[0]*sh[1],sh[2]))).transpose()
  V_tvsxy = (V.reshape((sh[0]*sh[1],sh[2]))).transpose()
  phiu_xyvsm = phi_u.reshape((sh[0]*sh[1],sh[2]))
  phiv_xyvsm = phi_v.reshape((sh[0]*sh[1],sh[2]))
  coe_mat = np.matmul(U_tvsxy,phiu_xyvsm)+np.matmul(V_tvsxy,phiv_xyvsm)
  #coe_mat = np.matmul((U.reshape((sh[0]*sh[1],sh[2]))).transpose(),phi_u.reshape((sh[0]*sh[1],sh[2])))+np.matmul((U.reshape((sh[0]*sh[1],sh[2]))).transpose(),phi_u.reshape((sh[0]*sh[1],sh[2])))
  coe_mat = coe_mat*area
  return coe_mat


def timcoeff(vel,phi):
  b = np.matmul(vel,phi)
  return b

##############################################################################
####### function for calculating reconstructed field based on no of modes ####
# 	coe is time x modes
#       phi is Y x X x modes
#       vel_reconst is Y x X x time
######################################################################

def reconst(coe,phi,st_mod,en_mod,Ny,Nx):

  #phi shape is points x modes
  #coe shape is time x modes 

  vel_reconst = np.zeros(np.shape(phi))

  num_tstp = np.shape(vel_reconst)[1]
  vel_reconst = vel_reconst.transpose() #time x points

  phi_T = np.transpose(phi)

  sh = np.shape(phi_T)

  for mode in range(st_mod,en_mod):
    for tstp in range(num_tstp):
        vel_reconst[tstp,:] += coe[tstp,mode]*phi_T[mode,:]

  vel = np.zeros((Ny,Nx,num_tstp))

  for step in range(num_tstp):
    vel[:,:,step] = np.reshape(vel_reconst[step,:],(Ny,Nx),order='C')

  return vel
######################################################################
######### function for calculating the norm for basis functions######
#    phi_u is basis function of u (Y x X x modes)
#    phi_v is basis function of v (Y x X x modes)
#    X is x meshgrid (Y x X)
#    Y is y meshgrid (Y x X)
######################################################################

def normlze(phi_vel):

  sh = np.shape(phi_vel)  #coded normc in equivalent such that sum of squares of the column (spatial modes) is 1

  for mode in range(sh[1]):
    phi_vel[:,mode] = phi_vel[:,mode]/np.sqrt(np.sum(phi_vel[:,mode]**2))

  return phi_vel







