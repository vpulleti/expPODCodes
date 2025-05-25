import sys
import numpy as np
import os
from pdb import set_trace as keyboard

numplane = 4
numrow = 8
numcol = 4

array = np.zeros((numplane,numrow,numcol))

for i in range(numplane):
  for j in range(numrow):
    for k in range(numcol):
      array[i,j,k] = (i+1)*(j+1)*(k+1)


p = np.fft.fft(array,axis=2)
print(array[:,0,0])
keyboard()
