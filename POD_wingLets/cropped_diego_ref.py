import time
start=time.time()
import numpy as np
import pickle
from pathlib import Path
import glob
import matplotlib.axis
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sc
from scipy import interpolate
from scipy import signal
from sklearn.neighbors import LocalOutlierFactor
from scipy import stats
from scipy.interpolate import griddata


import matplotlib.ticker

class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=-1, vmax=1):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)

plt.rcParams['font.family']='Dejavu serif'
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rc('axes',labelsize=14)
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)


D=120 # Rotor diameter
y_h=140 #BL and LLJ_POS hub height
y_H=200 #LLJ_HUB hub height

zeros = np.zeros((143,110))

### Inlet_BL ###

with open('Inlet_BLx','rb') as file:
    Inlet_BL_x=pickle.load(file)
with open('Inlet_BLy','rb') as file:
    Inlet_BL_y=pickle.load(file)
with open('Inlet_BLU','rb') as file:
    Inlet_BL_U=pickle.load(file)
with open('Inlet_BLV','rb') as file:
    Inlet_BL_V=pickle.load(file)
with open('Inlet_BLuu','rb') as file:
    Inlet_BL_uu=pickle.load(file)
with open('Inlet_BLvv','rb') as file:
    Inlet_BL_vv=pickle.load(file)
with open('Inlet_BLuv','rb') as file:
    Inlet_BL_uv=pickle.load(file)



### BL_w0_2D ###

with open('BL_w0_2Dx','rb') as file:
    BL_w0_2D_x=pickle.load(file)
with open('BL_w0_2Dy','rb') as file:
    BL_w0_2D_y=pickle.load(file)
with open('BL_w0_2DU','rb') as file:
    BL_w0_2D_U=pickle.load(file)
with open('BL_w0_2DV','rb') as file:
    BL_w0_2D_V=pickle.load(file)
with open('BL_w0_2Duu','rb') as file:
    BL_w0_2D_uu=pickle.load(file)
with open('BL_w0_2Dvv','rb') as file:
    BL_w0_2D_vv=pickle.load(file)
with open('BL_w0_2Duv','rb') as file:
    BL_w0_2D_uv=pickle.load(file)
   
### BL_w0_7D ###

with open('BL_w0_7Dx','rb') as file:
    BL_w0_7D_x=pickle.load(file)
with open('BL_w0_7Dy','rb') as file:
    BL_w0_7D_y=pickle.load(file)
with open('BL_w0_7DU','rb') as file:
    BL_w0_7D_U=pickle.load(file)
with open('BL_w0_7DV','rb') as file:
    BL_w0_7D_V=pickle.load(file)
with open('BL_w0_7Duu','rb') as file:
    BL_w0_7D_uu=pickle.load(file)
with open('BL_w0_7Dvv','rb') as file:
    BL_w0_7D_vv=pickle.load(file)
with open('BL_w0_7Duv','rb') as file:
    BL_w0_7D_uv=pickle.load(file)
   
       
### BL_w1_2D ###

with open('BL_w1_2Dx','rb') as file:
    BL_w1_2D_x=pickle.load(file)
with open('BL_w1_2Dy','rb') as file:
    BL_w1_2D_y=pickle.load(file)
with open('BL_w1_2DU','rb') as file:
    BL_w1_2D_U=pickle.load(file)
with open('BL_w1_2DV','rb') as file:
    BL_w1_2D_V=pickle.load(file)
with open('BL_w1_2Duu','rb') as file:
    BL_w1_2D_uu=pickle.load(file)
with open('BL_w1_2Dvv','rb') as file:
    BL_w1_2D_vv=pickle.load(file)
with open('BL_w1_2Duv','rb') as file:
    BL_w1_2D_uv=pickle.load(file)
   
### BL_w1_7D ###

with open('BL_w1_7Dx','rb') as file:
    BL_w1_7D_x=pickle.load(file)
with open('BL_w1_7Dy','rb') as file:
    BL_w1_7D_y=pickle.load(file)
with open('BL_w1_7DU','rb') as file:
    BL_w1_7D_U=pickle.load(file)
with open('BL_w1_7DV','rb') as file:
    BL_w1_7D_V=pickle.load(file)
with open('BL_w1_7Duu','rb') as file:
    BL_w1_7D_uu=pickle.load(file)
with open('BL_w1_7Dvv','rb') as file:
    BL_w1_7D_vv=pickle.load(file)
with open('BL_w1_7Duv','rb') as file:
    BL_w1_7D_uv=pickle.load(file)    
   
### BL_w2_2D ###

with open('BL_w2_2Dx','rb') as file:
    BL_w2_2D_x=pickle.load(file)
with open('BL_w2_2Dy','rb') as file:
    BL_w2_2D_y=pickle.load(file)
with open('BL_w2_2DU','rb') as file:
    BL_w2_2D_U=pickle.load(file)
with open('BL_w2_2DV','rb') as file:
    BL_w2_2D_V=pickle.load(file)
with open('BL_w2_2Duu','rb') as file:
    BL_w2_2D_uu=pickle.load(file)
with open('BL_w2_2Dvv','rb') as file:
    BL_w2_2D_vv=pickle.load(file)
with open('BL_w2_2Duv','rb') as file:
    BL_w2_2D_uv=pickle.load(file)
   
### BL_w2_7D ###

with open('BL_w2_7Dx','rb') as file:
    BL_w2_7D_x=pickle.load(file)
with open('BL_w2_7Dy','rb') as file:
    BL_w2_7D_y=pickle.load(file)
with open('BL_w2_7DU','rb') as file:
    BL_w2_7D_U=pickle.load(file)
with open('BL_w2_7DV','rb') as file:
    BL_w2_7D_V=pickle.load(file)
with open('BL_w2_7Duu','rb') as file:
    BL_w2_7D_uu=pickle.load(file)
with open('BL_w2_7Dvv','rb') as file:
    BL_w2_7D_vv=pickle.load(file)
with open('BL_w2_7Duv','rb') as file:
    BL_w2_7D_uv=pickle.load(file)    

### BL_w3_2D ###

with open('BL_w3_2Dx','rb') as file:
    BL_w3_2D_x=pickle.load(file)
with open('BL_w3_2Dy','rb') as file:
    BL_w3_2D_y=pickle.load(file)
with open('BL_w3_2DU','rb') as file:
    BL_w3_2D_U=pickle.load(file)
with open('BL_w3_2DV','rb') as file:
    BL_w3_2D_V=pickle.load(file)
with open('BL_w3_2Duu','rb') as file:
    BL_w3_2D_uu=pickle.load(file)
with open('BL_w3_2Dvv','rb') as file:
    BL_w3_2D_vv=pickle.load(file)
with open('BL_w3_2Duv','rb') as file:
    BL_w3_2D_uv=pickle.load(file)
   
### BL_w3_7D ###

with open('BL_w3_7Dx','rb') as file:
    BL_w3_7D_x=pickle.load(file)
with open('BL_w3_7Dy','rb') as file:
    BL_w3_7D_y=pickle.load(file)
with open('BL_w3_7DU','rb') as file:
    BL_w3_7D_U=pickle.load(file)
with open('BL_w3_7DV','rb') as file:
    BL_w3_7D_V=pickle.load(file)
with open('BL_w3_7Duu','rb') as file:
    BL_w3_7D_uu=pickle.load(file)
with open('BL_w3_7Dvv','rb') as file:
    BL_w3_7D_vv=pickle.load(file)
with open('BL_w3_7Duv','rb') as file:
    BL_w3_7D_uv=pickle.load(file)    



### LLJ_Inlet ###

with open('LLJ_Inletx','rb') as file:
    LLJ_Inlet_x=pickle.load(file)
with open('LLJ_Inlety','rb') as file:
    LLJ_Inlet_y=pickle.load(file)
with open('LLJ_InletU','rb') as file:
    LLJ_Inlet_U=pickle.load(file)
with open('LLJ_InletV','rb') as file:
    LLJ_Inlet_V=pickle.load(file)
with open('LLJ_Inletuu','rb') as file:
    LLJ_Inlet_uu=pickle.load(file)
with open('LLJ_Inletvv','rb') as file:
    LLJ_Inlet_vv=pickle.load(file)
with open('LLJ_Inletuv','rb') as file:
    LLJ_Inlet_uv=pickle.load(file)



### LLJ_POS_w0_2D ###

with open('LLJ_POS_w0_2Dx','rb') as file:
    LLJ_POS_w0_2D_x=pickle.load(file)
with open('LLJ_POS_w0_2Dy','rb') as file:
    LLJ_POS_w0_2D_y=pickle.load(file)
with open('LLJ_POS_w0_2DU','rb') as file:
    LLJ_POS_w0_2D_U=pickle.load(file)
with open('LLJ_POS_w0_2DV','rb') as file:
    LLJ_POS_w0_2D_V=pickle.load(file)
with open('LLJ_POS_w0_2Duu','rb') as file:
    LLJ_POS_w0_2D_uu=pickle.load(file)
with open('LLJ_POS_w0_2Dvv','rb') as file:
    LLJ_POS_w0_2D_vv=pickle.load(file)
with open('LLJ_POS_w0_2Duv','rb') as file:
    LLJ_POS_w0_2D_uv=pickle.load(file)
   
### LLJ_POS_w0_7D ###

with open('LLJ_POS_w0_7Dx','rb') as file:
    LLJ_POS_w0_7D_x=pickle.load(file)
with open('LLJ_POS_w0_7Dy','rb') as file:
    LLJ_POS_w0_7D_y=pickle.load(file)
with open('LLJ_POS_w0_7DU','rb') as file:
    LLJ_POS_w0_7D_U=pickle.load(file)
with open('LLJ_POS_w0_7DV','rb') as file:
    LLJ_POS_w0_7D_V=pickle.load(file)
with open('LLJ_POS_w0_7Duu','rb') as file:
    LLJ_POS_w0_7D_uu=pickle.load(file)
with open('LLJ_POS_w0_7Dvv','rb') as file:
    LLJ_POS_w0_7D_vv=pickle.load(file)
with open('LLJ_POS_w0_7Duv','rb') as file:
    LLJ_POS_w0_7D_uv=pickle.load(file)

### LLJ_POS_w1_2D ###

with open('LLJ_POS_w1_2Dx','rb') as file:
    LLJ_POS_w1_2D_x=pickle.load(file)
with open('LLJ_POS_w1_2Dy','rb') as file:
    LLJ_POS_w1_2D_y=pickle.load(file)
with open('LLJ_POS_w1_2DU','rb') as file:
    LLJ_POS_w1_2D_U=pickle.load(file)
with open('LLJ_POS_w1_2DV','rb') as file:
    LLJ_POS_w1_2D_V=pickle.load(file)
with open('LLJ_POS_w1_2Duu','rb') as file:
    LLJ_POS_w1_2D_uu=pickle.load(file)
with open('LLJ_POS_w1_2Dvv','rb') as file:
    LLJ_POS_w1_2D_vv=pickle.load(file)
with open('LLJ_POS_w1_2Duv','rb') as file:
    LLJ_POS_w1_2D_uv=pickle.load(file)
   
### LLJ_POS_w1_7D ###

with open('LLJ_POS_w1_7Dx','rb') as file:
    LLJ_POS_w1_7D_x=pickle.load(file)
with open('LLJ_POS_w1_7Dy','rb') as file:
    LLJ_POS_w1_7D_y=pickle.load(file)
with open('LLJ_POS_w1_7DU','rb') as file:
    LLJ_POS_w1_7D_U=pickle.load(file)
with open('LLJ_POS_w1_7DV','rb') as file:
    LLJ_POS_w1_7D_V=pickle.load(file)
with open('LLJ_POS_w1_7Duu','rb') as file:
    LLJ_POS_w1_7D_uu=pickle.load(file)
with open('LLJ_POS_w1_7Dvv','rb') as file:
    LLJ_POS_w1_7D_vv=pickle.load(file)
with open('LLJ_POS_w1_7Duv','rb') as file:
    LLJ_POS_w1_7D_uv=pickle.load(file)
   
### LLJ_POS_w2_2D ###

with open('LLJ_POS_w2_2Dx','rb') as file:
    LLJ_POS_w2_2D_x=pickle.load(file)
with open('LLJ_POS_w2_2Dy','rb') as file:
    LLJ_POS_w2_2D_y=pickle.load(file)
with open('LLJ_POS_w2_2DU','rb') as file:
    LLJ_POS_w2_2D_U=pickle.load(file)
with open('LLJ_POS_w2_2DV','rb') as file:
    LLJ_POS_w2_2D_V=pickle.load(file)
with open('LLJ_POS_w2_2Duu','rb') as file:
    LLJ_POS_w2_2D_uu=pickle.load(file)
with open('LLJ_POS_w2_2Dvv','rb') as file:
    LLJ_POS_w2_2D_vv=pickle.load(file)
with open('LLJ_POS_w2_2Duv','rb') as file:
    LLJ_POS_w2_2D_uv=pickle.load(file)
   
### LLJ_POS_w2_7D ###

#with open('LLJ_POS_w2_7Dx','rb') as file:
#    LLJ_POS_w2_7D_x=pickle.load(file)
#with open('LLJ_POS_w2_7Dy','rb') as file:
#    LLJ_POS_w2_7D_y=pickle.load(file)
#with open('LLJ_POS_w2_7DU','rb') as file:
#    LLJ_POS_w2_7D_U=pickle.load(file)
#with open('LLJ_POS_w2_7DV','rb') as file:
#    LLJ_POS_w2_7D_V=pickle.load(file)
#with open('LLJ_POS_w2_7Duu','rb') as file:
#    LLJ_POS_w2_7D_uu=pickle.load(file)
#with open('LLJ_POS_w2_7Dvv','rb') as file:
#    LLJ_POS_w2_7D_vv=pickle.load(file)
#with open('LLJ_POS_w2_7Duv','rb') as file:
#    LLJ_POS_w2_7D_uv=pickle.load(file)    

 ### LLJ_POS_w3_2D ###

with open('LLJ_POS_w3_2Dx','rb') as file:
    LLJ_POS_w3_2D_x=pickle.load(file)
with open('LLJ_POS_w3_2Dy','rb') as file:
    LLJ_POS_w3_2D_y=pickle.load(file)
with open('LLJ_POS_w3_2DU','rb') as file:
    LLJ_POS_w3_2D_U=pickle.load(file)
with open('LLJ_POS_w3_2DV','rb') as file:
    LLJ_POS_w3_2D_V=pickle.load(file)
with open('LLJ_POS_w3_2Duu','rb') as file:
    LLJ_POS_w3_2D_uu=pickle.load(file)
with open('LLJ_POS_w3_2Dvv','rb') as file:
    LLJ_POS_w3_2D_vv=pickle.load(file)
with open('LLJ_POS_w3_2Duv','rb') as file:
    LLJ_POS_w3_2D_uv=pickle.load(file)
   
## LLJ_POS_w3_7D ###

with open('LLJ_POS_w3_7Dx','rb') as file:
    LLJ_POS_w3_7D_x=pickle.load(file)
with open('LLJ_POS_w3_7Dy','rb') as file:
    LLJ_POS_w3_7D_y=pickle.load(file)
with open('LLJ_POS_w3_7DU','rb') as file:
    LLJ_POS_w3_7D_U=pickle.load(file)
with open('LLJ_POS_w3_7DV','rb') as file:
    LLJ_POS_w3_7D_V=pickle.load(file)
with open('LLJ_POS_w3_7Duu','rb') as file:
    LLJ_POS_w3_7D_uu=pickle.load(file)
with open('LLJ_POS_w3_7Dvv','rb') as file:
    LLJ_POS_w3_7D_vv=pickle.load(file)
with open('LLJ_POS_w3_7Duv','rb') as file:
    LLJ_POS_w3_7D_uv=pickle.load(file)    


### LLJ_HUB_w0_2D ###

with open('LLJ_HUB_w0_2Dx','rb') as file:
    LLJ_HUB_w0_2D_x=pickle.load(file)
with open('LLJ_HUB_w0_2Dy','rb') as file:
    LLJ_HUB_w0_2D_y=pickle.load(file)
with open('LLJ_HUB_w0_2DU','rb') as file:
    LLJ_HUB_w0_2D_U=pickle.load(file)
with open('LLJ_HUB_w0_2DV','rb') as file:
    LLJ_HUB_w0_2D_V=pickle.load(file)
with open('LLJ_HUB_w0_2Duu','rb') as file:
    LLJ_HUB_w0_2D_uu=pickle.load(file)
with open('LLJ_HUB_w0_2Dvv','rb') as file:
    LLJ_HUB_w0_2D_vv=pickle.load(file)
with open('LLJ_HUB_w0_2Duv','rb') as file:
    LLJ_HUB_w0_2D_uv=pickle.load(file)
   
### LLJ_HUB_w0_7D ###

with open('LLJ_HUB_w0_7Dx','rb') as file:
    LLJ_HUB_w0_7D_x=pickle.load(file)
with open('LLJ_HUB_w0_7Dy','rb') as file:
    LLJ_HUB_w0_7D_y=pickle.load(file)
with open('LLJ_HUB_w0_7DU','rb') as file:
    LLJ_HUB_w0_7D_U=pickle.load(file)
with open('LLJ_HUB_w0_7DV','rb') as file:
    LLJ_HUB_w0_7D_V=pickle.load(file)
with open('LLJ_HUB_w0_7Duu','rb') as file:
    LLJ_HUB_w0_7D_uu=pickle.load(file)
with open('LLJ_HUB_w0_7Dvv','rb') as file:
    LLJ_HUB_w0_7D_vv=pickle.load(file)
with open('LLJ_HUB_w0_7Duv','rb') as file:
    LLJ_HUB_w0_7D_uv=pickle.load(file)

### LLJ_HUB_w1_2D ###

with open('LLJ_HUB_w1_2Dx','rb') as file:
    LLJ_HUB_w1_2D_x=pickle.load(file)
with open('LLJ_HUB_w1_2Dy','rb') as file:
    LLJ_HUB_w1_2D_y=pickle.load(file)
with open('LLJ_HUB_w1_2DU','rb') as file:
    LLJ_HUB_w1_2D_U=pickle.load(file)
with open('LLJ_HUB_w1_2DV','rb') as file:
    LLJ_HUB_w1_2D_V=pickle.load(file)
with open('LLJ_HUB_w1_2Duu','rb') as file:
    LLJ_HUB_w1_2D_uu=pickle.load(file)
with open('LLJ_HUB_w1_2Dvv','rb') as file:
    LLJ_HUB_w1_2D_vv=pickle.load(file)
with open('LLJ_HUB_w1_2Duv','rb') as file:
    LLJ_HUB_w1_2D_uv=pickle.load(file)
   
### LLJ_HUB_w1_7D ###

with open('LLJ_HUB_w1_7Dx','rb') as file:
    LLJ_HUB_w1_7D_x=pickle.load(file)
with open('LLJ_HUB_w1_7Dy','rb') as file:
    LLJ_HUB_w1_7D_y=pickle.load(file)
with open('LLJ_HUB_w1_7DU','rb') as file:
    LLJ_HUB_w1_7D_U=pickle.load(file)
with open('LLJ_HUB_w1_7DV','rb') as file:
    LLJ_HUB_w1_7D_V=pickle.load(file)
with open('LLJ_HUB_w1_7Duu','rb') as file:
    LLJ_HUB_w1_7D_uu=pickle.load(file)
with open('LLJ_HUB_w1_7Dvv','rb') as file:
    LLJ_HUB_w1_7D_vv=pickle.load(file)
with open('LLJ_HUB_w1_7Duv','rb') as file:
    LLJ_HUB_w1_7D_uv=pickle.load(file)
   
### LLJ_HUB_w2_2D ###

with open('LLJ_HUB_w2_2Dx','rb') as file:
    LLJ_HUB_w2_2D_x=pickle.load(file)
with open('LLJ_HUB_w2_2Dy','rb') as file:
    LLJ_HUB_w2_2D_y=pickle.load(file)
with open('LLJ_HUB_w2_2DU','rb') as file:
    LLJ_HUB_w2_2D_U=pickle.load(file)
with open('LLJ_HUB_w2_2DV','rb') as file:
    LLJ_HUB_w2_2D_V=pickle.load(file)
with open('LLJ_HUB_w2_2Duu','rb') as file:
    LLJ_HUB_w2_2D_uu=pickle.load(file)
with open('LLJ_HUB_w2_2Dvv','rb') as file:
    LLJ_HUB_w2_2D_vv=pickle.load(file)
with open('LLJ_HUB_w2_2Duv','rb') as file:
    LLJ_HUB_w2_2D_uv=pickle.load(file)
   
### LLJ_HUB_w2_7D ###

with open('LLJ_HUB_w2_7Dx','rb') as file:
    LLJ_HUB_w2_7D_x=pickle.load(file)
with open('LLJ_HUB_w2_7Dy','rb') as file:
    LLJ_HUB_w2_7D_y=pickle.load(file)
with open('LLJ_HUB_w2_7DU','rb') as file:
    LLJ_HUB_w2_7D_U=pickle.load(file)
with open('LLJ_HUB_w2_7DV','rb') as file:
    LLJ_HUB_w2_7D_V=pickle.load(file)
with open('LLJ_HUB_w2_7Duu','rb') as file:
    LLJ_HUB_w2_7D_uu=pickle.load(file)
with open('LLJ_HUB_w2_7Dvv','rb') as file:
    LLJ_HUB_w2_7D_vv=pickle.load(file)
with open('LLJ_HUB_w2_7Duv','rb') as file:
    LLJ_HUB_w2_7D_uv=pickle.load(file)    


### LLJ_HUB_w3_2D ###

with open('LLJ_HUB_w3_2Dx','rb') as file:
    LLJ_HUB_w3_2D_x=pickle.load(file)
with open('LLJ_HUB_w3_2Dy','rb') as file:
    LLJ_HUB_w3_2D_y=pickle.load(file)
with open('LLJ_HUB_w3_2DU','rb') as file:
    LLJ_HUB_w3_2D_U=pickle.load(file)
with open('LLJ_HUB_w3_2DV','rb') as file:
    LLJ_HUB_w3_2D_V=pickle.load(file)
with open('LLJ_HUB_w3_2Duu','rb') as file:
    LLJ_HUB_w3_2D_uu=pickle.load(file)
with open('LLJ_HUB_w3_2Dvv','rb') as file:
    LLJ_HUB_w3_2D_vv=pickle.load(file)
with open('LLJ_HUB_w3_2Duv','rb') as file:
    LLJ_HUB_w3_2D_uv=pickle.load(file)
   
### LLJ_HUB_w3_7D ###

with open('LLJ_HUB_w3_7Dx','rb') as file:
    LLJ_HUB_w3_7D_x=pickle.load(file)
with open('LLJ_HUB_w3_7Dy','rb') as file:
    LLJ_HUB_w3_7D_y=pickle.load(file)
with open('LLJ_HUB_w3_7DU','rb') as file:
    LLJ_HUB_w3_7D_U=pickle.load(file)
with open('LLJ_HUB_w3_7DV','rb') as file:
    LLJ_HUB_w3_7D_V=pickle.load(file)
with open('LLJ_HUB_w3_7Duu','rb') as file:
    LLJ_HUB_w3_7D_uu=pickle.load(file)
with open('LLJ_HUB_w3_7Dvv','rb') as file:
    LLJ_HUB_w3_7D_vv=pickle.load(file)
with open('LLJ_HUB_w3_7Duv','rb') as file:
    LLJ_HUB_w3_7D_uv=pickle.load(file)      
   


## Normalization
#x = x/D

#y = y/H

##Delta U


#Domain reshape

x = np.reshape(BL_w0_2D_x,(154,159))
y = np.reshape(BL_w0_2D_y,(154,159))

U_BL_H = 7.85767360
U_LLJ_HUB_H = 8.37171
U_LLJ_POS_H = 6.74958

#### Domain correction (check with raw images)
   
y = y + 372.2
#y = y + 370
#x = x + (314.34+(5*D))



#Domain reshape

U_Inlet_BL = np.reshape(Inlet_BL_U,(154,159))
V_Inlet_BL = np.reshape(Inlet_BL_V,(154,159))
uu_Inlet_BL = np.reshape(Inlet_BL_uu,(154,159))
vv_Inlet_BL = np.reshape(Inlet_BL_vv,(154,159))
uv_Inlet_BL = np.reshape(Inlet_BL_uv,(154,159))

U_BL_w0_2D = np.reshape(BL_w0_2D_U,(154,159))
V_BL_w0_2D = np.reshape(BL_w0_2D_V,(154,159))
uu_BL_w0_2D = np.reshape(BL_w0_2D_uu,(154,159))
vv_BL_w0_2D = np.reshape(BL_w0_2D_vv,(154,159))
uv_BL_w0_2D = np.reshape(BL_w0_2D_uv,(154,159))

U_BL_w0_7D = np.reshape(BL_w0_7D_U,(154,159))
V_BL_w0_7D = np.reshape(BL_w0_7D_V,(154,159))
uu_BL_w0_7D = np.reshape(BL_w0_7D_uu,(154,159))
vv_BL_w0_7D = np.reshape(BL_w0_7D_vv,(154,159))
uv_BL_w0_7D = np.reshape(BL_w0_7D_uv,(154,159))

U_BL_w1_2D = np.reshape(BL_w1_2D_U,(154,159))
V_BL_w1_2D = np.reshape(BL_w1_2D_V,(154,159))
uu_BL_w1_2D = np.reshape(BL_w1_2D_uu,(154,159))
vv_BL_w1_2D = np.reshape(BL_w1_2D_vv,(154,159))
uv_BL_w1_2D = np.reshape(BL_w1_2D_uv,(154,159))

U_BL_w1_7D = np.reshape(BL_w1_7D_U,(154,159))
V_BL_w1_7D = np.reshape(BL_w1_7D_V,(154,159))
uu_BL_w1_7D = np.reshape(BL_w1_7D_uu,(154,159))
vv_BL_w1_7D = np.reshape(BL_w1_7D_vv,(154,159))
uv_BL_w1_7D = np.reshape(BL_w1_7D_uv,(154,159))

U_BL_w2_2D = np.reshape(BL_w2_2D_U,(154,159))
V_BL_w2_2D = np.reshape(BL_w2_2D_V,(154,159))
uu_BL_w2_2D = np.reshape(BL_w2_2D_uu,(154,159))
vv_BL_w2_2D = np.reshape(BL_w2_2D_vv,(154,159))
uv_BL_w2_2D = np.reshape(BL_w2_2D_uv,(154,159))

U_BL_w2_7D = np.reshape(BL_w2_7D_U,(154,159))
V_BL_w2_7D = np.reshape(BL_w2_7D_V,(154,159))
uu_BL_w2_7D = np.reshape(BL_w2_7D_uu,(154,159))
vv_BL_w2_7D = np.reshape(BL_w2_7D_vv,(154,159))
uv_BL_w2_7D = np.reshape(BL_w2_7D_uv,(154,159))

U_BL_w3_2D = np.reshape(BL_w3_2D_U,(154,159))
V_BL_w3_2D = np.reshape(BL_w3_2D_V,(154,159))
uu_BL_w3_2D = np.reshape(BL_w3_2D_uu,(154,159))
vv_BL_w3_2D = np.reshape(BL_w3_2D_vv,(154,159))
uv_BL_w3_2D = np.reshape(BL_w3_2D_uv,(154,159))

U_BL_w3_7D = np.reshape(BL_w3_7D_U,(154,159))
V_BL_w3_7D = np.reshape(BL_w3_7D_V,(154,159))
uu_BL_w3_7D = np.reshape(BL_w3_7D_uu,(154,159))
vv_BL_w3_7D = np.reshape(BL_w3_7D_vv,(154,159))
uv_BL_w3_7D = np.reshape(BL_w3_7D_uv,(154,159))



U_LLJ_Inlet = np.reshape(LLJ_Inlet_U,(154,159))
V_LLJ_Inlet = np.reshape(LLJ_Inlet_V,(154,159))
uu_LLJ_Inlet = np.reshape(LLJ_Inlet_uu,(154,159))
vv_LLJ_Inlet = np.reshape(LLJ_Inlet_vv,(154,159))
uv_LLJ_Inlet = np.reshape(LLJ_Inlet_uv,(154,159))


U_LLJ_POS_w0_2D = np.reshape(LLJ_POS_w0_2D_U,(154,159))
V_LLJ_POS_w0_2D = np.reshape(LLJ_POS_w0_2D_V,(154,159))
uu_LLJ_POS_w0_2D = np.reshape(LLJ_POS_w0_2D_uu,(154,159))
vv_LLJ_POS_w0_2D = np.reshape(LLJ_POS_w0_2D_vv,(154,159))
uv_LLJ_POS_w0_2D = np.reshape(LLJ_POS_w0_2D_uv,(154,159))

U_LLJ_POS_w0_7D = np.reshape(LLJ_POS_w0_7D_U,(154,159))
V_LLJ_POS_w0_7D = np.reshape(LLJ_POS_w0_7D_V,(154,159))
uu_LLJ_POS_w0_7D = np.reshape(LLJ_POS_w0_7D_uu,(154,159))
vv_LLJ_POS_w0_7D = np.reshape(LLJ_POS_w0_7D_vv,(154,159))
uv_LLJ_POS_w0_7D = np.reshape(LLJ_POS_w0_7D_uv,(154,159))


U_LLJ_POS_w1_2D = np.reshape(LLJ_POS_w1_2D_U,(154,159))
V_LLJ_POS_w1_2D = np.reshape(LLJ_POS_w1_2D_V,(154,159))
uu_LLJ_POS_w1_2D = np.reshape(LLJ_POS_w1_2D_uu,(154,159))
vv_LLJ_POS_w1_2D = np.reshape(LLJ_POS_w1_2D_vv,(154,159))
uv_LLJ_POS_w1_2D = np.reshape(LLJ_POS_w1_2D_uv,(154,159))

U_LLJ_POS_w1_7D = np.reshape(LLJ_POS_w1_7D_U,(154,159))
V_LLJ_POS_w1_7D = np.reshape(LLJ_POS_w1_7D_V,(154,159))
uu_LLJ_POS_w1_7D = np.reshape(LLJ_POS_w1_7D_uu,(154,159))
vv_LLJ_POS_w1_7D = np.reshape(LLJ_POS_w1_7D_vv,(154,159))
uv_LLJ_POS_w1_7D = np.reshape(LLJ_POS_w1_7D_uv,(154,159))


U_LLJ_POS_w2_2D = np.reshape(LLJ_POS_w2_2D_U,(154,159))
V_LLJ_POS_w2_2D = np.reshape(LLJ_POS_w2_2D_V,(154,159))
uu_LLJ_POS_w2_2D = np.reshape(LLJ_POS_w2_2D_uu,(154,159))
vv_LLJ_POS_w2_2D = np.reshape(LLJ_POS_w2_2D_vv,(154,159))
uv_LLJ_POS_w2_2D = np.reshape(LLJ_POS_w2_2D_uv,(154,159))

#U_LLJ_POS_w2_7D = np.reshape(LLJ_POS_w2_7D_U,(154,159))
#V_LLJ_POS_w2_7D = np.reshape(LLJ_POS_w2_7D_V,(154,159))
#uu_LLJ_POS_w2_7D = np.reshape(LLJ_POS_w2_7D_uu,(154,159))
#vv_LLJ_POS_w2_7D = np.reshape(LLJ_POS_w2_7D_vv,(154,159))
#uv_LLJ_POS_w2_7D = np.reshape(LLJ_POS_w2_7D_uv,(154,159))


U_LLJ_POS_w3_2D = np.reshape(LLJ_POS_w3_2D_U,(154,159))
V_LLJ_POS_w3_2D = np.reshape(LLJ_POS_w3_2D_V,(154,159))
uu_LLJ_POS_w3_2D = np.reshape(LLJ_POS_w3_2D_uu,(154,159))
vv_LLJ_POS_w3_2D = np.reshape(LLJ_POS_w3_2D_vv,(154,159))
uv_LLJ_POS_w3_2D = np.reshape(LLJ_POS_w3_2D_uv,(154,159))

U_LLJ_POS_w3_7D = np.reshape(LLJ_POS_w3_7D_U,(154,159))
V_LLJ_POS_w3_7D = np.reshape(LLJ_POS_w3_7D_V,(154,159))
uu_LLJ_POS_w3_7D = np.reshape(LLJ_POS_w3_7D_uu,(154,159))
vv_LLJ_POS_w3_7D = np.reshape(LLJ_POS_w3_7D_vv,(154,159))
uv_LLJ_POS_w3_7D = np.reshape(LLJ_POS_w3_7D_uv,(154,159))


U_LLJ_HUB_w0_2D = np.reshape(LLJ_HUB_w0_2D_U,(154,159))
V_LLJ_HUB_w0_2D = np.reshape(LLJ_HUB_w0_2D_V,(154,159))
uu_LLJ_HUB_w0_2D = np.reshape(LLJ_HUB_w0_2D_uu,(154,159))
vv_LLJ_HUB_w0_2D = np.reshape(LLJ_HUB_w0_2D_vv,(154,159))
uv_LLJ_HUB_w0_2D = np.reshape(LLJ_HUB_w0_2D_uv,(154,159))

U_LLJ_HUB_w0_7D = np.reshape(LLJ_HUB_w0_7D_U,(154,159))
V_LLJ_HUB_w0_7D = np.reshape(LLJ_HUB_w0_7D_V,(154,159))
uu_LLJ_HUB_w0_7D = np.reshape(LLJ_HUB_w0_7D_uu,(154,159))
vv_LLJ_HUB_w0_7D = np.reshape(LLJ_HUB_w0_7D_vv,(154,159))
uv_LLJ_HUB_w0_7D = np.reshape(LLJ_HUB_w0_7D_uv,(154,159))


U_LLJ_HUB_w1_2D = np.reshape(LLJ_HUB_w1_2D_U,(154,159))
V_LLJ_HUB_w1_2D = np.reshape(LLJ_HUB_w1_2D_V,(154,159))
uu_LLJ_HUB_w1_2D = np.reshape(LLJ_HUB_w1_2D_uu,(154,159))
vv_LLJ_HUB_w1_2D = np.reshape(LLJ_HUB_w1_2D_vv,(154,159))
uv_LLJ_HUB_w1_2D = np.reshape(LLJ_HUB_w1_2D_uv,(154,159))

U_LLJ_HUB_w1_7D = np.reshape(LLJ_HUB_w1_7D_U,(154,159))
V_LLJ_HUB_w1_7D = np.reshape(LLJ_HUB_w1_7D_V,(154,159))
uu_LLJ_HUB_w1_7D = np.reshape(LLJ_HUB_w1_7D_uu,(154,159))
vv_LLJ_HUB_w1_7D = np.reshape(LLJ_HUB_w1_7D_vv,(154,159))
uv_LLJ_HUB_w1_7D = np.reshape(LLJ_HUB_w1_7D_uv,(154,159))


U_LLJ_HUB_w2_2D = np.reshape(LLJ_HUB_w2_2D_U,(154,159))
V_LLJ_HUB_w2_2D = np.reshape(LLJ_HUB_w2_2D_V,(154,159))
uu_LLJ_HUB_w2_2D = np.reshape(LLJ_HUB_w2_2D_uu,(154,159))
vv_LLJ_HUB_w2_2D = np.reshape(LLJ_HUB_w2_2D_vv,(154,159))
uv_LLJ_HUB_w2_2D = np.reshape(LLJ_HUB_w2_2D_uv,(154,159))

U_LLJ_HUB_w2_7D = np.reshape(LLJ_HUB_w2_7D_U,(154,159))
V_LLJ_HUB_w2_7D = np.reshape(LLJ_HUB_w2_7D_V,(154,159))
uu_LLJ_HUB_w2_7D = np.reshape(LLJ_HUB_w2_7D_uu,(154,159))
vv_LLJ_HUB_w2_7D = np.reshape(LLJ_HUB_w2_7D_vv,(154,159))
uv_LLJ_HUB_w2_7D = np.reshape(LLJ_HUB_w2_7D_uv,(154,159))


U_LLJ_HUB_w3_2D = np.reshape(LLJ_HUB_w3_2D_U,(154,159))
V_LLJ_HUB_w3_2D = np.reshape(LLJ_HUB_w3_2D_V,(154,159))
uu_LLJ_HUB_w3_2D = np.reshape(LLJ_HUB_w3_2D_uu,(154,159))
vv_LLJ_HUB_w3_2D = np.reshape(LLJ_HUB_w3_2D_vv,(154,159))
uv_LLJ_HUB_w3_2D = np.reshape(LLJ_HUB_w3_2D_uv,(154,159))

U_LLJ_HUB_w3_7D = np.reshape(LLJ_HUB_w3_7D_U,(154,159))
V_LLJ_HUB_w3_7D = np.reshape(LLJ_HUB_w3_7D_V,(154,159))
uu_LLJ_HUB_w3_7D = np.reshape(LLJ_HUB_w3_7D_uu,(154,159))
vv_LLJ_HUB_w3_7D = np.reshape(LLJ_HUB_w3_7D_vv,(154,159))
uv_LLJ_HUB_w3_7D = np.reshape(LLJ_HUB_w3_7D_uv,(154,159))

### Crop Domain ###

x = x[5:148,20:130]
y = y[5:148,20:130]

### Domain Normalization ###

x1 = x/D
x2 = x1+5

y1 = (y-y_h)/D
y2 = (y-y_H)/D

### Crop Domain ###

U_Inlet_BL = U_Inlet_BL[5:148,20:130]
V_Inlet_BL = V_Inlet_BL[5:148,20:130]
uu_Inlet_BL = uu_Inlet_BL[5:148,20:130]
vv_Inlet_BL = vv_Inlet_BL[5:148,20:130]
uv_Inlet_BL = uv_Inlet_BL[5:148,20:130]

U_BL_w0_2D = U_BL_w0_2D[5:148,20:130]
V_BL_w0_2D = V_BL_w0_2D[5:148,20:130]
uu_BL_w0_2D = uu_BL_w0_2D[5:148,20:130]
vv_BL_w0_2D = vv_BL_w0_2D[5:148,20:130]
uv_BL_w0_2D = uv_BL_w0_2D[5:148,20:130]

U_BL_w0_7D = U_BL_w0_7D[5:148,20:130]
V_BL_w0_7D = V_BL_w0_7D[5:148,20:130]
uu_BL_w0_7D = uu_BL_w0_7D[5:148,20:130]
vv_BL_w0_7D = vv_BL_w0_7D[5:148,20:130]
uv_BL_w0_7D = uv_BL_w0_7D[5:148,20:130]

U_BL_w1_2D = U_BL_w1_2D[5:148,20:130]
V_BL_w1_2D = V_BL_w1_2D[5:148,20:130]
uu_BL_w1_2D = uu_BL_w1_2D[5:148,20:130]
vv_BL_w1_2D = vv_BL_w1_2D[5:148,20:130]
uv_BL_w1_2D = uv_BL_w1_2D[5:148,20:130]

U_BL_w1_7D = U_BL_w1_7D[5:148,20:130]
V_BL_w1_7D = V_BL_w1_7D[5:148,20:130]
uu_BL_w1_7D = uu_BL_w1_7D[5:148,20:130]
vv_BL_w1_7D = vv_BL_w1_7D[5:148,20:130]
uv_BL_w1_7D = uv_BL_w1_7D[5:148,20:130]

U_BL_w2_2D = U_BL_w2_2D[5:148,20:130]
V_BL_w2_2D = V_BL_w2_2D[5:148,20:130]
uu_BL_w2_2D = uu_BL_w2_2D[5:148,20:130]
vv_BL_w2_2D = vv_BL_w2_2D[5:148,20:130]
uv_BL_w2_2D = uv_BL_w2_2D[5:148,20:130]

U_BL_w2_7D = U_BL_w2_7D[5:148,20:130]
V_BL_w2_7D = V_BL_w2_7D[5:148,20:130]
uu_BL_w2_7D = uu_BL_w2_7D[5:148,20:130]
vv_BL_w2_7D = vv_BL_w2_7D[5:148,20:130]
uv_BL_w2_7D = uv_BL_w2_7D[5:148,20:130]

U_BL_w3_2D = U_BL_w3_2D[5:148,20:130]
V_BL_w3_2D = V_BL_w3_2D[5:148,20:130]
uu_BL_w3_2D = uu_BL_w3_2D[5:148,20:130]
vv_BL_w3_2D = vv_BL_w3_2D[5:148,20:130]
uv_BL_w3_2D = uv_BL_w3_2D[5:148,20:130]

U_BL_w3_7D = U_BL_w3_7D[5:148,20:130]
V_BL_w3_7D = V_BL_w3_7D[5:148,20:130]
uu_BL_w3_7D = uu_BL_w3_7D[5:148,20:130]
vv_BL_w3_7D = vv_BL_w3_7D[5:148,20:130]
uv_BL_w3_7D = uv_BL_w3_7D[5:148,20:130]

U_LLJ_Inlet = U_LLJ_Inlet[5:148,20:130]
V_LLJ_Inlet = V_LLJ_Inlet[5:148,20:130]
uu_LLJ_Inlet = uu_LLJ_Inlet[5:148,20:130]
vv_LLJ_Inlet = vv_LLJ_Inlet[5:148,20:130]
uv_LLJ_Inlet = uv_LLJ_Inlet[5:148,20:130]

U_LLJ_POS_w0_2D = U_LLJ_POS_w0_2D[5:148,20:130]
V_LLJ_POS_w0_2D = V_LLJ_POS_w0_2D[5:148,20:130]
uu_LLJ_POS_w0_2D = uu_LLJ_POS_w0_2D[5:148,20:130]
vv_LLJ_POS_w0_2D = vv_LLJ_POS_w0_2D[5:148,20:130]
uv_LLJ_POS_w0_2D = uv_LLJ_POS_w0_2D[5:148,20:130]

U_LLJ_POS_w0_7D = U_LLJ_POS_w0_7D[5:148,20:130]
V_LLJ_POS_w0_7D = V_LLJ_POS_w0_7D[5:148,20:130]
uu_LLJ_POS_w0_7D = uu_LLJ_POS_w0_7D[5:148,20:130]
vv_LLJ_POS_w0_7D = vv_LLJ_POS_w0_7D[5:148,20:130]
uv_LLJ_POS_w0_7D = uv_LLJ_POS_w0_7D[5:148,20:130]

U_LLJ_POS_w1_2D = U_LLJ_POS_w1_2D[5:148,20:130]
V_LLJ_POS_w1_2D = V_LLJ_POS_w1_2D[5:148,20:130]
uu_LLJ_POS_w1_2D = uu_LLJ_POS_w1_2D[5:148,20:130]
vv_LLJ_POS_w1_2D = vv_LLJ_POS_w1_2D[5:148,20:130]
uv_LLJ_POS_w1_2D = uv_LLJ_POS_w1_2D[5:148,20:130]

U_LLJ_POS_w1_7D = U_LLJ_POS_w1_7D[5:148,20:130]
V_LLJ_POS_w1_7D = V_LLJ_POS_w1_7D[5:148,20:130]
uu_LLJ_POS_w1_7D = uu_LLJ_POS_w1_7D[5:148,20:130]
vv_LLJ_POS_w1_7D = vv_LLJ_POS_w1_7D[5:148,20:130]
uv_LLJ_POS_w1_7D = uv_LLJ_POS_w1_7D[5:148,20:130]


U_LLJ_POS_w2_2D = U_LLJ_POS_w2_2D[5:148,20:130]
V_LLJ_POS_w2_2D = V_LLJ_POS_w2_2D[5:148,20:130]
uu_LLJ_POS_w2_2D = uu_LLJ_POS_w2_2D[5:148,20:130]
vv_LLJ_POS_w2_2D = vv_LLJ_POS_w2_2D[5:148,20:130]
uv_LLJ_POS_w2_2D = uv_LLJ_POS_w2_2D[5:148,20:130]


#U_LLJ_POS_w2_7D = U_LLJ_POS_w2_7D[5:148,20:130]
#V_LLJ_POS_w2_7D = V_LLJ_POS_w2_7D[5:148,20:130]
#uu_LLJ_POS_w2_7D = uu_LLJ_POS_w2_7D[5:148,20:130]
#vv_LLJ_POS_w2_7D = vv_LLJ_POS_w2_7D[5:148,20:130]
#uv_LLJ_POS_w2_7D = uv_LLJ_POS_w2_7D[5:148,20:130]

U_LLJ_POS_w3_2D = U_LLJ_POS_w3_2D[5:148,20:130]
V_LLJ_POS_w3_2D = V_LLJ_POS_w3_2D[5:148,20:130]
uu_LLJ_POS_w3_2D = uu_LLJ_POS_w3_2D[5:148,20:130]
vv_LLJ_POS_w3_2D = vv_LLJ_POS_w3_2D[5:148,20:130]
uv_LLJ_POS_w3_2D = uv_LLJ_POS_w3_2D[5:148,20:130]

U_LLJ_POS_w3_7D = U_LLJ_POS_w3_7D[5:148,20:130]
V_LLJ_POS_w3_7D = V_LLJ_POS_w3_7D[5:148,20:130]
uu_LLJ_POS_w3_7D = uu_LLJ_POS_w3_7D[5:148,20:130]
vv_LLJ_POS_w3_7D = vv_LLJ_POS_w3_7D[5:148,20:130]
uv_LLJ_POS_w3_7D = uv_LLJ_POS_w3_7D[5:148,20:130]



U_LLJ_HUB_w0_2D = U_LLJ_HUB_w0_2D[5:148,20:130]
V_LLJ_HUB_w0_2D = V_LLJ_HUB_w0_2D[5:148,20:130]
uu_LLJ_HUB_w0_2D = uu_LLJ_HUB_w0_2D[5:148,20:130]
vv_LLJ_HUB_w0_2D = vv_LLJ_HUB_w0_2D[5:148,20:130]
uv_LLJ_HUB_w0_2D = uv_LLJ_HUB_w0_2D[5:148,20:130]

U_LLJ_HUB_w0_7D = U_LLJ_HUB_w0_7D[5:148,20:130]
V_LLJ_HUB_w0_7D = V_LLJ_HUB_w0_7D[5:148,20:130]
uu_LLJ_HUB_w0_7D = uu_LLJ_HUB_w0_7D[5:148,20:130]
vv_LLJ_HUB_w0_7D = vv_LLJ_HUB_w0_7D[5:148,20:130]
uv_LLJ_HUB_w0_7D = uv_LLJ_HUB_w0_7D[5:148,20:130]

U_LLJ_HUB_w1_2D = U_LLJ_HUB_w1_2D[5:148,20:130]
V_LLJ_HUB_w1_2D = V_LLJ_HUB_w1_2D[5:148,20:130]
uu_LLJ_HUB_w1_2D = uu_LLJ_HUB_w1_2D[5:148,20:130]
vv_LLJ_HUB_w1_2D = vv_LLJ_HUB_w1_2D[5:148,20:130]
uv_LLJ_HUB_w1_2D = uv_LLJ_HUB_w1_2D[5:148,20:130]

U_LLJ_HUB_w1_7D = U_LLJ_HUB_w1_7D[5:148,20:130]
V_LLJ_HUB_w1_7D = V_LLJ_HUB_w1_7D[5:148,20:130]
uu_LLJ_HUB_w1_7D = uu_LLJ_HUB_w1_7D[5:148,20:130]
vv_LLJ_HUB_w1_7D = vv_LLJ_HUB_w1_7D[5:148,20:130]
uv_LLJ_HUB_w1_7D = uv_LLJ_HUB_w1_7D[5:148,20:130]


U_LLJ_HUB_w2_2D = U_LLJ_HUB_w2_2D[5:148,20:130]
V_LLJ_HUB_w2_2D = V_LLJ_HUB_w2_2D[5:148,20:130]
uu_LLJ_HUB_w2_2D = uu_LLJ_HUB_w2_2D[5:148,20:130]
vv_LLJ_HUB_w2_2D = vv_LLJ_HUB_w2_2D[5:148,20:130]
uv_LLJ_HUB_w2_2D = uv_LLJ_HUB_w2_2D[5:148,20:130]

U_LLJ_HUB_w2_7D = U_LLJ_HUB_w2_7D[5:148,20:130]
V_LLJ_HUB_w2_7D = V_LLJ_HUB_w2_7D[5:148,20:130]
uu_LLJ_HUB_w2_7D = uu_LLJ_HUB_w2_7D[5:148,20:130]
vv_LLJ_HUB_w2_7D = vv_LLJ_HUB_w2_7D[5:148,20:130]
uv_LLJ_HUB_w2_7D = uv_LLJ_HUB_w2_7D[5:148,20:130]

U_LLJ_HUB_w3_2D = U_LLJ_HUB_w3_2D[5:148,20:130]
V_LLJ_HUB_w3_2D = V_LLJ_HUB_w3_2D[5:148,20:130]
uu_LLJ_HUB_w3_2D = uu_LLJ_HUB_w3_2D[5:148,20:130]
vv_LLJ_HUB_w3_2D = vv_LLJ_HUB_w3_2D[5:148,20:130]
uv_LLJ_HUB_w3_2D = uv_LLJ_HUB_w3_2D[5:148,20:130]

U_LLJ_HUB_w3_7D = U_LLJ_HUB_w3_7D[5:148,20:130]
V_LLJ_HUB_w3_7D = V_LLJ_HUB_w3_7D[5:148,20:130]
uu_LLJ_HUB_w3_7D = uu_LLJ_HUB_w3_7D[5:148,20:130]
vv_LLJ_HUB_w3_7D = vv_LLJ_HUB_w3_7D[5:148,20:130]
uv_LLJ_HUB_w3_7D = uv_LLJ_HUB_w3_7D[5:148,20:130]


### Remove Outliers ###


P10 = np.percentile(uv_BL_w0_7D.flatten(),10)
P90 = np.percentile(uv_BL_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w0_7D[(uv_BL_w0_7D<lower_lim)|(uv_BL_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w1_2D.flatten(),10)
P90 = np.percentile(uv_BL_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w1_2D[(uv_BL_w1_2D<lower_lim)|(uv_BL_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w1_7D.flatten(),10)
P90 = np.percentile(uv_BL_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w1_7D[(uv_BL_w1_7D<lower_lim)|(uv_BL_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w2_2D.flatten(),10)
P90 = np.percentile(uv_BL_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w2_2D[(uv_BL_w2_2D<lower_lim)|(uv_BL_w2_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w2_7D.flatten(),10)
P90 = np.percentile(uv_BL_w2_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w2_7D[(uv_BL_w2_7D<lower_lim)|(uv_BL_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_BL_w3_2D.flatten(),10)
P90 = np.percentile(uv_BL_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w3_2D[(uv_BL_w3_2D<lower_lim)|(uv_BL_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w3_7D.flatten(),10)
P90 = np.percentile(uv_BL_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w3_7D[(uv_BL_w3_7D<lower_lim)|(uv_BL_w3_7D>upper_lim)] = np.nan


P10 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_2D[(uv_LLJ_HUB_w0_2D<lower_lim)|(uv_LLJ_HUB_w0_2D>upper_lim)] = np.nan
#df = pd.DataFrame(uv_LLJ_HUB_w0_2D)
#df.interpolate(method='nearest',axis=0,inplace=True)
#uv_LLJ_HUB_w0_2D = df.to_numpy()


P10 = np.percentile(uv_LLJ_POS_w0_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w0_7D[(uv_LLJ_POS_w0_7D<lower_lim)|(uv_LLJ_POS_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w1_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w1_2D[(uv_LLJ_POS_w1_2D<lower_lim)|(uv_LLJ_POS_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w1_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w1_7D[(uv_LLJ_POS_w1_7D<lower_lim)|(uv_LLJ_POS_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w2_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w2_2D[(uv_LLJ_POS_w2_2D<lower_lim)|(uv_LLJ_POS_w2_2D>upper_lim)] = np.nan

#P10 = np.percentile(uv_LLJ_POS_w2_7D.flatten(),10)
#P90 = np.percentile(uv_LLJ_POS_w2_7D.flatten(),90)
#IQR = P90-P10
#lower_lim = P10-1.0*IQR
#upper_lim = P90+1.0*IQR
#uv_LLJ_POS_w2_7D[(uv_LLJ_POS_w2_7D<lower_lim)|(uv_LLJ_POS_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_LLJ_POS_w3_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w3_2D[(uv_LLJ_POS_w3_2D<lower_lim)|(uv_LLJ_POS_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w3_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w3_7D[(uv_LLJ_POS_w3_7D<lower_lim)|(uv_LLJ_POS_w3_7D>upper_lim)] = np.nan




P10 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_2D[(uv_LLJ_HUB_w0_2D<lower_lim)|(uv_LLJ_HUB_w0_2D>upper_lim)] = np.nan
#df = pd.DataFrame(uv_LLJ_HUB_w0_2D)
#df.interpolate(method='nearest',axis=0,inplace=True)
#uv_LLJ_HUB_w0_2D = df.to_numpy()


P10 = np.percentile(uv_LLJ_HUB_w0_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_7D[(uv_LLJ_HUB_w0_7D<lower_lim)|(uv_LLJ_HUB_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w1_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w1_2D[(uv_LLJ_HUB_w1_2D<lower_lim)|(uv_LLJ_HUB_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w1_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w1_7D[(uv_LLJ_HUB_w1_7D<lower_lim)|(uv_LLJ_HUB_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w2_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w2_2D[(uv_LLJ_HUB_w2_2D<lower_lim)|(uv_LLJ_HUB_w2_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w2_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w2_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w2_7D[(uv_LLJ_HUB_w2_7D<lower_lim)|(uv_LLJ_HUB_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_LLJ_HUB_w3_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w3_2D[(uv_LLJ_HUB_w3_2D<lower_lim)|(uv_LLJ_HUB_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w3_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w3_7D[(uv_LLJ_HUB_w3_7D<lower_lim)|(uv_LLJ_HUB_w3_7D>upper_lim)] = np.nan



P10 = np.percentile(uv_BL_w0_7D.flatten(),10)
P90 = np.percentile(uv_BL_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w0_7D[(uv_BL_w0_7D<lower_lim)|(uv_BL_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w1_2D.flatten(),10)
P90 = np.percentile(uv_BL_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w1_2D[(uv_BL_w1_2D<lower_lim)|(uv_BL_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w1_7D.flatten(),10)
P90 = np.percentile(uv_BL_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w1_7D[(uv_BL_w1_7D<lower_lim)|(uv_BL_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w2_2D.flatten(),10)
P90 = np.percentile(uv_BL_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w2_2D[(uv_BL_w2_2D<lower_lim)|(uv_BL_w2_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w2_7D.flatten(),10)
P90 = np.percentile(uv_BL_w2_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w2_7D[(uv_BL_w2_7D<lower_lim)|(uv_BL_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_BL_w3_2D.flatten(),10)
P90 = np.percentile(uv_BL_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w3_2D[(uv_BL_w3_2D<lower_lim)|(uv_BL_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_BL_w3_7D.flatten(),10)
P90 = np.percentile(uv_BL_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_BL_w3_7D[(uv_BL_w3_7D<lower_lim)|(uv_BL_w3_7D>upper_lim)] = np.nan


P10 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_2D[(uv_LLJ_HUB_w0_2D<lower_lim)|(uv_LLJ_HUB_w0_2D>upper_lim)] = np.nan
#df = pd.DataFrame(uv_LLJ_HUB_w0_2D)
#df.interpolate(method='nearest',axis=0,inplace=True)
#uv_LLJ_HUB_w0_2D = df.to_numpy()


P10 = np.percentile(uv_LLJ_POS_w0_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w0_7D[(uv_LLJ_POS_w0_7D<lower_lim)|(uv_LLJ_POS_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w1_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w1_2D[(uv_LLJ_POS_w1_2D<lower_lim)|(uv_LLJ_POS_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w1_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w1_7D[(uv_LLJ_POS_w1_7D<lower_lim)|(uv_LLJ_POS_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w2_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w2_2D[(uv_LLJ_POS_w2_2D<lower_lim)|(uv_LLJ_POS_w2_2D>upper_lim)] = np.nan

#P10 = np.percentile(uv_LLJ_POS_w2_7D.flatten(),10)
#P90 = np.percentile(uv_LLJ_POS_w2_7D.flatten(),90)
#IQR = P90-P10
#lower_lim = P10-1.0*IQR
#upper_lim = P90+1.0*IQR
#uv_LLJ_POS_w2_7D[(uv_LLJ_POS_w2_7D<lower_lim)|(uv_LLJ_POS_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_LLJ_POS_w3_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w3_2D[(uv_LLJ_POS_w3_2D<lower_lim)|(uv_LLJ_POS_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_POS_w3_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_POS_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_POS_w3_7D[(uv_LLJ_POS_w3_7D<lower_lim)|(uv_LLJ_POS_w3_7D>upper_lim)] = np.nan




P10 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_2D[(uv_LLJ_HUB_w0_2D<lower_lim)|(uv_LLJ_HUB_w0_2D>upper_lim)] = np.nan
#df = pd.DataFrame(uv_LLJ_HUB_w0_2D)
#df.interpolate(method='nearest',axis=0,inplace=True)
#uv_LLJ_HUB_w0_2D = df.to_numpy()


P10 = np.percentile(uv_LLJ_HUB_w0_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w0_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w0_7D[(uv_LLJ_HUB_w0_7D<lower_lim)|(uv_LLJ_HUB_w0_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w1_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w1_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w1_2D[(uv_LLJ_HUB_w1_2D<lower_lim)|(uv_LLJ_HUB_w1_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w1_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w1_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w1_7D[(uv_LLJ_HUB_w1_7D<lower_lim)|(uv_LLJ_HUB_w1_7D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w2_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w2_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w2_2D[(uv_LLJ_HUB_w2_2D<lower_lim)|(uv_LLJ_HUB_w2_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w2_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w2_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w2_7D[(uv_LLJ_HUB_w2_7D<lower_lim)|(uv_LLJ_HUB_w2_7D>upper_lim)] = np.nan
   
P10 = np.percentile(uv_LLJ_HUB_w3_2D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w3_2D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w3_2D[(uv_LLJ_HUB_w3_2D<lower_lim)|(uv_LLJ_HUB_w3_2D>upper_lim)] = np.nan

P10 = np.percentile(uv_LLJ_HUB_w3_7D.flatten(),10)
P90 = np.percentile(uv_LLJ_HUB_w3_7D.flatten(),90)
IQR = P90-P10
lower_lim = P10-1.0*IQR
upper_lim = P90+1.0*IQR
uv_LLJ_HUB_w3_7D[(uv_LLJ_HUB_w3_7D<lower_lim)|(uv_LLJ_HUB_w3_7D>upper_lim)] = np.nan




### Inlet Profiles ###


U_Inlet_BL_profile = np.mean(U_Inlet_BL[:,50:60],axis=1)
U_LLJ_Inlet_profile = np.mean(U_LLJ_Inlet[:,50:60],axis=1)

TKE_Inlet_BL_profile = np.mean(0.5*(uu_Inlet_BL[:,50:60]+vv_Inlet_BL[:,50:60]),axis=1)
TKE_LLJ_Inlet_profile = np.mean(0.5*(uu_LLJ_Inlet[:,50:60]+vv_LLJ_Inlet[:,50:60]),axis=1)

U_BL = np.reshape(np.repeat(U_Inlet_BL_profile,110),(143,110))
U_LLJ = np.reshape(np.repeat(U_LLJ_Inlet_profile,110),(143,110))
y_profile = y[:,0]


### Velocity Deficit###

DeltaU_BL_w0_2D = (np.abs(U_BL_w0_2D-U_BL))/U_BL
DeltaU_BL_w1_2D = (np.abs(U_BL_w1_2D-U_BL))/U_BL
DeltaU_BL_w2_2D = (np.abs(U_BL_w2_2D-U_BL))/U_BL
DeltaU_BL_w3_2D = (np.abs(U_BL_w3_2D-U_BL))/U_BL

DeltaU_BL_w0_7D = (np.abs(U_BL_w0_7D-U_BL))/U_BL
DeltaU_BL_w1_7D = (np.abs(U_BL_w1_7D-U_BL))/U_BL
DeltaU_BL_w2_7D = (np.abs(U_BL_w2_7D-U_BL))/U_BL
DeltaU_BL_w3_7D = (np.abs(U_BL_w3_7D-U_BL))/U_BL

DeltaU_LLJ_POS_w0_2D = (np.abs(U_LLJ_POS_w0_2D-U_LLJ))/U_BL
DeltaU_LLJ_POS_w1_2D = (np.abs(U_LLJ_POS_w1_2D-U_LLJ))/U_BL
DeltaU_LLJ_POS_w2_2D = (np.abs(U_LLJ_POS_w2_2D-U_LLJ))/U_BL
DeltaU_LLJ_POS_w3_2D = (np.abs(U_LLJ_POS_w3_2D-U_LLJ))/U_BL

DeltaU_LLJ_POS_w0_7D = (np.abs(U_LLJ_POS_w0_7D-U_LLJ))/U_BL
DeltaU_LLJ_POS_w1_7D = (np.abs(U_LLJ_POS_w1_7D-U_LLJ))/U_BL
#DeltaU_LLJ_POS_w2_7D = (np.abs(U_LLJ_POS_w2_7D-U_LLJ_Inlet))/U_BL
DeltaU_LLJ_POS_w3_7D = (np.abs(U_LLJ_POS_w3_7D-U_LLJ))/U_BL

DeltaU_LLJ_HUB_w0_2D = (np.abs(U_LLJ_HUB_w0_2D-U_LLJ))/U_BL
DeltaU_LLJ_HUB_w1_2D = (np.abs(U_LLJ_HUB_w1_2D-U_LLJ))/U_BL
DeltaU_LLJ_HUB_w2_2D = (np.abs(U_LLJ_HUB_w2_2D-U_LLJ))/U_BL
DeltaU_LLJ_HUB_w3_2D = (np.abs(U_LLJ_HUB_w3_2D-U_LLJ))/U_BL

DeltaU_LLJ_HUB_w0_7D = (np.abs(U_LLJ_HUB_w0_7D-U_LLJ))/U_BL
DeltaU_LLJ_HUB_w1_7D = (np.abs(U_LLJ_HUB_w1_7D-U_LLJ))/U_BL
DeltaU_LLJ_HUB_w2_7D = (np.abs(U_LLJ_HUB_w2_7D-U_LLJ_Inlet))/U_BL
DeltaU_LLJ_HUB_w3_7D = (np.abs(U_LLJ_HUB_w3_7D-U_LLJ))/U_BL



### Tubulent Fluxes Normalized by the ###

flux_BL_w0_2D = (-U_BL_w0_2D*uv_BL_w0_2D)/U_BL**3
flux_BL_w1_2D = (-U_BL_w1_2D*uv_BL_w1_2D)/U_BL**3
flux_BL_w2_2D = (-U_BL_w2_2D*uv_BL_w2_2D)/U_BL**3
flux_BL_w3_2D = (-U_BL_w3_2D*uv_BL_w3_2D)/U_BL**3

flux_BL_w0_7D = (-U_BL_w0_7D*uv_BL_w0_7D)/U_BL**3
flux_BL_w1_7D = (-U_BL_w1_7D*uv_BL_w1_7D)/U_BL**3
flux_BL_w2_7D = (-U_BL_w2_7D*uv_BL_w2_7D)/U_BL**3
flux_BL_w3_7D = (-U_BL_w3_7D*uv_BL_w3_7D)/U_BL**3

flux_LLJ_POS_w0_2D = (-U_LLJ_POS_w0_2D*uv_LLJ_POS_w0_2D)/U_LLJ**3
flux_LLJ_POS_w1_2D = (-U_LLJ_POS_w1_2D*uv_LLJ_POS_w1_2D)/U_LLJ**3
flux_LLJ_POS_w2_2D = (-U_LLJ_POS_w2_2D*uv_LLJ_POS_w2_2D)/U_LLJ**3
flux_LLJ_POS_w3_2D = (-U_LLJ_POS_w3_2D*uv_LLJ_POS_w3_2D)/U_LLJ**3
#
flux_LLJ_POS_w0_7D = (-U_LLJ_POS_w0_7D*uv_LLJ_POS_w0_7D)/U_LLJ**3
flux_LLJ_POS_w1_7D = (-U_LLJ_POS_w1_7D*uv_LLJ_POS_w1_7D)/U_LLJ**3
#flux_LLJ_POS_w2_7D = (-U_LLJ_POS_w2_7D*uv_LLJ_POS_w2_7D)/U_LLJ**3
flux_LLJ_POS_w3_7D = (-U_LLJ_POS_w3_7D*uv_LLJ_POS_w3_7D)/U_LLJ**3


flux_LLJ_HUB_w0_2D = (-U_LLJ_HUB_w0_2D*uv_LLJ_HUB_w0_2D)/U_LLJ**3
flux_LLJ_HUB_w1_2D = (-U_LLJ_HUB_w1_2D*uv_LLJ_HUB_w1_2D)/U_LLJ**3
flux_LLJ_HUB_w2_2D = (-U_LLJ_HUB_w2_2D*uv_LLJ_HUB_w2_2D)/U_LLJ**3
flux_LLJ_HUB_w3_2D = (-U_LLJ_HUB_w3_2D*uv_LLJ_HUB_w3_2D)/U_LLJ**3
#
flux_LLJ_HUB_w0_7D = (-U_LLJ_HUB_w0_7D*uv_LLJ_HUB_w0_7D)/U_LLJ**3
flux_LLJ_HUB_w1_7D = (-U_LLJ_HUB_w1_7D*uv_LLJ_HUB_w1_7D)/U_LLJ**3
flux_LLJ_HUB_w2_7D = (-U_LLJ_HUB_w2_7D*uv_LLJ_HUB_w2_7D)/U_LLJ**3
flux_LLJ_HUB_w3_7D = (-U_LLJ_HUB_w3_7D*uv_LLJ_HUB_w3_7D)/U_LLJ**3





### Integral Turbulent Fluxes ###
sum_flux_BL_w0_2D = np.zeros(110)
sum_flux_BL_w1_2D = np.zeros(110)
sum_flux_BL_w2_2D = np.zeros(110)
sum_flux_BL_w3_2D = np.zeros(110)
sum_flux_BL_w0_7D = np.zeros(110)
sum_flux_BL_w1_7D = np.zeros(110)
sum_flux_BL_w2_7D = np.zeros(110)
sum_flux_BL_w3_7D = np.zeros(110)

sum_flux_LLJ_POS_w0_2D = np.zeros(110)
sum_flux_LLJ_POS_w1_2D = np.zeros(110)
sum_flux_LLJ_POS_w2_2D = np.zeros(110)
sum_flux_LLJ_POS_w3_2D = np.zeros(110)
sum_flux_LLJ_POS_w0_7D = np.zeros(110)
sum_flux_LLJ_POS_w1_7D = np.zeros(110)
sum_flux_LLJ_POS_w2_7D = np.zeros(110)
sum_flux_LLJ_POS_w3_7D = np.zeros(110)

sum_flux_LLJ_HUB_w0_2D = np.zeros(110)
sum_flux_LLJ_HUB_w1_2D = np.zeros(110)
sum_flux_LLJ_HUB_w2_2D = np.zeros(110)
sum_flux_LLJ_HUB_w3_2D = np.zeros(110)
sum_flux_LLJ_HUB_w0_7D = np.zeros(110)
sum_flux_LLJ_HUB_w1_7D = np.zeros(110)
sum_flux_LLJ_HUB_w2_7D = np.zeros(110)
sum_flux_LLJ_HUB_w3_7D = np.zeros(110)


for i in range(110):
    sum_flux_BL_w0_2D[i]=np.sum(flux_BL_w0_2D[63:111,i])
    sum_flux_BL_w1_2D[i]=np.sum(flux_BL_w1_2D[63:111,i])
    sum_flux_BL_w2_2D[i]=np.sum(flux_BL_w2_2D[63:111,i])
    sum_flux_BL_w3_2D[i]=np.sum(flux_BL_w3_2D[63:111,i])
    sum_flux_BL_w0_7D[i]=np.sum(flux_BL_w0_7D[63:111,i])
    sum_flux_BL_w1_7D[i]=np.sum(flux_BL_w1_7D[63:111,i])
    sum_flux_BL_w2_7D[i]=np.sum(flux_BL_w2_7D[63:111,i])
    sum_flux_BL_w3_7D[i]=np.sum(flux_BL_w3_7D[63:111,i])
    sum_flux_LLJ_POS_w0_2D[i]=np.sum(flux_LLJ_POS_w0_2D[63:111,i])
    sum_flux_LLJ_POS_w1_2D[i]=np.sum(flux_LLJ_POS_w1_2D[63:111,i])
    sum_flux_LLJ_POS_w2_2D[i]=np.sum(flux_LLJ_POS_w2_2D[63:111,i])
    sum_flux_LLJ_POS_w3_2D[i]=np.sum(flux_LLJ_POS_w3_2D[63:111,i])
    sum_flux_LLJ_POS_w0_7D[i]=np.sum(flux_LLJ_POS_w0_7D[63:111,i])
    sum_flux_LLJ_POS_w1_7D[i]=np.sum(flux_LLJ_POS_w1_7D[63:111,i])
    #sum_flux_LLJ_POS_w2_7D[i]=np.sum(flux_LLJ_POS_w2_7D[63:111,i])
    sum_flux_LLJ_POS_w3_7D[i]=np.sum(flux_LLJ_POS_w3_7D[63:111,i])
    sum_flux_LLJ_HUB_w0_2D[i]=np.sum(flux_LLJ_HUB_w0_2D[39:87,i])
    sum_flux_LLJ_HUB_w1_2D[i]=np.sum(flux_LLJ_HUB_w1_2D[39:87,i])
    sum_flux_LLJ_HUB_w2_2D[i]=np.sum(flux_LLJ_HUB_w2_2D[39:87,i])
    sum_flux_LLJ_HUB_w3_2D[i]=np.sum(flux_LLJ_HUB_w3_2D[39:87,i])
    sum_flux_LLJ_HUB_w0_7D[i]=np.sum(flux_LLJ_HUB_w0_7D[39:87,i])
    sum_flux_LLJ_HUB_w1_7D[i]=np.sum(flux_LLJ_HUB_w1_7D[39:87,i])
    sum_flux_LLJ_HUB_w2_7D[i]=np.sum(flux_LLJ_HUB_w2_7D[39:87,i])
    sum_flux_LLJ_HUB_w3_7D[i]=np.sum(flux_LLJ_HUB_w3_7D[39:87,i])



### Filter flux ###
#
#df = pd.DataFrame(sum_flux_BL_w0_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w0_2D = (df.to_numpy()).flatten()
#sum_flux_BL_w0_2D = signal.savgol_filter(sum_flux_BL_w0_2D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w0_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w0_7D = (df.to_numpy()).flatten()
#sum_flux_BL_w0_7D = signal.savgol_filter(sum_flux_BL_w0_7D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w1_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w1_2D = (df.to_numpy()).flatten()
#sum_flux_BL_w1_2D = signal.savgol_filter(sum_flux_BL_w1_2D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w1_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w1_7D = (df.to_numpy()).flatten()
#sum_flux_BL_w1_7D = signal.savgol_filter(sum_flux_BL_w1_7D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w2_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w2_2D = (df.to_numpy()).flatten()
#sum_flux_BL_w2_2D = signal.savgol_filter(sum_flux_BL_w2_2D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w2_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w2_7D = (df.to_numpy()).flatten()
#sum_flux_BL_w2_7D = signal.savgol_filter(sum_flux_BL_w2_7D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w3_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w3_2D = (df.to_numpy()).flatten()
#sum_flux_BL_w3_2D = signal.savgol_filter(sum_flux_BL_w3_2D,31,1)
#
#df = pd.DataFrame(sum_flux_BL_w3_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_BL_w3_7D = (df.to_numpy()).flatten()
#sum_flux_BL_w3_7D = signal.savgol_filter(sum_flux_BL_w3_7D,31,1)
#
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w0_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w0_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w0_2D = signal.savgol_filter(sum_flux_LLJ_POS_w0_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w0_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w0_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w0_7D = signal.savgol_filter(sum_flux_LLJ_POS_w0_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w1_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w1_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w1_2D = signal.savgol_filter(sum_flux_LLJ_POS_w1_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w1_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w1_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w1_7D = signal.savgol_filter(sum_flux_LLJ_POS_w1_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w2_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w2_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w2_2D = signal.savgol_filter(sum_flux_LLJ_POS_w2_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w2_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w2_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w2_7D = signal.savgol_filter(sum_flux_LLJ_POS_w2_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w3_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w3_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w3_2D = signal.savgol_filter(sum_flux_LLJ_POS_w3_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_POS_w3_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_POS_w3_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_POS_w3_7D = signal.savgol_filter(sum_flux_LLJ_POS_w3_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w0_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w0_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w0_2D = signal.savgol_filter(sum_flux_LLJ_HUB_w0_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w0_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w0_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w0_7D = signal.savgol_filter(sum_flux_LLJ_HUB_w0_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w1_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w1_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w1_2D = signal.savgol_filter(sum_flux_LLJ_HUB_w1_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w1_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w1_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w1_7D = signal.savgol_filter(sum_flux_LLJ_HUB_w1_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w2_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w2_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w2_2D = signal.savgol_filter(sum_flux_LLJ_HUB_w2_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w2_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w2_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w2_7D = signal.savgol_filter(sum_flux_LLJ_HUB_w2_7D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w3_2D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w3_2D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w3_2D = signal.savgol_filter(sum_flux_LLJ_HUB_w3_2D,31,1)
#
#df = pd.DataFrame(sum_flux_LLJ_HUB_w3_7D)
#df.interpolate(method='linear',inplace=True)
#sum_flux_LLJ_HUB_w3_7D = (df.to_numpy()).flatten()
#sum_flux_LLJ_HUB_w3_7D = signal.savgol_filter(sum_flux_LLJ_HUB_w3_7D,31,1)


### Figures ###

fig=plt.figure(figsize=(10, 4))
###
ax1=fig.add_subplot(131)
##
ax1.plot(x1[0,:],sum_flux_BL_w0_2D,label='0% R',color='k')
ax1.plot(x1[0,:],sum_flux_BL_w1_2D,label='9% R',color='C0')
ax1.plot(x1[0,:],sum_flux_BL_w2_2D,label='10% R',color='C1')
ax1.plot(x1[0,:],sum_flux_BL_w3_2D,label='11% R',color='C2')

ax1.plot(x2[0,:],sum_flux_BL_w0_7D,color='k')
ax1.plot(x2[0,:],sum_flux_BL_w1_7D,color='C0')
ax1.plot(x2[0,:],sum_flux_BL_w2_7D,color='C1')
ax1.plot(x2[0,:],sum_flux_BL_w3_7D,color='C2')

ax1.set_ylim(-0.01,0.11)    
ax1.set_ylabel(r'$\int_{bt}^{tt}{{-U\overline{u^\prime v^\prime}}/{U_{inlet}}^3}$'
               ,fontsize=14)
ax1.ticklabel_format(axis='y',style='sci',scilimits=(-2,-2),useMathText=True)
ax1.set_yticks([0,0.02,0.04,0.06,0.08,0.1])
ax1.set_xlabel('$x/D$',fontsize=14)
ax1.set_xlim(0,8)
ax1.legend(fontsize=12,loc='lower right',title='Winglet Span')
ax1.set_title('(a)',size='16')
 
ax2=fig.add_subplot(132)  
#
ax2.plot(x1[0,:],sum_flux_LLJ_POS_w0_2D,label='0% R',color='k')
ax2.plot(x1[0,:],sum_flux_LLJ_POS_w1_2D,label='9% R',color='C0')
ax2.plot(x1[0,:],sum_flux_LLJ_POS_w2_2D,label='10% R',color='C1')
ax2.plot(x1[0,:],sum_flux_LLJ_POS_w3_2D,label='11% R',color='C2')

ax2.plot(x2[0,:],sum_flux_LLJ_POS_w0_7D,color='k')
ax2.plot(x2[0,:],sum_flux_LLJ_POS_w1_7D,color='C0')
#ax2.plot(x2[0,:],sum_flux_LLJ_POS_w2_7D,color='C1')
ax2.plot(x2[0,:],sum_flux_LLJ_POS_w3_7D,color='C2')

ax2.set_ylim(-0.01,0.11)
ax2.set_yticks([])
ax2.ticklabel_format(axis='y',style='sci',scilimits=(-2,-2),useMathText=True)
ax2.set_yticks([0,0.02,0.04,0.06,0.08,0.1])
ax2.set_xlabel('$x/D$',fontsize=14)
#ax2.legend(fontsize=12,loc='lower right')
ax2.set_xlim(0,8)
ax2.set_title('(b)',size='16')
 

ax3=fig.add_subplot(133)  

ax3.plot(x1[0,:],sum_flux_LLJ_HUB_w0_2D,label='w0',color='k')
ax3.plot(x1[0,:],sum_flux_LLJ_HUB_w1_2D,label='w1',color='C0')
ax3.plot(x1[0,:],sum_flux_LLJ_HUB_w2_2D,label='w2',color='C1')
ax3.plot(x1[0,:],sum_flux_LLJ_HUB_w3_2D,label='w3',color='C2')

ax3.plot(x2[0,:],sum_flux_LLJ_HUB_w0_7D,color='k')
ax3.plot(x2[0,:],sum_flux_LLJ_HUB_w1_7D,color='C0')
ax3.plot(x2[0,:],sum_flux_LLJ_HUB_w2_7D,color='C1')
ax3.plot(x2[0,:],sum_flux_LLJ_HUB_w3_7D,color='C2')

ax3.set_ylim(-0.01,0.11)
ax3.set_yticks([])
ax3.ticklabel_format(axis='y',style='sci',scilimits=(-2,-2),useMathText=True)
ax3.set_yticks([0,0.02,0.04,0.06,0.08,0.1])
ax3.set_xlabel('$x/D$',fontsize=14)
ax3.set_xlim(0,8)  
ax3.set_title('(c)',size='16')
 
plt.savefig('sum_flux.png',dpi=300)
plt.savefig('sum_flux.eps',dpi=300)
plt.savefig('sum_flux.svg',dpi=300)

#
#plt.ylim(-0.02,0.12)
##plt.yticks([0,0.04,0.08,0.12])
#plt.xlabel('$x/D$')
#plt.ylabel(r'$\int{{-U\overline{u^\prime v^\prime}}/{U_{inlet}}^3}$')
#
#plt.legend(fontsize=12)  
#plt.tight_layout()


#plt.savefig('sum_flux_BL.png',dpi=300)
#plt.savefig('sum_flux_BL.eps',dpi=300)
#plt.savefig('sum_flux_BL.svg',dpi=300)

#plt.savefig('sum_flux_LLJ_POS.png',dpi=300)
#plt.savefig('sum_flux_LLJ_POS.eps',dpi=300)
#plt.savefig('sum_flux_LLJ_POS.svg',dpi=300)
#
#plt.savefig('sum_flux_LLJ_HUB.png',dpi=300)
#plt.savefig('sum_flux_LLJ_HUB.eps',dpi=300)
#plt.savefig('sum_flux_LLJ_HUB.svg',dpi=300)


#### Turbulent Vertical Flux ####
#
#fig1=plt.figure(figsize=(6, 10))
###
#ax1=fig1.add_subplot(421)
#ax1.contourf(x1,y1,flux_BL_w0_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig1.add_subplot(422)
#ax2.contourf(x2,y1,flux_BL_w0_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig1.add_subplot(423)
#ax3.contourf(x1,y1,flux_BL_w1_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig1.add_subplot(424)
#ax4.contourf(x2,y1,flux_BL_w1_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig1.add_subplot(425)
#ax5.contourf(x1,y1,flux_BL_w2_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig1.add_subplot(426)
#ax6.contourf(x2,y1,flux_BL_w2_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig1.add_subplot(427)
#ax7.contourf(x1,y1,flux_BL_w3_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig1.add_subplot(428)
#ax8.contourf(x2,y1,flux_BL_w3_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('flux_BL.png',dpi=300)
#plt.savefig('flux_BL.eps',dpi=300)
#plt.savefig('flux_BL.svg',dpi=300)
#
#
#fig2=plt.figure(figsize=(6, 10))
###
#ax1=fig2.add_subplot(421)
#ax1.contourf(x1,y1,flux_LLJ_POS_w0_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig2.add_subplot(422)
#ax2.contourf(x2,y1,flux_LLJ_POS_w0_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig2.add_subplot(423)
#ax3.contourf(x1,y1,flux_LLJ_POS_w1_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig2.add_subplot(424)
#ax4.contourf(x2,y1,flux_LLJ_POS_w1_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig2.add_subplot(425)
#ax5.contourf(x1,y1,flux_LLJ_POS_w2_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig2.add_subplot(426)
##ax6.contourf(x2,y1,flux_LLJ_POS_w2_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax6.contourf(x2,y1,zeros,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig2.add_subplot(427)
#ax7.contourf(x1,y1,flux_LLJ_POS_w3_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig2.add_subplot(428)
#ax8.contourf(x2,y1,flux_LLJ_POS_w3_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('flux_LLJ_POS.png',dpi=300)
#plt.savefig('flux_LLJ_POS.eps',dpi=300)
#plt.savefig('flux_LLJ_POS.svg',dpi=300)
#
#
#
#fig3=plt.figure(figsize=(6, 10))
###
#ax1=fig3.add_subplot(421)
#ax1.contourf(x1,y2,flux_LLJ_HUB_w0_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig3.add_subplot(422)
#ax2.contourf(x2,y2,flux_LLJ_HUB_w0_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig3.add_subplot(423)
#ax3.contourf(x1,y2,flux_LLJ_HUB_w1_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig3.add_subplot(424)
#ax4.contourf(x2,y2,flux_LLJ_HUB_w1_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig3.add_subplot(425)
#ax5.contourf(x1,y2,flux_LLJ_HUB_w2_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig3.add_subplot(426)
#ax6.contourf(x2,y2,flux_LLJ_HUB_w2_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
##ax6.contourf(x2,y2,zeros,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig3.add_subplot(427)
#ax7.contourf(x1,y2,flux_LLJ_HUB_w3_2D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
##ax7.contourf(x2,y2,zeros,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig3.add_subplot(428)
#ax8.contourf(x2,y2,flux_LLJ_HUB_w3_7D,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
##ax8.contourf(x2,y2,zeros,cmap='seismic',levels=np.linspace(-0.01,0.01,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('flux_LLJ_HUB.png',dpi=300)
#plt.savefig('flux_LLJ_HUB.eps',dpi=300)
#plt.savefig('flux_LLJ_HUB.svg',dpi=300)
#
#
#### Velocity Deficit ####
#
#fig4=plt.figure(figsize=(6, 10))
###
#ax1=fig4.add_subplot(421)
#ax1.contourf(x1,y1,DeltaU_BL_w0_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig4.add_subplot(422)
#ax2.contourf(x2,y1,DeltaU_BL_w0_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig4.add_subplot(423)
#ax3.contourf(x1,y1,DeltaU_BL_w1_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig4.add_subplot(424)
#ax4.contourf(x2,y1,DeltaU_BL_w1_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig4.add_subplot(425)
#ax5.contourf(x1,y1,DeltaU_BL_w2_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig4.add_subplot(426)
#ax6.contourf(x2,y1,DeltaU_BL_w2_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig4.add_subplot(427)
#ax7.contourf(x1,y1,DeltaU_BL_w3_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig4.add_subplot(428)
#ax8.contourf(x2,y1,DeltaU_BL_w3_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('DeltaU_BL.png',dpi=300)
#plt.savefig('DeltaU_BL.eps',dpi=300)
#plt.savefig('DeltaU_BL.svg',dpi=300)
#
#
#fig5=plt.figure(figsize=(6, 10))
###
#ax1=fig5.add_subplot(421)
#ax1.contourf(x1,y1,DeltaU_LLJ_POS_w0_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig5.add_subplot(422)
#ax2.contourf(x2,y1,DeltaU_LLJ_POS_w0_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig5.add_subplot(423)
#ax3.contourf(x1,y1,DeltaU_LLJ_POS_w1_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig5.add_subplot(424)
#ax4.contourf(x2,y1,DeltaU_LLJ_POS_w1_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig5.add_subplot(425)
#ax5.contourf(x1,y1,DeltaU_LLJ_POS_w2_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig5.add_subplot(426)
##ax6.contourf(x2,y1,DeltaU_LLJ_POS_w2_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax6.contourf(x2,y1,zeros,cmap='jet',levels=np.linspace(0.0,0.01,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig5.add_subplot(427)
#ax7.contourf(x1,y1,DeltaU_LLJ_POS_w3_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig5.add_subplot(428)
#ax8.contourf(x2,y1,DeltaU_LLJ_POS_w3_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('DeltaU_LLJ_POS.png',dpi=300)
#plt.savefig('DeltaU_LLJ_POS.eps',dpi=300)
#plt.savefig('DeltaU_LLJ_POS.svg',dpi=300)
#
#
#
#fig6=plt.figure(figsize=(6, 10))
###
#ax1=fig6.add_subplot(421)
#ax1.contourf(x1,y2,DeltaU_LLJ_HUB_w0_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax1.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax1.axis('scaled')
#ax1.set_ylim(-0.75,0.75)
#ax1.set_ylabel('$(y-y_H)/D$')
#ax1.set_yticks([-0.5,0.0,0.5])
#ax1.set_title('(a)',size='12')
##
##
#ax2=fig6.add_subplot(422)
#ax2.contourf(x2,y2,DeltaU_LLJ_HUB_w0_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax2.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax2.axis('scaled')
#ax2.set_ylim(-0.75,0.75)
#ax2.set_yticks([])
#ax2.set_title('(b)',size='12')
##
##
#ax3=fig6.add_subplot(423)
#ax3.contourf(x1,y2,DeltaU_LLJ_HUB_w1_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax3.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax3.axis('scaled')
#ax3.set_ylim(-0.75,0.75)
#ax3.set_ylabel('$(y-y_H)/D$')
#ax3.set_yticks([-0.5,0.0,0.5])
#ax3.set_title('(c)',size='12')
##
##
#ax4=fig6.add_subplot(424)
#ax4.contourf(x2,y2,DeltaU_LLJ_HUB_w1_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax4.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax4.axis('scaled')
#ax4.set_ylim(-0.75,0.75)
#ax4.set_yticks([])
#ax4.set_title('(d)',size='12')
##
##
#ax5=fig6.add_subplot(425)
#ax5.contourf(x1,y2,DeltaU_LLJ_HUB_w2_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax5.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax5.axis('scaled')
#ax5.set_ylim(-0.75,0.75)
#ax5.set_ylabel('$(y-y_H)/D$')
#ax5.set_yticks([-0.5,0.0,0.5])
#ax5.set_title('(e)',size='12')
##
##
#ax6=fig6.add_subplot(426)
#ax6.contourf(x2,y2,DeltaU_LLJ_HUB_w2_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
##ax6.contourf(x2,y2,zeros,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax6.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax6.axis('scaled')
#ax6.set_ylim(-0.75,0.75)
#ax6.set_yticks([])
#ax6.set_title('(f)',size='12')
##
##
#ax7=fig6.add_subplot(427)
#ax7.contourf(x1,y2,DeltaU_LLJ_HUB_w3_2D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
##ax7.contourf(x2,y2,zeros,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax7.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax7.axis('scaled')
#ax7.set_ylim(-0.75,0.75)
#ax7.set_xlabel('$x/D$')
#ax7.set_ylabel('$(y-y_H)/D$')
#ax7.set_yticks([-0.5,0.0,0.5])
#ax7.set_title('(g)',size='12')
##
##
#ax8=fig6.add_subplot(428)
#ax8.contourf(x2,y2,DeltaU_LLJ_HUB_w3_7D,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
##ax8.contourf(x2,y2,zeros,cmap='jet',levels=np.linspace(0.0,0.6,50),extend='both')
#ax8.axhline(y=-0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axhline(y=0.5,linestyle='--',color='k',linewidth=0.8)
#ax8.axis('scaled')
#ax8.set_ylim(-0.75,0.75)
#ax8.set_xlabel('$x/D$')
#ax8.set_yticks([])
#ax8.set_title('(h)',size='12')
##
##
#plt.savefig('DeltaU_LLJ_HUB.png',dpi=300)
#plt.savefig('DeltaU_LLJ_HUB.eps',dpi=300)
#plt.savefig('DeltaU_LLJ_HUB.svg',dpi=300)


#X_ = X[0,:]
#Y_ = np.interp(X_,x_,y_)
#
#
#inlet_U = np.reshape(hill_1_base_U,(154,206))
#inlet_U = np.delete(inlet_U,(0,1,205),axis=1)
#
#inlet_U[inlet_U==0]=['nan']
#
#inlet_V = np.reshape(hill_1_base_V,(154,206))
#inlet_V = np.delete(inlet_V,(0,1,205),axis=1)
#
#U_hub = np.zeros_like(X_)
#
#for i in range(len(X_)):
#    U_hub[i] = np.interp((Y_[i]+1),(Y[:,i]),(inlet_U[:,i]))
#    
#
#hill_1_one_turbine_U = np.reshape(hill_1_one_turbine_U,(154,206))
#hill_1_one_turbine_U = np.delete(hill_1_one_turbine_U,(0,1,205),axis=1)
#
#hill_1_one_turbine_V = np.reshape(hill_1_one_turbine_V,(154,206))
#
#hill_1_one_turbine_V = np.delete(hill_1_one_turbine_V,(0,1,205),axis=1)
#
#
#hill_1_two_turbine_U = np.reshape(hill_1_two_turbine_U,(154,206))
#hill_1_two_turbine_U = np.delete(hill_1_two_turbine_U,(0,1,205),axis=1)
#
#hill_1_two_turbine_V = np.reshape(hill_1_two_turbine_V,(154,206))
#hill_1_two_turbine_V = np.delete(hill_1_two_turbine_V,(0,1,205),axis=1)
#
#hill_1_one_turbine_uv = np.reshape(hill_1_one_turbine_uv,(154,206))
#hill_1_one_turbine_uv = np.delete(hill_1_one_turbine_uv,(0,1,205),axis=1)
#
#hill_1_two_turbine_uv = np.reshape(hill_1_two_turbine_uv,(154,206))
#hill_1_two_turbine_uv = np.delete(hill_1_two_turbine_uv,(0,1,205),axis=1)
#
#deltaU_one_turbine = np.abs(inlet_U - hill_1_one_turbine_U)/inlet_U
#deltaU_two_turbine = np.abs(inlet_U - hill_1_two_turbine_U)/inlet_U


##Gradients
#
#dx=(abs(X[0,0]-X[0,1])*D)/1000
#dy=(abs(Y[0,0]-Y[1,0])*D)/1000
##
#
#grady_U_base = np.gradient(inlet_U,dy,axis=0)
#grady_U_one_turbine = np.gradient(hill_1_one_turbine_U,dy,axis=0)
#grady_U_two_turbine = np.gradient(hill_1_two_turbine_U,dy,axis=0)
#
#grady_uv_one_turbine = np.gradient(hill_1_one_turbine_uv,dy,axis=0)
#grady_uv_two_turbine = np.gradient(hill_1_two_turbine_uv,dy,axis=0)
#
#

##Gradients
#
#dx=(abs(X[0,0]-X[0,1])*D)/1000
#dy=(abs(Y[0,0]-Y[1,0])*D)/1000
##
#
#grady_U_base = np.gradient(inlet_U,dy,axis=0)
#grady_U_one_turbine = np.gradient(hill_1_one_turbine_U,dy,axis=0)
#grady_U_two_turbine = np.gradient(hill_1_two_turbine_U,dy,axis=0)
#
#grady_uv_one_turbine = np.gradient(hill_1_one_turbine_uv,dy,axis=0)
#grady_uv_two_turbine = np.gradient(hill_1_two_turbine_uv,dy,axis=0)
#
#
#adv_one_turbine = (-hill_1_one_turbine_V*grady_U_one_turbine)*((D/1000)/inlet_U**2)
#adv_two_turbine = (-hill_1_two_turbine_V*grady_U_two_turbine)*((D/1000)/inlet_U**2)
#
#tur_one_turbine = (-grady_uv_one_turbine)*((D/1000)/inlet_U**2)
#tur_two_turbine = (-grady_uv_two_turbine)*((D/1000)/inlet_U**2)
#
#
#turflux_one_turbine = (-hill_1_one_turbine_U*hill_1_one_turbine_uv)*(1/inlet_U**3)
#turflux_two_turbine = (-hill_1_two_turbine_U*hill_1_two_turbine_uv)*(1/inlet_U**3)
#
#advflux_one_turbine = (0.5*((hill_1_one_turbine_U)**2)*hill_1_one_turbine_V)*(1/inlet_U**3)
#advflux_two_turbine = (0.5*((hill_1_two_turbine_U)**2)*hill_1_two_turbine_V)*(1/inlet_U**3)
#

### Momentum Integration ###

#y_ground = np.interp(X[0,:],x_,y_)
#y_aux = np.arange(0.5,1.505,0.005)
#
#
#sum_adv_one_turbine = np.zeros_like(X[0,:])
#sum_adv_two_turbine = np.zeros_like(X[0,:])
#sum_tur_one_turbine = np.zeros_like(X[0,:])
#sum_tur_two_turbine = np.zeros_like(X[0,:])
##
#for i in range (len(X[0,:])):
#    f = interpolate.interp1d(Y[:,i],adv_one_turbine[:,i])
#    sum_adv_one_turbine[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],adv_two_turbine[:,i])
#    sum_adv_two_turbine[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],tur_one_turbine[:,i])
#    sum_tur_one_turbine[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],tur_two_turbine[:,i])
#    sum_tur_two_turbine[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#
#np.savetxt('hill1_sum_adv_one_turbine.txt',sum_adv_one_turbine)
#np.savetxt('hill1_sum_adv_two_turbine.txt',sum_adv_two_turbine)
#np.savetxt('hill1_sum_tur_one_turbine.txt',sum_tur_one_turbine)
#np.savetxt('hill1_sum_tur_two_turbine.txt',sum_tur_two_turbine)
#np.savetxt('hill1_x.txt',X[0,:])



### MKE ###

#grady_U2_one_turbine = np.gradient(hill_1_one_turbine_U**2,dy,axis=0)
#grady_U2_two_turbine = np.gradient(hill_1_two_turbine_U**2,dy,axis=0)
#
#adv_e_one = (0.5*hill_1_one_turbine_V*grady_U2_one_turbine)*((D/1000)/inlet_U**3)
#adv_e_two = (0.5*hill_1_two_turbine_V*grady_U2_two_turbine)*((D/1000)/inlet_U**3)
#
#tur_e_one = (-grady_U_one_turbine*hill_1_one_turbine_uv)*((D/1000)/inlet_U**3)
#tur_e_two = (-grady_U_two_turbine*hill_1_two_turbine_uv)*((D/1000)/inlet_U**3)
#
#
#y_ground = np.interp(X[0,:],x_,y_)
#y_aux = np.arange(0.5,1.505,0.005)
#
#sum_adv_e_one = np.zeros_like(X[0,:])
#sum_adv_e_two = np.zeros_like(X[0,:])
#sum_tur_e_one = np.zeros_like(X[0,:])
#sum_tur_e_two = np.zeros_like(X[0,:])
#
#for i in range (len(X[0,:])):
#    f = interpolate.interp1d(Y[:,i],adv_e_one[:,i])
#    sum_adv_e_one[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],adv_e_two[:,i])
#    sum_adv_e_two[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],tur_e_one[:,i])
#    sum_tur_e_one[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#    f = interpolate.interp1d(Y[:,i],tur_e_two[:,i])
#    sum_tur_e_two[i] = np.trapz(f(y_ground[i]+y_aux),dx=0.005)
#
#
#np.savetxt('hill1_sum_adv_e_one.txt',sum_adv_e_one)
#np.savetxt('hill1_sum_adv_e_two.txt',sum_adv_e_two)
#np.savetxt('hill1_sum_tur_e_one.txt',sum_tur_e_one)
#np.savetxt('hill1_sum_tur_e_two.txt',sum_tur_e_two)
#np.savetxt('hill1_x.txt',X[0,:])

### MKE ###

##### Interpolated data at x/D = (5,6,7,7.5,8,9,10)###
#
#x_d = [5,6,7,7.5,8,9,10]
#
#X_d = np.tile(x_d,(len(X[:,0]),1))
#
#y_ground = np.interp(x_d,x_,y_)
#
#U_base_profile = np.zeros([len(X[:,0]),len(x_d)])
#U_one_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#U_two_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#uv_one_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#uv_two_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#y_profile = Y[:,0]
#
#
#adv_one_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#adv_two_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#tur_one_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#tur_two_turbine_profile = np.zeros([len(X[:,0]),len(x_d)])
#
##
#for i in range(len(X[:,0])):
#   f = interpolate.interp1d(X[i,:],inlet_U[i,:])
#   U_base_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],hill_1_one_turbine_U[i,:])
#   U_one_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],hill_1_two_turbine_U[i,:])
#   U_two_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],hill_1_one_turbine_uv[i,:])
#   uv_one_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],hill_1_two_turbine_uv[i,:])
#   uv_two_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],adv_one_turbine[i,:])
#   adv_one_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],adv_two_turbine[i,:])
#   adv_two_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],tur_one_turbine[i,:])
#   tur_one_turbine_profile[i,:] = f(X_d[i,:])
#   f = interpolate.interp1d(X[i,:],tur_two_turbine[i,:])
#   tur_two_turbine_profile[i,:] = f(X_d[i,:])
#
#
#np.savetxt('hill1_U_base_profile.txt',U_base_profile)
#np.savetxt('hill1_U_one_turbine_profile.txt',U_one_turbine_profile)
#np.savetxt('hill1_U_two_turbine_profile.txt',U_two_turbine_profile)
#np.savetxt('hill1_uv_one_turbine_profile.txt',uv_one_turbine_profile)
#np.savetxt('hill1_uv_two_turbine_profile.txt',uv_two_turbine_profile)
#np.savetxt('hill1_y_profile.txt',y_profile)
#np.savetxt('hill1_y_ground.txt',y_ground)
#np.savetxt('hill1_x_ground.txt',x_d)
#
#np.savetxt('hill1_adv_one_turbine_profile.txt',adv_one_turbine_profile)
#np.savetxt('hill1_adv_two_turbine_profile.txt',adv_two_turbine_profile)
#np.savetxt('hill1_tur_one_turbine_profile.txt',tur_one_turbine_profile)
#np.savetxt('hill1_tur_two_turbine_profile.txt',tur_two_turbine_profile)
#
#y_ground_ = np.interp(X[0,:],x_,y_)
#y_aux = np.arange(0.4,1.6,0.01)
#deltaU_1tur = np.zeros_like(y_aux)
#w_1tur = np.zeros_like(X[0,:])
#deltaU_2tur = np.zeros_like(y_aux)
#w_2tur = np.zeros_like(X[0,:])
#
#
#for i in range(len(X[0,:])):
#   f = interpolate.interp1d(Y[:,i],deltaU_one_turbine[:,i])
#   w_1tur[i] = y_aux[np.argmax(np.abs(f(y_aux+y_ground_[i])))]
#   f = interpolate.interp1d(Y[:,i],deltaU_two_turbine[:,i])
#   w_2tur[i] = y_aux[np.argmax(np.abs(f(y_aux+y_ground_[i])))]
 

### save data ###
#
#with open('X_hill1','wb') as file:
#    pickle.dump(X,file)
#with open('Y_hill1','wb') as file:
#    pickle.dump(Y,file)
#
#### hill shape ###    
#with open('xhill1','wb') as file:
#    pickle.dump(x_,file)
#with open('yhill1','wb') as file:
#    pickle.dump(y_,file)
#
#with open('U_one_hill1','wb') as file:
#    pickle.dump(hill_1_one_turbine_U,file)
#with open('U_two_hill1','wb') as file:
#    pickle.dump(hill_1_two_turbine_U,file)
#    
#with open('V_one_hill1','wb') as file:
#    pickle.dump(hill_1_one_turbine_V,file)    
#with open('V_two_hill1','wb') as file:
#    pickle.dump(hill_1_two_turbine_V,file)
#
#with open('deltaU_one_hill1','wb') as file:
#    pickle.dump(deltaU_one_turbine,file)
#with open('deltaU_two_hill1','wb') as file:
#    pickle.dump(deltaU_two_turbine,file)
#    
#with open('adv_one_hill1','wb') as file:
#    pickle.dump(adv_one_turbine,file)
#with open('adv_two_hill1','wb') as file:
#    pickle.dump(adv_two_turbine,file)
#
#with open('tur_one_hill1','wb') as file:
#    pickle.dump(tur_one_turbine,file)
#with open('tur_two_hill1','wb') as file:
#    pickle.dump(tur_two_turbine,file)
#
#with open('turflux_one_hill1','wb') as file:
#    pickle.dump(turflux_one_turbine,file)
#with open('turflux_two_hill1','wb') as file:
#    pickle.dump(turflux_two_turbine,file)
#
#
#with open('advflux_one_hill1','wb') as file:
#    pickle.dump(advflux_one_turbine,file)
#with open('advflux_two_hill1','wb') as file:
#    pickle.dump(advflux_two_turbine,file)

####PLOTS
#
### Velocity Deficit one turbine

#plt.streamplot(X,Y,hill_1_one_turbine_U,hill_1_one_turbine_V,density=1,
#               linewidth=0.6,color='k',cmap='None',arrowsize=1,arrowstyle='->')
#plt.contourf(X,Y,np.abs(deltaU_one_turbine),cmap='plasma',levels=np.linspace(0,0.2,50))
#plt.plot(x_,y_,color='k')
##plt.plot(X_,Y_,color='r')
#plt.plot(x_,y_+0.5,color='w',linestyle='--')
#plt.plot(x_,y_+1.5,color='w',linestyle='--')
##plt.plot(X[0,:],w_1tur+y_ground_,color='w',linestyle=':')
#plt.axis('scaled')
#plt.axis([5,10,0,3.5])
#cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(0,0.2,5),
#format=OOMFormatter(-1, mathText=True)).set_label(r'$|U-U_{nw}| \times U_{nw}^{-1}$')
##plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0,3.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
####
#plt.savefig('deltaU_hill1_1T.png',dpi=300)
#plt.savefig('deltaU_hill1_1T.eps',dpi=300)
#plt.savefig('deltaU_hill1_1T.svg',dpi=300)
#
### Velocity Deficit two turbine
#
#plt.streamplot(X,Y,hill_1_two_turbine_U,hill_1_two_turbine_V,density=1,
#               linewidth=0.6,color='k',cmap='None',arrowsize=1,arrowstyle='->')
#plt.contourf(X,Y,deltaU_two_turbine,cmap='jet',levels=np.linspace(0,0.8,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='w',linestyle='--')
#plt.plot(x_,y_+1.5,color='w',linestyle='--')
##plt.plot(X[0,:],w_2tur+y_ground_,color='w',linestyle=':')
#plt.axis('scaled')
#plt.axis([6,10,0,3.5])
#cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(0,0.8,5),
#format=OOMFormatter(-1, mathText=True)).set_label(r'$|U-U_{nw}| \times U_{nw}^{-1}$')
#plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
#plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
#plt.savefig('deltaU_hill1_2T.png',dpi=300)
#plt.savefig('deltaU_hill1_2T.eps',dpi=300)
#plt.savefig('deltaU_hill1_2T.svg',dpi=300)

#plt.savefig('deltaU_hill1_2T_colorbar.png',dpi=300)
#plt.savefig('deltaU_hill1_2T_colorbar.eps',dpi=300)
#plt.savefig('deltaU_hill1_2T_colorbar.svg',dpi=300)

#
#
#

### Velocity one turbine

#plt.contourf(X,Y,hill_1_one_turbine_U/U_hub,cmap='jet',levels=np.linspace(0.0,1.5,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.axis([4.9,10.1,0.0,2.25])
#cb=plt.colorbar(orientation='vertical',ticks=np.linspace(0.0,1.5,7),
#format=OOMFormatter(0, mathText=True)).set_label(r'$U \times U_{H}^{-1}$')
#plt.xlabel('$x/D$')
#plt.yticks([1.0,2.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
#plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
##
#plt.savefig('U_hill1_1T.png',dpi=300)
#plt.savefig('U_hill1_1T.eps',dpi=300)
#plt.savefig('U_hill1_1T.svg',dpi=300)

#

### Velocity two turbine

#plt.contourf(X,Y,hill_1_two_turbine_U/U_hub,cmap='jet',levels=np.linspace(0.0,1.5,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.axis([5.9,10.1,0.0,2.25])
#cb=plt.colorbar(orientation='vertical',ticks=np.linspace(0.0,1.5,6),
#format=OOMFormatter(0, mathText=True)).set_label(r'$U \times U_{H}^{-1}$')
#plt.xlabel('$x/D$')
#plt.yticks([1.0,2.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
#plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
##
#plt.savefig('U_hill1_2T.png',dpi=300)
#plt.savefig('U_hill1_2T.eps',dpi=300)
#plt.savefig('U_hill1_2T.svg',dpi=300)


#### Velocity base
##
#plt.contourf(X,Y,inlet_U,cmap='jet',levels=np.linspace(0.0,12,50))
#plt.plot(x_,y_,color='k')
#plt.axis('scaled')
#plt.axis([4.9,10.1,0.0,2.25])
#cb=plt.colorbar(orientation='vertical',ticks=np.linspace(0.0,1.5,7),
#format=OOMFormatter(0, mathText=True)).set_label(r'$U \times U_{H}^{-1}$')
#plt.xlabel('$x/D$')
#plt.yticks([1.0,2.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
#plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
##
#plt.savefig('U_hill1_base.png',dpi=300)
#plt.savefig('U_hill1_base.eps',dpi=300)
#plt.savefig('U_hill1_base.svg',dpi=300)



### Kinetic Energy Flux - U <u'v'> One Turbine

#plt.contourf(X,Y,((-hill_1_one_turbine_U*hill_1_one_turbine_uv)/inlet_U**3),cmap='bwr',levels=np.linspace(-0.01,0.01,50))
#plt.axis('scaled')
#
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.axis([0,10.1,0.0,5.25])
##cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.02,0.02,5),
##format=OOMFormatter(-2, mathText=True)).set_label(r'$-U \overline{u^\prime v^\prime} \times U_{H}^{-3}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1.0,2.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
###
#plt.savefig('Uuv_hill_1_1T.png',dpi=300)
#plt.savefig('Uuv_hill_1_1T.eps',dpi=300)
#plt.savefig('Uuv_hill_1_1T.svg',dpi=300)


### Kinetic Energy Flux - U <u'v'> Two Turbine
#
#plt.contourf(X,Y,((-hill_1_two_turbine_U*hill_1_two_turbine_uv)/inlet_U**3),cmap='bwr',levels=np.linspace(-0.01,0.01,50))
#plt.axis('scaled')
#
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.axis([0,10.1,0.0,5.25])
##cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.02,0.02,5),
##format=OOMFormatter(-2, mathText=True)).set_label(r'$-U \overline{u^\prime v^\prime} \times U_{H}^{-3}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1.0,2.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
##
###
#plt.savefig('Uuv_hill_1_2T.png',dpi=300)
#plt.savefig('Uuv_hill_1_2T.eps',dpi=300)
#plt.savefig('Uuv_hill_1_2T.svg',dpi=300)


#### -V dUdy  ### One Turbine
#
### bl
#
#plt.contourf(X,Y,(-hill_1_one_turbine_V*grady_U_one_turbine)*((D/1000)/inlet_U**2),cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.axis('scaled')
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis([0,10.1,0,5.25])
#cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.15,0.15,7),
#format=OOMFormatter(-1, mathText=True)).set_label(r'$-V \partial{U} / \partial{y} \times D U_{nw}^{-2}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1,2.0,3.0,4.0,5.0])
#plt.xticks([5,6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
#plt.savefig('adv_hill1_one_turbine.png',dpi=300)
#plt.savefig('adv_hill1_one_turbine.eps',dpi=300)
#plt.savefig('adv_hill1_one_turbine.svg',dpi=300)
#

#### -V dUdy  ### Two turbine
#
### bl
###
#plt.contourf(X,Y,(-hill_1_two_turbine_V*grady_U_two_turbine)*((D/1000)/inlet_U**2),cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.axis('scaled')
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis([0,10.1,0,5.25])
##cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.15,0.15,7),
##format=OOMFormatter(-1, mathText=True)).set_label(r'$-V \partial{U} / \partial{y} \times D U_{base}^{-2}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1,2.0,3.0,4.0,5.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
#plt.savefig('adv_hill1_two_turbine.png',dpi=300)
#plt.savefig('adv_hill1_two_turbine.eps',dpi=300)
#plt.savefig('adv_hill1_two_turbine.svg',dpi=300)



#### duvdy  ### One Turbine
#
### bl

#plt.contourf(X,Y,(-grady_uv_one_turbine)*((D/1000)/inlet_U**2),cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.axis('scaled')
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis([0,10.1,0,5.25])
#cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.15,0.15,7),
#format=OOMFormatter(-1, mathText=True)).set_label(r'$-\partial{\overline{u^\prime v^\prime}} / \partial{y} \times D U_{nw}^{-2}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1,2.0,3.0,4.0,5.0])
#plt.xticks([5,6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
#plt.savefig('tur_hill_1_one_turbine.png',dpi=300)
#plt.savefig('tur_hill_1_one_turbine.eps',dpi=300)
#plt.savefig('tur_hill_1_one_turbine.svg',dpi=300)


#### -duvdy  ### Two turbine
#
#### bl
####
#plt.contourf(X,Y,(-grady_uv_two_turbine)*((D/1000)/inlet_U**2),cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.axis('scaled')
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis([0,10.1,0,5.25])
##cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.15,0.15,7),
##format=OOMFormatter(-1, mathText=True)).set_label(r'$-V \partial{U} / \partial{y} \times D U_{base}^{-2}$')
#
##plt.xlabel('$x/D$')
#plt.yticks([0,1,2.0,3.0,4.0,5.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
##plt.ylabel('$y/y_{H}$')
#plt.tight_layout()
#
#plt.savefig('tur_hill_1_two_turbine.png',dpi=300)
#plt.savefig('tur_hill_1_two_turbine.eps',dpi=300)
#plt.savefig('tur_hill_1_two_turbine.svg',dpi=300)



######## Streamwise Velocity deficit at 3D#######
#
#f=plt.figure()
#ax=f.add_subplot(111)
#plt.plot(deltaU_bl_3D/U_bl,Z_bl,label='Turbulent BL',linestyle=':',color='k')
#plt.plot(deltaU_jet_ps_3D/U_jet_ps,Z_jet_ps,label='Jet at Top Tip')
#plt.plot(deltaU_jet_mid_3D/U_jet_mid,Z_jet_mid,label='Jet at Midspan')
#plt.xlim(0,0.35)
#plt.ylim(-1,1)
#plt.axhline(y=0.5,linestyle='--',color='k')
#plt.axhline(y=-0.5,linestyle='--',color='k')
#plt.xlabel('$\Delta U/U_{H}$')
#plt.xticks([0.0,0.1,0.2,0.3])
#plt.yticks([-1,-0.5,0,0.5,1])
#plt.ylabel('$(y-H)/L$')
##plt.legend(prop={'size':14})
#plt.tight_layout()
#
#plt.savefig('deltaU_3D.png',dpi=300)
#plt.savefig('deltaU_3D.eps',dpi=300)
#plt.savefig('deltaU_3D.svg',dpi=300)


### Velocity deficit along the centerline###

#f=plt.figure()
#ax=f.add_subplot(111)
#plt.plot(X,deltaU_bl_hub/U_bl,label='Turbulent BL',color='k',linestyle=':')
#plt.plot(X,deltaU_jet_ps_hub/U_jet_ps,label='LLJ peak at top tip',color='C0')
#plt.plot(X,deltaU_jet_mid_hub/U_jet_mid,label='LLJ peak at midspan',color='C1')
#plt.xlabel(r'$ x/D $')
#plt.ylabel(r'$\Delta{U}/U_{H}$')
#plt.yticks([0.2,0.3,0.4])
#plt.legend(prop={'size':14})
#plt.tight_layout()
#
#plt.savefig('deltaU_hub.png',dpi=300)
#plt.savefig('deltaU_hub.eps',dpi=300)
#plt.savefig('deltaU_hub.svg',dpi=300)

### Advection one turbine ###

#plt.contourf(X,Y,adv_e_one,cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.xlim(5,10)
#plt.ylim(0,3.5)
##cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.15,0.15,7),
##format=OOMFormatter(-1, mathText=True)).set_label(r'$1/2 V \partial{U^2}/\partial{y} \times D U_{base}^{-3}$')
###plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0,3.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
###plt.ylabel('$y/y_{H}$')
##plt.tight_layout()
#######
#plt.savefig('adv_e_one.png',dpi=300)
#plt.savefig('adv_e_one.eps',dpi=300)
#plt.savefig('adv_e_one.svg',dpi=300)


### Advection two turbine ###
#
#plt.contourf(X,Y,adv_e_two,cmap='seismic',levels=np.linspace(-0.15,0.15,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.xlim(6,10)
#plt.ylim(0,3.5)
##cb=plt.colorbar(orientation='vertical',ticks=np.linspace(-0.15,0.15,7),
##format=OOMFormatter(-1, mathText=True)).set_label(r'$\frac{1}{2} \,  V  \, \frac{\partial{U^2}}{\partial{y}} \, \times \, \frac{D}{U_{local}^{-3}} $')
###plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0,3.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
###plt.ylabel('$y/y_{H}$')
##plt.tight_layout()
#######
##plt.savefig('adv_e_two.png',dpi=300)
##plt.savefig('adv_e_two.eps',dpi=300)
##plt.savefig('adv_e_two.svg',dpi=300)
#
#plt.savefig('adv_e_two_colorbar.png',dpi=300)
#plt.savefig('adv_e_two_colorbar.eps',dpi=300)
#plt.savefig('adv_e_two_colorbar.svg',dpi=300)

### Turbulence one turbine ###
#
#plt.contourf(X,Y,tur_e_one,cmap='seismic',levels=np.linspace(-0.015,0.015,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.xlim(5,10)
#plt.ylim(0,3.5)
###cb=plt.colorbar(orientation='horizontal',ticks=np.linspace(-0.25,0.25,7),
###format=OOMFormatter(-1, mathText=True)).set_label(r'$1/2 V \partial{U^2}/\partial{y} \times D U_{base}^{-3}$')
####plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0,3.0])
#plt.xticks([5.0,6.0,7.0,8.0,9.0,10.0])
####plt.ylabel('$y/y_{H}$')
###plt.tight_layout()
########
#plt.savefig('tur_e_one.png',dpi=300)
#plt.savefig('tur_e_one.eps',dpi=300)
#plt.savefig('tur_e_one.svg',dpi=300)

### Turbulence two turbine ###

#plt.contourf(X,Y,tur_e_two,cmap='seismic',levels=np.linspace(-0.015,0.015,50))
#plt.plot(x_,y_,color='k')
#plt.plot(x_,y_+0.5,color='k',linestyle='--')
#plt.plot(x_,y_+1.5,color='k',linestyle='--')
#plt.axis('scaled')
#plt.xlim(6,10)
#plt.ylim(0,3.5)
##cb=plt.colorbar(orientation='vertical',ticks=np.linspace(-0.015,0.015,7),
##format=OOMFormatter(-2, mathText=True)).set_label(r'$ U \, \frac{\partial{\,\overline{u^\prime v^\prime}}}{\partial{y}} \, \times \, \frac{D}{U_{local}^{-3}} $')
####plt.xlabel('$x/D$')
#plt.yticks([0.0,1.0,2.0,3.0])
#plt.xticks([6.0,7.0,8.0,9.0,10.0])
###plt.ylabel('$y/y_{H}$')
##plt.tight_layout()
#######
#plt.savefig('tur_e_two.png',dpi=300)
#plt.savefig('tur_e_two.eps',dpi=300)
#plt.savefig('tur_e_two.svg',dpi=300)
#
#plt.savefig('two_e_two_colorbar.png',dpi=300)
#plt.savefig('two_e_two_colorbar.eps',dpi=300)
#plt.savefig('two_e_two_colorbar.svg',dpi=300)

