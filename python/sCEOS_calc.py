#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:01:04 2021

@author: kanferg
"""
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageEnhance
#from PIL import fromarray
from skimage import data,io
from skimage.filters import threshold_otsu, threshold_local, threshold_local
from skimage.morphology import convex_hull_image
import matplotlib.pyplot as plt
import scipy as ndimage
from scipy.ndimage.morphology import binary_opening
from skimage.morphology import disk
from skimage.segmentation import watershed
from skimage import data
from skimage.filters import rank
from skimage.util import img_as_ubyte
from skimage import data, util
from skimage.measure import label
from skimage.measure import perimeter
from skimage import measure
from skimage import exposure, io, util
import os
import pandas as pd
from pandas import DataFrame
from scipy.ndimage import label, generate_binary_structure
from scipy.ndimage.morphology import binary_fill_holes
import cv2


os.chdir('/Users/kanferg/Desktop/Gil_LabWork/ANNA-PALM/Simulator/Similtor_genration/Noise_simulator')
#ret,orig_100=cv2.imreadmulti('100.tif', [], cv2.IMREAD_ANYCOLOR)
ret,orig_100=cv2.imreadmulti('100.tif', [], cv2.IMREAD_ANYDEPTH)
np.shape(orig_100)
print(orig_100)
temp_flat=orig_100[1].ravel()
print(temp_flat)
np.min(temp_flat)




###################
# Calculate offset
##################
#1) Photon per pixel = (electrons (grey-value / conversion coefficent 0.49 for the flash 4)/0.83 (QE)))   
#2) Oi = summary of photon per pixel accrose frame / Total number of frame 

# First turn all to large data frame using panda
M=3000
h=64
l=64
CC=0.49
QE=0.83
data = pd.DataFrame()
for i in range(M):
#    temp_flat=orig_100[0].ravel()
#    temp_flat=temp_flat/256
#    temp_flat_e=(temp_flat/0.49)/0.83
    df = pd.DataFrame(data=orig_100[i].flatten())
    df = df.T
    data=data.append(df)
    print(i)
data.reset_index(drop=True,inplace=True)
# float all the values
data=data/256
# getting electrons
electron=(data/CC)/QE
Offset_df=electron.sum(axis=0)
Offset_df=Offset_df/3000
Offset_df=Offset_df.T

# convert df to np array 
Offset_mat_1d=Offset_df.values
offset_mat_2d=np.reshape(Offset_mat_1d,(64,64))

#####################
# Calculate Variance
#####################
# calculate for 60,000 frame (in this example I have only 3000)
#electron**2 - offset**2
electro_sq=electron**2
Offset_mat_1d_sq=Offset_mat_1d**2

Var_df = pd.DataFrame()
for i in range(M):
    temp=electro_sq.loc[i]-Offset_mat_1d_sq   # can it be minus ? chack
    temp=temp.T
    Var_df=Var_df.append(temp)
    print(i)
Var_df_sum=Var_df.sum(axis=0)
Var_df_sum=Var_df_sum/3000
Var_df_sum_1d=Var_df_sum.values
Var_df_sum_2d=np.reshape(Var_df_sum_1d,(64,64))

#####################
# Calculate Gain
#####################
#K total number illumination (21 to 200 photon per pixel) we have 5 illumination
#vi in k variance of illumination k minus the previus veriation at K
#Di in K is the mean of the electron count in sequance k minus oi that was calculated previsuly




