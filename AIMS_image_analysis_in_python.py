


import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image
from numpy import asarray
from skimage import data


import os


img=io.imread('wella1A1_0004.tif')

########normalisation

#frame<-readImage("C:/outproc/color.tif",type = "tif")
#GFP<-frame[1:512,1:512,1]
#dapi<-frame[1:512,1:512,2]
## minmax normalization
#minDapi<-min(as.vector(dapi))
#maxDapi<-max(as.vector(dapi))
#minGFP<-min(as.vector(GFP))
#maxGFP<-max(as.vector(GFP))
#dapin<-normalize(dapi, ft=c(0,1),c(minDapi,maxDapi))
#GFPn<-normalize(GFP, ft=c(0,1) ,c(minGFP,maxGFP))

os.chdir('/Users/gilkanfer/iCloud Drive (Archive)/Desktop/gil/Gil studies/NIH_PostDoc/Python/Computur_vision/Images/')
image = Image.open('wella1A1_0004.tif')
# summarize some details about the image
print(image.format)
print(image.mode)
print(image.size)
# show the image
image.show()
red, green, blue = image.split()
red.show()


data = image.getdata()


pixels = asarray(image)
