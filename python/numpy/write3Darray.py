fin_zoom=fin[13:16,17:20,2]
os.chdir('/Users/kanferg/Desktop/Gil_LabWork/ANNA-PALM/Simulator/Similtor_genration/Data')
path="/Users/kanferg/Desktop/Gil_LabWork/ANNA-PALM/Simulator/Similtor_genration/Data/temp"

import imageio

#os.mkdir(path) 
fin_temp=fin[0:32,0:32,0:5]
for i in range(6):
    newimage = fin_temp[:, :, i]
    imageio.imwrite("/Users/kanferg/Desktop/Gil_LabWork/ANNA-PALM/Simulator/Similtor_genration/Data/temp/image%d.tiff" %i, newimage)

