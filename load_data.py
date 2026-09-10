##Load Airfoil Data and Blade Data to use in next steps

import numpy as np
import pandas as pd

#Load blade data and airfoil data from txt files
blade_dat = pd.read_fwf(r'data/bladedat.txt',header = None)
blade_dat.columns =['r','c','beta','t/c']


airfoil_data = np.zeros((6, 105,4))  # 6 airfoils,, 105 data points each,4 columns (alpha, C_l(alpha), C_d(alpha), C_m(alpha))
#each airfoil has one thickness/chord ratio t/c 

# airfoil_data = pd.read_fwf(r'data/cylinder.txt',header = None)
airfoil_data[0] = pd.read_csv(r'data/cylinder.txt',header = None,sep = None)
airfoil_data[1]= pd.read_csv(r'data/FFA-W3-600.txt',header = None,sep = None)
airfoil_data[2]= pd.read_csv(r'data/FFA-W3-480.txt',header = None,sep = None)
airfoil_data[3]= pd.read_csv(r'data/FFA-W3-360.txt',header = None,sep = None)
airfoil_data[4]= pd.read_csv(r'data/FFA-W3-301.txt',header = None,sep = None)
airfoil_data[5]= pd.read_csv(r'data/FFA-W3-241.txt',header = None,sep = None)


print(airfoil_data[3])



