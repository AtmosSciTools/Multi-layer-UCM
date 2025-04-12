####
####
import pandas as pd
import sys, glob, os



odir = 'ucm_data/ics/'
if not os.path.exists(odir): os.makedirs(odir)

do = pd.DataFrame([-0.75, 0.75, 3., 5., 7., 9., 11., 14., 17., 21., 25., 28.])
ofile = odir + 'grid.dat'
do.to_csv(ofile, index=None, header=None)
    
    







            





