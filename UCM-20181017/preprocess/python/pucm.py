###
import glob, os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



idir = '../kugahara/'
files = glob.glob(idir+'kugahara*.txt')

sd = pd.read_csv('../kugahara/kugahara_0901_sd.txt', delim_whitespace=True, header=None)



var = pd.read_csv('../kugahara/kugahara_0901_vars.txt', delim_whitespace=True, header=None)


var.columns = ['sd','lw','temp','qv','u','ah']

for c in  var.columns:
    var[c].plot()
    plt.show()




    
    
    
    
    



