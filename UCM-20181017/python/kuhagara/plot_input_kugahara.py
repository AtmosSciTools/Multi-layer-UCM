##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, glob
import xarray as xr
#import matplotlib.patches as Patch
#import matplotlib.collections as collections
from matplotlib.lines import Line2D

def set_fig (figsize=(4.5,3), axe = [0.1,0.15,0.8,0.75]):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axe)
    return fig, ax

utype = 0
f = open("../../namelist/URBPARM.TBL.kugahara")
lines = f.readlines()
hl = float([l for l in lines if "BUILDING_HEIGHT:" in l][0].split(':')[1].split(",")[utype])
wl = float([l for l in lines if "ROOF_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
rl = float([l for l in lines if "ROAD_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
x1 = np.array([wl**2, 2*wl*hl+rl**2] + 4*[wl*hl]) 
x2 = x1/(rl+wl)**2
frc_u = float([l for l in lines if "FRC_URB:" in l][0].split(':')[1].split(",")[utype]) 
print lab[i], hl, wl, rl, frc_u
ah = float([l for l in lines if "AH:" in l][0].split(':')[1].split(",")[utype])  
ahd = np.array([float(a) for a in [l for l in lines if "AH_DIURNAL_PROFILE:" in l][0].split(':')[1].split()])


cc = ['#1f77b4', '#d62728','#2ca02c',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']  
    
    
    
    
df = pd.read_csv('../../preprocess/20180718/data/ucm_data/kugahara_topvars_20050901.csv')
dp = pd.read_csv('../../preprocess/20180718/data/Kugahara-full-datasets/tmp_profile_20050901.csv', index_col=0, parse_dates=True)
dp.columns = [float(c.split('m')[0][4:]) for c in dp.columns]


fig = plt.figure(figsize=(7.5,5))
xl, yl = [0.11,0.57,0.11,.57], [0.58,0.58,.1,.1]
leg = [["Sdown", "Ldown"],["Temperature"],["Wind speed"],["Specific humidity"]]
ylab = ["($W/m^2$)","($^oC$)","($m/s$)","($W/m^2$)"]
ylim = [(-50,900),(20,33),(0,5),(10,20)]
tit = ["a)", "b)", "c)", "d)"]
for i in range(4):
    ax = plt.axes((xl[i],yl[i],0.35,.38))
    if i == 0:
        df[["swtop"]][:24].plot(ax=ax, color=cc[1],lw=1, ls ="-")
        df[["lwtop"]][:24].plot(ax=ax, color=cc[1],lw=1, ls ="--")
        ax.set_ylim(-50,900)
        ax.legend(leg[i] )
    if i == 1: 
        df[["ttop"]][:24].plot(ax=ax, color=cc[0],lw=1)
        ax.set_ylim(20,33)
        ax.legend(["Temperature"], loc=1)
        ax1 = ax.twinx()
        df[["qtop"]][:24].plot(ax=ax1,color="g",lw=1)
        ax1.set_ylim(10,20)
        ax1.set_ylabel("($g/kg$)")
        ax1.legend(["Specific humidity"], loc=4)
    if i == 2: 
        df[["wstop"]][:24].plot(ax=ax,color=cc[4],lw=1)
        ax.set_ylim(0,5)
        ax.legend(leg[i] )
    if i == 3: 
        ax.plot(range(24),ahd*ah,color="k",lw=1)
        ax.set_ylim(0,27)
        ax.legend(["AH"] )
    #ax.legend(leg[i] )
    ax.set_ylabel(ylab[i])
    #ax.set_ylim( ylim[i] )
    ax.set_xlabel("LST ($hour$)")
    ax.text(0.0,1.02,tit[i], transform = ax.transAxes, fontsize=15)

fig.savefig("input_data.png", dpi=300)





    
    
     
 