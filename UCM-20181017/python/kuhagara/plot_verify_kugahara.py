##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, glob
import xarray as xr
#import matplotlib.patches as Patch
#import matplotlib.collections as collections
from matplotlib.lines import Line2D

def set_fig (figsize=(5,3.3), axe = [0.15,0.15,0.8,0.75]):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axe)
    return fig, ax

new_colors = ['#1f77b4','#2ca02c', '#d62728','#ff7f0e']
    
ifile = 'mucmout_2005_09_01.nc'
ifile = '../../mucmout_2005_09_01.nc'



ds = xr.open_dataset(ifile)
df = pd.read_csv('../../preprocess/20180718/data/ucm_data/kugahara_topvars_20050901.csv')
dp = pd.read_csv('../../preprocess/20180718/data/Kugahara-full-datasets/tmp_profile_20050901.csv', index_col=0, parse_dates=True)
dp.columns = [float(c.split('m')[0][4:]) for c in dp.columns]


t2 = pd.to_datetime(ds.time.values[-1]).strftime('%Y-%m-%d %H')
t1 = pd.to_datetime(ds.time.values[-1] - np.timedelta64(24, 'h')).strftime('%Y-%m-%d %H')
d = ds.sel(time=slice(t1,t2))

do = pd.DataFrame({'obs':dp[0.75].values,'mucm':d.TEMP[:-1:6,1].values}, index=dp.index)

#==============================
#
#
fig, ax = set_fig()
ax.text(0.0,1.02,"", transform = ax.transAxes, fontsize=15)
ax.plot(range(24),do.iloc[:,0].values, c = 'k', label = 'MUCM',lw=1.2)
ax.plot(range(24),do.iloc[:,1].values,
                  c = 'none', 
                  marker = 'o', 
                  markersize = 6, 
                  markeredgecolor = 'k', 
                  markeredgewidth = 0.7,
                  markerfacecolor = 'k',
                  #alpha=0.6,
                  label = 'OBS')
plt.grid(True,ls="--")
plt.legend(numpoints = 1, fontsize=12)
ax.set_ylabel('($^oC$)')
ax.set_ylim(21,35)
ax.set_xlabel("LST ($hour$)")
fig.savefig("validate_t0p75m.png", dpi=300)

#==============================
#
#



fig = plt.figure(figsize=(7,4))
xlim = [[20,36],[0,4],[10,20]]

text = ["temperature", "wind speed", "q vapor"]
tit = ["a)","b)"]
ylab = ["Temperature ($^oC$)","Wind speed ($m/s$)"]
for ii, (var, varo) in enumerate(zip(["TEMP", "U", "QV"], ['tmp','wps','hum'])[:2]):
    dp = pd.read_csv('../../preprocess/20180718/data/Kugahara-full-datasets/'+varo+'_profile_20050901.csv', index_col=0, parse_dates=True)
    if varo == 'tmp': dp.columns = [float(c.split('m')[0][4:]) for c in dp.columns]
    if varo == 'hum': dp.columns = [float(c.split('m')[0][3:]) for c in dp.columns]
    if varo == 'wps': dp.columns = [float(c.split('m')[0][1:]) for c in dp.columns]

                 
    ax = plt.axes((0.1+ii*0.4,0.15,0.3,.8))
    
    tld = pd.to_datetime(d.time.values)
    tldr = tld.round("1min")
    cc = ['b', 'g', 'r', 'orange']
    leg1 = []
    leg2 = []
    z = ds.z_dim[1:]
    for i, h in enumerate([0,6, 12, 18][:]):
        print ii, h
        ind = np.argwhere((tldr.hour == h) & (tldr.minute == 0))[0,0]
        t  = d.sel(time=tld[ind])[var][1:]
        
        if var == "QV": tu = t.values * 0.001
        if var == "TEMP": tu = t.values-t.z_dim.values*0.004
        if var == "U":
            tu = (d.sel(time=tld[ind])["U"][1:]**2 + d.sel(time=tld[ind])["V"][1:]**2)**0.5
        
        to = dp[dp.index.hour == h].values[0]
        ax.plot(tu, z, label = str(h)+" h", c = new_colors[i],lw=1.2)
        ax.scatter(to, dp.columns, label = str(h)+" h", c = 'none', s = 40,
                  marker = 'o', edgecolors = new_colors[i])
        leg1.append(Line2D( [0],[0],color=new_colors[i]))
        leg2.append(Line2D( [0],[0],color="none", marker = 'o', markeredgecolor=new_colors[i], label = str(h)+' h' ))
    ax.set_ylim(0,30)
    ax.set_xlim(xlim[ii])
    if ii == 1:
        l = ax.legend(handles = leg1 + leg2, ncol=2, columnspacing=0, numpoints = 1, title = 'MUCM OBS',
                  fontsize=10, frameon=True, bbox_to_anchor=[0.1,0.17,1.55,0.2])
        l.get_title().set_fontsize(9)
        l._legend_box.align = 'left'
    ax.set_xlabel(ylab[ii])
    ax.set_ylabel('Height ($m$)')
    ax.text(0,1.01,tit[ii],transform = ax.transAxes, fontsize=15)

fig.savefig("validate_vert.png", dpi=300)





utype = 0
f = open("../../namelist/URBPARM.TBL.kugahara")
lines = f.readlines()
hl = float([l for l in lines if "BUILDING_HEIGHT:" in l][0].split(':')[1].split(",")[utype])
wl = float([l for l in lines if "ROOF_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
rl = float([l for l in lines if "ROAD_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
x1 = np.array([wl**2, 2*wl*hl+rl**2] + 4*[wl*hl]) 
x2 = x1/(rl+wl)**2
frc_u = float([l for l in lines if "FRC_URB:" in l][0].split(':')[1].split(",")[utype]) 
print  hl, wl, rl, frc_u
    

hbus = ["RNET_UCM","H_UCM","GR_UCM","LE_UCM"]
hbvs = ["RNET_SLB","H_SLB","GR_SLB","LE_SLB"]

fig, ax = set_fig ()

new_colors = ['#1f77b4', '#d62728','#ff7f0e', '#2ca02c']       
mk = ["x", "o", "v", ""]
for ii, (v1, v2) in enumerate(zip(hbus,hbvs)):
    dd = (d[v1]*x2).sum(dim="surface")*frc_u + d[v2]*(1.-frc_u)
    ax.plot(range(24),dd[:-1:6].values,lw=1.2,marker=mk[ii],ms=4, c = new_colors[ii])
ax.set_ylabel(d[v1].units)
ax.legend(["Rnet", "H", "Gr", "Le"])
ax.set_ylabel('Heat flux ($W/m^{2}$)')
#ax.text(0.,1.02,lab[i], transform = ax.transAxes, fontsize=13)
ax.set_ylim(-110,750)
ax.set_xlim(-0.5,23.5)
ax.set_xlabel("LST ($hour$)")
ax.axhline(linewidth=1, color='k',ls="--")
fig.savefig("Head_budget.png", dpi=300)








    
    
    
    
    
    
     
 