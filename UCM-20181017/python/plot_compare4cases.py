##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, glob
import xarray as xr



def set_fig (figsize=(5,3.5), axe = [0.15,0.2,0.8,0.75]):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(axe)
    return fig, ax

o = pd.read_csv("tokyo_amedas.csv", index_col=0, header=None)

dd = []
fig, ax = set_fig ()


lab = ["u case 1", "u case 2", "u case 3", "slab"]
c = ["steelblue", "orange", "r", "g"]

dat = []
dd = []
for i, ifi in enumerate( ["case1", "case2", "case3","slab3" ]): 
    ifile = ifi+'_mucmout_2014_08_01.nc'

    ds = xr.open_dataset(ifile)
    t2 = pd.to_datetime(ds.time.values[-1]).strftime('%Y-%m-%d %H')
    t1 = pd.to_datetime(ds.time.values[-1] - np.timedelta64(24, 'h')).strftime('%Y-%m-%d %H')
    d = ds.sel(time=slice(t1,t2))
    dd.append(d)
    dat.append( pd.Series( d["T2"][6::6].values,
                           index = pd.to_datetime(d["T2"][6::6].time.values).round("1min"), name = lab[i]))

    d["T2"][6::6].plot(ax=ax,label=lab[i], c = c[i])
    
ax.plot(d["T2"][6::6].time, o[1][1:].values,c = 'none', 
                  marker = 'o', 
                  markersize = 7, 
                  markeredgecolor = 'k',
                  alpha = 1.,
                  markeredgewidth = 1.,
                  markerfacecolor = 'None',
                  label = 'OBS')


dat.append(pd.Series(o[1][1:].values,index=pd.to_datetime(d["T2"][6::6].time.values).round("1min"), name = "obs"))
ax.set_ylabel("Temperature at 2 m ($C^oC$)")
ax.set_xlabel("2014")
plt.legend()
ax.set_ylim(22,40)
fig.savefig("fig/T2_comparison.png", dpi=300)
pd.concat(dat,axis=1).to_csv("fig/T2_time_series.csv")




for i, d in enumerate( dd ): 

    d1 = d    

    
    fig, ax = set_fig ()
    tld = pd.to_datetime(d.time.values)
    tldr = tld.round("1min")
    z_dim = d.z_dim[1:]
    
    for h in [0, 6, 12, 18]:
        ind = np.argwhere((tldr.hour == h) & (tldr.minute == 0))[0,0]
        t = d.sel(time=tld[ind])["TEMP"][1:]
        ax.plot(t.values, z_dim, label = str(h)+" h")
    
    
    ax.legend()
    ax.set_xlabel("Potential temperature ($^oC$)")
    ax.set_ylabel('Height ($m$)')
    ax.set_xlim(24,40)
    ax.text(0.,1.02,lab[i], transform = ax.transAxes, fontsize=13)
    fig.savefig("fig/Pot_temp_"+lab[i].replace(' ','_')+".png", dpi=300)
    
    
    #ax.set_ylim(0,50)
    
    
    if "RNET_UCM" in d.variables.keys():
        
        utype = i 
        f = open("../URBPARM.TBL")
        lines = f.readlines()
        hl = float([l for l in lines if "BUILDING_HEIGHT:" in l][0].split(':')[1].split(",")[utype])
        wl = float([l for l in lines if "ROOF_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
        rl = float([l for l in lines if "ROAD_WIDTH:" in l][0].split(':')[1].split(",")[utype])  
        x1 = np.array([wl**2, 2*wl*hl+rl**2] + 4*[wl*hl]) 
        x2 = x1/(rl+wl)**2
        frc_u = float([l for l in lines if "FRC_URB:" in l][0].split(':')[1].split(",")[utype]) 
        print lab[i], hl, wl, rl, frc_u
            
        
        hbus = ["RNET_UCM","H_UCM","GR_UCM","LE_UCM"]
        hbvs = ["RNET_SLB","H_SLB","GR_SLB","LE_SLB"]
        
        fig, ax = set_fig ()
        for v1, v2 in zip(hbus,hbvs):
            dd = (d[v1]*x2).sum(dim="surface")*frc_u + d[v2]*(1.-frc_u)
            dd[::6].plot(ax=ax,label=v1)
        ax.set_ylabel(d[v1].units)
        plt.legend()
        ax.set_ylabel('Heat flux ($Wm^{-2}$)')
        ax.text(0.,1.02,lab[i], transform = ax.transAxes, fontsize=13)
        ax.set_ylim(-200,800)
        fig.savefig("fig/Head_budget_"+lab[i].replace(' ','_')+".png", dpi=300)
    
    else:
        
        hbvs = ["RNET_SLB","H_SLB","GR_SLB","LE_SLB"]
        fig, ax = set_fig ()
        for var in hbvs:
            d1[var].plot(ax=ax,label=var)
            ax.set_ylabel(d[var].units)
            plt.legend()
        ax.set_ylabel('Heat flux ($Wm^{-2}$)')
        plt.legend()
        ax.set_ylim(-200,800)
        ax.text(0.,1.02,lab[i], transform = ax.transAxes, fontsize=13)
        fig.savefig("fig/Head_budget_"+lab[i].replace(' ','_')+".png", dpi=300)
    
















     
 