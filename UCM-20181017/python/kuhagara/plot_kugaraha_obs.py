##
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys, glob
import xarray as xr

def set_fig ():
    fig = plt.figure(figsize=(5,3.5))
    ax = fig.add_axes([0.1,0.15,0.8,0.75])
    return fig, ax

    
    
    
ifile = '../mucmout_2005_09_01.nc'


ds = xr.open_dataset(ifile)
df = pd.read_csv('../preprocess/20180718/data/ucm_data/kugahara_topvars_20050901.csv')
dp = pd.read_csv('../preprocess/20180718/data/Kugahara-full-datasets/tmp_profile_20050901.csv', index_col=0, parse_dates=True)
dp.columns = [float(c.split('m')[0][4:]) for c in dp.columns]
#t75 = dp["Tair0.75m (oC)"]   
#t75.plot()
dp[[0.75,11,28]].plot()



















sys.exit()

t2 = pd.to_datetime(ds.time.values[-1]).strftime('%Y-%m-%d %H')
t1 = pd.to_datetime(ds.time.values[-1] - np.timedelta64(24, 'h')).strftime('%Y-%m-%d %H')
d = ds.sel(time=slice(t1,t2))


fig, ax = set_fig ()
for var in ["RNET_SLB","H_SLB","GR_SLB","LE_SLB"]:
    d[var].plot(ax=ax,label=var)
ax.set_ylabel(d[var].units)
plt.legend()



for var in ["T2"]:
    fig, ax = set_fig ()
    d[var].plot(ax=ax,label=var)
    ax.set_ylabel(d[var].units)
    plt.legend()



for var in ["TEMP", "QV", "KM", "U"]:
    fig, ax = set_fig ()
    tld = pd.to_datetime(d.time.values)
    for h in [9, 15, 21]:
        print h
        ind = np.argwhere((tld.hour == h) & (tld.minute == 0))[0,0]
        print ind
        t = d.sel(time=tld[ind])[var][1:]
        ax.plot(t.values, t.z_dim, label = str(h)+" h")
    ax.legend()
    ax.set_xlabel(t.units)
    ax.set_ylabel('Height (m)')


    
    
    
    
    
    
    
    
sys.exit()








plt.show()
ds.T2.plot()
plt.show()
#x = pd.DataFrame(ds.RNET_UCM.values,index=ds.time,columns=['R','G','1','2','3','4'])
#x.plot()
ds.TAU_T_SLB.plot()
plt.show()
tt = ds.TEMP[ [ 479, 539, 599, 659 ], 1: ]

for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.TEMP[i,1:],ds.z_dim[1:])
plt.show()

for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.KH[i,1:],ds.z_dim[1:])
plt.show()

for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.KM[i,1:],ds.z_dim[1:])
plt.show()

for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.U[i,1:],ds.z_dim[1:])
plt.show()
for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.V[i,1:],ds.z_dim[1:])
plt.show()

for i in [ 479, 539, 599, 659 ]:
    plt.plot(ds.QV[i,1:],ds.z_dim[1:])
plt.show()



x = pd.DataFrame(ds.TS_UCM.values,index=ds.time,columns=['R','G','1','2','3','4'])
x.plot()


sys.exit()

var = [u'z_dim',
 u'time',
 u'TEMP',
 u'QV',
 u'U',
 u'V',
 u'KM',
 u'KH',
 u'T2',
 u'Q2',
 u'U10',
 u'V10',
 u'ELEVATION',
 u'ADI',
 u'COSADI',
 u'SINADI',
 u'BRiNu',
 u'CM_SLB',
 u'CH_SLB',
 u'TAU_U_SLB',
 u'TAU_V_SLB',
 u'TAU_T_SLB',
 u'RNET_SLB',
 u'H_SLB',
 u'LE_SLB',
 u'GR_SLB',
 u'TS_SLB',
 u'THETA_T',
 u'QV_T',
 u'U_T',
 u'V_T',
 u'SW_B',
 u'LW_B',
 u'COSZ',
 u'OMG',
 u'DECLIN',
 u'AZTH_C',
 u'AZTH_S',
 u'SOL_ELEV']

x = [u'z_dim',
     u'time',
     u'TEMP',
     u'QV',
     u'U',
     u'V',
     u'KM',
     u'KH',
     u'T2',
     u'Q2',
     u'U10',
     u'V10',
     u'ELEVATION',
     u'ADI',
     u'COSADI',
     u'SINADI',
     u'BRiNu',
     u'CM_SLB',
     u'CH_SLB',
     u'TAU_U_SLB',
     u'TAU_V_SLB',
     u'TAU_T_SLB',
     u'RNET_SLB',
     u'H_SLB',
     u'LE_SLB',
     u'GR_SLB',
     u'TS_SLB',
     u'RNET_UCM',
     u'H_UCM',
     u'LE_UCM',
     u'GR_UCM',
     u'TS_UCM',
     u'CM_UCM',
     u'CH_UCM',
     u'SD_UCM',
     u'Ri_UCM',
     u'SHADE_UCM',
     u'SWNET',
     u'LWNET',
     u'LWDOWN',
     u'LWUP',
     u'TAU_U_UCM',
     u'TAU_V_UCM',
     u'TAU_T_UCM',
     u'THETA_T',
     u'QV_T',
     u'U_T',
     u'V_T',
     u'SW_B',
     u'LW_B',
     u'COSZ',
     u'OMG',
     u'DECLIN',
     u'AZTH_C',
     u'AZTH_S',
     u'SOL_ELEV']
     
 