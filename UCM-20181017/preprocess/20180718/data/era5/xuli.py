#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 09:55:23 2018

@author: doan
"""

import xarray as xr
import sys


ifile = "/Users/doan/Downloads/adaptor.mars.internal-1545788751.2480798-20108-5-e2fffe7a-e902-40cc-9be0-06a523bb384b.nc"
ifile0 = "/Users/doan/Downloads/adaptor.mars.internal-1545793459.27678-5302-13-e9a41bca-8379-4778-93ab-893b2013a306.nc"

ds = xr.open_mfdataset([ifile0, ifile])
#xr.open_dataset(ifile)
jj = np.argwhere(((ds.latitude > 30) & (ds.latitude < 40)).values)[:,0]
ii = np.argwhere(((ds.longitude > 135) & (ds.longitude < 145)).values)[:,0]
d = ds.isel(latitude=jj,longitude=ii)
d.coords["time"] = pd.to_datetime(d.coords["time"].values) + pd.Timedelta("9 hours")



import matplotlib.pyplot as plt
import numpy as np

import cartopy
import cartopy.crs as ccrs

def divergence(F):
    """ compute the divergence of n-D scalar field `F` """
    return reduce(np.add,np.gradient(F))


shape=(20, 30)
crs = ccrs.NorthPolarStereo()
#scale = 1e7
#x = np.linspace(-scale, scale, shape[1])
#y = np.linspace(-scale, scale, shape[0])

#x2d, y2d = np.meshgrid(x, y)
#u = 10 * np.cos(2 * x2d / scale + 3 * y2d / scale)
#v = 20 * np.cos(6 * x2d / scale)

cmap = plt.cm.get_cmap("RdBu_r")
#levels = np.arange(0,15,0.1)
levels = np.arange(-5,5,0.1)
for i in range(10,35):
    fig = plt.figure(figsize=(6,5))
    ax = plt.axes( (0.1,0.1,0.7,0.8), projection=ccrs.PlateCarree())
    ax.coastlines('10m')
    ax.set_extent([138.5, 141.5, 34.5, 37.5], ccrs.PlateCarree())
    x, y , u, v = d.longitude.values, d.latitude.values, d.u.values[i], d.v.values[i]
    ws = divergence(u) + divergence(v) 
    #ws = (u**2 + v**2)**0.5
    c = ax.contourf(x,y,ws, levels = levels, cmap=cmap)
    q = ax.quiver(x, y, u, v, units='x', pivot='tip', width=0.022,
               scale=1 / 0.05) #, transform=vector_crs)
    clb = plt.colorbar(c, shrink=0.85) # draw colorbar
    clb.set_label('ms-1', labelpad=-10, y=1.05, rotation=0)
    #ax1.quiverkey(q, X=0.3, Y=1.1, U=10,
    #         label='Quiver key, length = 10', labelpos='E',transform = ax.transAxes)
    
    text = (pd.to_datetime(d.time[i].values)).strftime("%Y-%m-%d %H:%M")
    ax.text(0,1.02,text,transform = ax.transAxes,fontsize=14)
    
    #plt.show()
    fig.savefig("fig-3/"+(pd.to_datetime(d.time[i].values)).strftime("%Y-%m-%d_%H")+".png", dpi=150)
    plt.close()


'''
ax2 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
plt.title('The same vector field regridded')
ax2.coastlines('50m')
ax2.set_extent([-45, 55, 20, 80], ccrs.PlateCarree())
ax2.quiver(x, y, u, v, transform=vector_crs, regrid_shape=20)
'''
plt.show()




