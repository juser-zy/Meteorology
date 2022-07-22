# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 14:56:59 2021

@author: lenovo
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt 
import numpy as np
import netCDF4 as nc
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader

file = r'ERA5_China.nc'
nc_data_object = nc.Dataset(file)

lons = nc_data_object.variables['longitude'][:]
lats = nc_data_object.variables['latitude'][:]
tcwv = np.asarray(nc_data_object.variables['tcwv'])[0]

lons_grid, lats_grid = np.meshgrid(lons, lats)

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 10
plt.rcParams['font.weight'] = 'bold'

china = shpreader.Reader('china1.shp').geometries()
proj = ccrs.PlateCarree()
fig = plt.figure(figsize=(6, 4.5))
ax = fig.add_subplot(1, 1 , 1, projection=proj)
ax.add_geometries(china, ccrs.PlateCarree(), facecolor='none', edgecolor='k', linewidth=1, zorder=1)

import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
cmap = mpl.cm.RdBu
newcolors=cmap(np.linspace(0, 1, 256))
newcmp = ListedColormap(newcolors[30:226])

cf = ax.contourf(lons_grid, lats_grid, tcwv, levels=50, cmap=newcmp, 
                 transform=proj, extend='both') 
cb = fig.colorbar(cf, shrink=0.70, orientation='horizontal', pad=0.08)
cb.ax.set_xlabel('TCWV(kg/$cm^2$)', fontweight='bold')
cb.ax.tick_params(which='major', direction='in', length=3)
gl = ax.gridlines(alpha=0.5, linestyle='--', draw_labels=True, 
                   dms=True, x_inline=False, y_inline=False)
gl.right_labels = 0
gl.top_labels = 0
ax.add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=1)

ax_sub = fig.add_axes([0.6615,0.3,0.12,0.2], projection=proj)  # [*left*, *bottom*, *width*,*height*]
ax_sub.set_extent([105,125,0,25], crs=ccrs.PlateCarree())
ax_sub.add_feature(cfeature.COASTLINE.with_scale('50m'))
china2 = shpreader.Reader('china1.shp').geometries()
ax_sub.add_geometries(china2, ccrs.PlateCarree(), facecolor='none', edgecolor='r', linewidth=1, zorder=1)

import maskout
clip=maskout.shp2clip(cf, ax, r'F:\Py_project\geoVisiual\hls_cartopy_china_nanhai\china0') # 白化

