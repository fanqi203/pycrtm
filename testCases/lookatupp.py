import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
x=xr.open_dataset("profiles/profile2d4_2019_dorain_gfs.nc")
hydro=x.water_cloud+x.ice_cloud+x.snow_cloud+x.rain_cloud+x.graupel_cloud
hydrosum=hydro.sum(dim='lv_ISBL0')
plt.pcolor(np.where((hydrosum>0) & (hydrosum>0.001),hydrosum,np.nan))
plt.pcolor(x.sfctemp.where(x.sfctemp>300).values)
ic=x.ice_cloud+x.water_cloud
ic_sum=ic.sum(dim='lv_ISBL0').values
plt.pcolor(np.where(ic_sum>0,ic_sum,np.nan),cmap='gist_ncar')
plt.pcolor(np.where(ic_sum>0.1,ic_sum,np.nan),cmap='gist_ncar')
plt.pcolor(np.where(ic_sum>0.01,ic_sum,np.nan),cmap='gist_ncar')
plt.show()

