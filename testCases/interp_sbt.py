#import ipdb
#ipdb.set_trace()
import xarray as xr
import xesmf as xe
import os

fin="profiles/profile2d4_2019_dorain_gfs.nc"
gfs=xr.open_dataset(fin)
lat=gfs.lat_0
lon=gfs.lon_0
ds_out=gfs.rename({"lat_0": 'lat', "lon_0": 'lon'})

fsbt="profiles/sbtdorian.nc"
fsbt="GFSPRS.006"
fsbt=os.getenv('new_sbt')
ds_in=xr.open_dataset(fsbt,engine='pynio').rename({"lat_0": 'lat', "lon_0": 'lon'})

ds_in_ch13=ds_in.VAR_3_192_58_P0_L8_GGA0
ds_in_sbt=ds_in.SBT124_P0_L8_GGA0
#ds_in2=xr.open_dataset(fsbt,engine='pynio').VAR_3_192_59_P0_L8_GGA0

regridder = xe.Regridder(ds_in, ds_out, 'bilinear',reuse_weights=True)
dr_out_ch13 = regridder(ds_in_ch13)
dr_out_sbt=regridder(ds_in_sbt)
dr_out_ch13.to_netcdf("sbt_2021_regrid.nc")
dr_out_sbt.to_netcdf("sbt_2021_sbt_regrid.nc")
