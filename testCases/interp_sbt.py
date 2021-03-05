import ipdb
ipdb.set_trace()
import xarray as xr
import xesmf as xe
fin="profiles/profile2d4_2019_dorain_gfs.nc"
gfs=xr.open_dataset(fin)
lat=gfs.lat_0
lon=gfs.lon_0
ds_out=gfs.rename({"lat_0": 'lat', "lon_0": 'lon'})

fsbt="profiles/sbtdorian.nc"
ds_in=xr.open_dataset(fsbt).SBT124_P0_L8_GLC0

regridder = xe.Regridder(ds_in, ds_out, 'bilinear',reuse_weights=True)
dr_out = regridder(ds_in)
dr_out.to_netcdf("sbt_dorian_regrid.nc")
