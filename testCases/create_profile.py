#!/usr/bin/env python3
import configparser
import os, sys
from matplotlib import pyplot as plt
import xarray as xr
thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(thisDir)
sys.path.insert(0,parentDir)

from metpy.calc import *
from metpy.units import units
import metpy
from skyfield import api
from skyfield.api import EarthSatellite, Topos, load
import numpy as np
from datetime import datetime, timedelta


f1="/Users/sbao/PycharmProjects/nggps_py/pycrtm/gfs.t06z.pgrb2.0p25.f120"
f1="/Users/sbao/PycharmProjects/nggps_py/pycrtm/gfs.t00z.pgrb2.0p25.f000.2020"

n_clouds=5
n_aerosol=0

def angle_2d(lat,lon,y,m,d,h):
# given surface lat lon and time, calculate the its angles involving satellite and Sun 

    
    # information about satellites (chose GOES 16)
    sats = api.load.tle('https://celestrak.com/NORAD/elements/goes.txt')
    satellite = sats['GOES 16 [+]']

    # planets 
    planets = load('de421.bsp')
    sun = planets['sun']
    earth= planets['earth']

    # create array to hold angles 
    angles=np.zeros((len(lat),len(lon),5))

    for i in range(len(lat)):
        for j in range(len(lon)):
            #time 
            ts = api.load.timescale()
            tm = ts.tt(y,m,d,h,0,0)
            
            # call the surface station "boston"
            boston=Topos(str(lat[i])+' N', str(lon[j])+' W',elevation_m=0.0)

            # the angle between the two vectors: earth center to satellite and earth center to observer
            theta = satellite.at(tm).separation_from(boston.at(tm))

            # geometry 
            difference = satellite - boston
            geometry = difference.at(tm).altaz()

            # angles involving satellite 
            scan = np.round(180 - (90 + geometry[0].degrees + theta.degrees), 2)
            zenith = np.round(geometry[0].degrees, 2)
            azimuth = np.round(geometry[1].degrees, 2)
            angles[i,j,0]=zenith
            angles[i,j,1]=azimuth
            angles[i,j,2]=scan

            # angles involving the Sun 
            observer = earth + boston
            geo2 = observer.at(tm).observe(sun).apparent().altaz()
            zenithsun = np.round(geo2[0].degrees, 2)
            azimuthsun = np.round(geo2[1].degrees, 2)
            angles[i,j,3]=zenithsun
            angles[i,j,4]=azimuthsun
    return angles

def create_profile2d(f):
    # constants 
    Rd = 287.0
    Rv = 461.0
    fv = Rv / Rd - 1

    # open GRiB2 file 
    gfs = xr.open_dataset(f1, engine='pynio')
    y0=np.arange(20,60,1)
    x0=np.arange(200,300,1)
    n=len(y0)*len(x0)

    # p levels  
    pint=gfs.lv_ISBL0

    # heights 
    hgt=gfs.HGT_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')

    # y m d h
    init_time = hgt.initial_time

    month=int(init_time.split('/')[0])
    day = int(init_time.split('/')[1])
    year =int(init_time.split('/')[2].split(" ")[0])
    hour =int(init_time.split('/')[2].split(" ")[1][1:3])

    fcst_time = hgt.forecast_time[0]

    valid_time = datetime(year,month,day,hour)+timedelta(hours=int(fcst_time))



    # water vapor levels 
    # qp=gfs.lv_ISBL5

    # relative humidity
    q=gfs.RH_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    #  interploate q to all levels (upper levels has no moisture)
    qpint=q.interp(lv_ISBL5=pint,kwargs={"fill_value": 0.0})

    # temp 
    t = gfs.TMP_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')

    # mixing ratio 
    mix = mixing_ratio_from_relative_humidity(pint.broadcast_like(t),t,qpint)
    mix=xr.DataArray(mix,coords=t.coords,dims=t.dims)
    mix.attrs['units']='kg/kg'
    mixavg=mix.rolling(lv_ISBL0=2,center=True).mean()[1:]

    # get the t and p for the layers (instead of levels)
    tavg=t.rolling(lv_ISBL0=2,center=True).mean()[1:]
    dp=xr.DataArray(np.diff(pint),dims=pint.dims).broadcast_like(tavg)
    dz=xr.DataArray(np.diff(hgt,axis=0),dims=hgt.dims).broadcast_like(tavg)
    R=287.05
    g=9.8
    pavg=-R*tavg*(1+mixavg*fv)*dp/dz/g
    rho=-(1+mixavg*fv)*dp/dz/g

    # ozone and the five types of "clouds"
    o3=gfs.O3MR_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    o3_pavg=o3.interp(lv_ISBL12=pavg.coords['lv_ISBL0'],kwargs={"fill_value": 0.0})

    cld=gfs.CLWMR_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    cld_pavg=cld.interp(lv_ISBL7=pavg.coords['lv_ISBL0'],kwargs={"fill_value": 0.0})
    cld_wc=cld_pavg*rho*dz
    print(type(cld_pavg))
    print(type(cld_wc))

    ice_cld=gfs.ICMR_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    ice_cld_pavg = ice_cld.interp(lv_ISBL7=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    ice_wc=ice_cld_pavg*rho*dz
    
    rain_cld=gfs.RWMR_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    rain_cld_pavg = rain_cld.interp(lv_ISBL7=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    rain_wc=rain_cld_pavg*rho*dz

    
    snow_cld=gfs.SNMR_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    snow_cld_pavg = snow_cld.interp(lv_ISBL7=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    snow_wc=snow_cld_pavg*rho*dz
    
    grp_cld=gfs.GRLE_P0_L100_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    grp_cld_pavg=grp_cld.interp(lv_ISBL7=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    grp_wc=grp_cld_pavg*rho*dz
    
    year=valid_time.year
    month=valid_time.month
    day=valid_time.day
    hour=valid_time.hour

    angles = angle_2d(t.coords['lat_0'].values, t.coords['lon_0'].values, year,month,day,hour)
    xangles = xr.DataArray(angles, name="xangles", dims=[t.dims[1], t.dims[2], 'angles'],
                          coords=[t.coords['lat_0'], t.coords['lon_0'], np.arange(0, 5)])

    datetimes=xr.DataArray(np.zeros((len(y0),len(x0),6)),dims=[t.dims[1],t.dims[2],'ymd'],coords=[t.coords['lat_0'],t.coords['lon_0'],np.arange(0,6)])
    datetimes[:, :, 0] = year
    datetimes[:, :, 1] = month
    datetimes[:, :, 2] = day
    datetimes[:, :, 3] = hour
    datetimes[:, :, 4] = 0
    datetimes[:, :, 5] = 0

    # u and v 10m
    u10=gfs.UGRD_P0_L103_GLL0.sel(lv_HTGL8=10.0, lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    v10=gfs.VGRD_P0_L103_GLL0.sel(lv_HTGL8=10.0, lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    #u10=u10*units.meter_per_second
    #v10=v10*units.meter_per_second
    speed10m = np.sqrt(u10 * u10 + v10 * v10)
    dir10m=np.arctan2(v10, u10)/np.pi*180

    # land mask 
    lm=gfs.LAND_P0_L1_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    sm=1-lm
    snow=gfs.CSNOW_P0_L1_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')
    ice=gfs.ICEC_P0_L1_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')

    # surface temp 
    sfctemp=gfs.TMP_P0_L1_GLL0.sel(lat_0=np.arange(20,60,1), lon_0=np.arange(200,300,1),method='nearest')

    # surface type 
    sfctype=xr.DataArray(np.zeros((len(y0),len(x0),6)),dims=[t.dims[1],t.dims[2],'type'],coords=[t.coords['lat_0'],t.coords['lon_0'],np.arange(0,6)])

    # land type 
    landtype=xr.DataArray(np.zeros((len(y0),len(x0))),dims=[t.dims[1],t.dims[2]],coords=[t.coords['lat_0'],t.coords['lon_0']])

    # read surface type database 
    f = xr.open_dataset("MCD12C1.006.ncml.nc")
    modis_data = f.Majority_Land_Cover_Type_1[18][::-1, ::1]
    lat = np.arange(-90, 90, 0.05)
    lon = np.arange(-180, 180, 0.05)
    lon = np.where(lon<0,lon+360,lon)
    modis_data = modis_data.assign_coords({"Latitude": lat})
    modis_data = modis_data.assign_coords({"Longitude": lon})

    # interpolate database to model grid 
    landtype=modis_data.interp(Latitude=t.coords['lat_0'],Longitude=t.coords['lon_0'])

    sfctype[:,:,0]=landtype[:,:]
    sfctype[:,:,1:5]=0

    datetimes.name="valid_time"
    pint.name="plevel"
    pavg.name="player"
    tavg.name="temp"
    mixavg.name="moisture"
    o3_pavg.name="o3"
    cld_wc.name="water_cloud"
    ice_wc.name='ice_cloud'
    snow_wc.name = 'snow_cloud'
    rain_wc.name = 'rain_cloud'
    grp_wc.name = 'graupel_cloud'
    lm.name="landmask"
    sm.name="seamask"
    sfctemp.name="sfctemp"
    snow.name="snow_cover"
    ice.name="ice_cover"
    sfctype.name="land_type"
    speed10m.name='wind_speed'
    dir10m.name="wind_dir"
    all_data=xr.merge([xangles,datetimes,pint,pavg,tavg,mixavg,o3_pavg,cld_wc,ice_wc,snow_wc,rain_wc,grp_wc,lm,sfctemp,snow,ice,sfctype,speed10m,dir10m])
    print("combined")
    all_data.to_netcdf("profile2d4_2020.nc","w")


if __name__ == "__main__":
    create_profile2d(f1)

