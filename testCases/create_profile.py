#!/usr/bin/env python3
import ipdb
#ipdb.set_trace()
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

fobs="obs_grid.nc"
f1="tar/gfs.t18z.pgrb2.0p25.f003"
f1="tar/gfs.t12z.pgrb2.0p25.f006"
f1="/work/noaa/dtc-hwrf/sbao/EMC_post/DOMAINPATH/postprd/GFSPRS.006.iveg2.grb2"
n_clouds=5
n_aerosol=1


def find_var(f,string):
    for i in f:
        if (string in i and i.startswith(string[0:3])):
            found=i
    return found


def effr_cld(cld):
    min_qc=1.e-7
    gfdl_rhor=1000.
    gfdl_ccn=1.0e8
    gfdl_rewmin=5.0
    gfdl_rewmax=10.0
    cld_positive=xr.where(cld>min_qc, cld, min_qc)
    result=np.exp(1.0/3.0*np.log((3.*cld_positive)/(4.*np.pi*gfdl_rhor*gfdl_ccn)))*1.0e6
    result=xr.where(result<gfdl_rewmin, gfdl_rewmin, result)
    result=xr.where(result>gfdl_rewmax, gfdl_rewmax, result)
    result=xr.where(cld>min_qc, result, 0)
    return result*2

def effr_ice(ice,temp):
    min_qi=1.e-8
    gfdl_tice = 273.16
    gfdl_beta = 1.22
    gfdl_reimin = 10.0
    gfdl_reimax = 150.0
    result=xr.full_like(ice,0.0)
    ice_positive=xr.where(ice>min_qi, ice, min_qi)
    result=xr.where(temp[1:]-gfdl_tice>=-30.0, gfdl_beta / 9.387 * np.exp ((1 - 0.969) * np.log (1.0e3 * ice_positive)) * 1.0e3, result)
    result=xr.where(temp[1:]-gfdl_tice<-30.0, gfdl_beta / 9.208 * np.exp ((1 - 0.945) * np.log (1.0e3 * ice_positive)) * 1.0e3, result)
    result=xr.where(temp[1:]-gfdl_tice<-40.0, gfdl_beta / 9.337 * np.exp ((1 - 0.920) * np.log (1.0e3 * ice_positive)) * 1.0e3, result)
    result=xr.where(temp[1:]-gfdl_tice<-50.0, gfdl_beta / 9.917 * np.exp ((1 - 0.891) * np.log (1.0e3 * ice_positive)) * 1.0e3, result)
    result=xr.where(result<gfdl_reimin, gfdl_reimin, result)
    result=xr.where(result>gfdl_reimax, gfdl_reimax, result)
    result=xr.where(ice>min_qi, result, 0)
    return result*2

def effr_rain(rain):
    gfdl_rhor=1000.
    gfdl_n0r=8.e6
    min_qr=1.e-7
    gfdl_gammar=17.837789
    gfdl_alphar=0.8
    gfdl_rermin=0.0
    gfdl_rermax=10000.0
    rain_positive=xr.where(rain>min_qr, rain, min_qr)
    lamdbar=np.exp(0.25*np.log(np.pi*gfdl_rhor*gfdl_n0r/rain_positive))
    result=0.5*np.exp(np.log(gfdl_gammar /6.0)/gfdl_alphar)/lamdbar*1.0e6
    result=xr.where(result<gfdl_rermin, gfdl_rermin, result)
    result=xr.where(result>gfdl_rermax, gfdl_rermax, result)
    result=xr.where(rain>min_qr, result, 0)
    return result
        

def effr_snow(snow):
    gfdl_rhos=100.
    gfdl_n0s=3.e6
    min_qs=1.e-8
    gfdl_gammas=8.2850630
    gfdl_alphas=0.25
    gfdl_resmin=0.0
    gfdl_resmax=10000.0
    snow_positive=xr.where(snow>min_qs, snow, min_qs)
    lambdas=np.exp(0.25*np.log(np.pi*gfdl_rhos*gfdl_n0s/snow_positive))
    result=0.5*np.exp(np.log(gfdl_gammas /6.0)/gfdl_alphas)/lambdas*1.0e6
    result=xr.where(result<gfdl_resmin, gfdl_resmin, result)
    result=xr.where(result>gfdl_resmax, gfdl_resmax, result)
    result=xr.where(snow>min_qs, result, 0)
    return result

def effr_grp(grp):
    gfdl_rhog=400.
    gfdl_n0g=4.e6
    min_qg=1.e-7
    gfdl_gammag=11.631769
    gfdl_alphag=0.5
    gfdl_regmin=0.0
    gfdl_regmax=10000.0
    grp_positive=xr.where(grp>min_qg, grp, min_qg)
    lambdag=np.exp(0.25*np.log(np.pi*gfdl_rhog*gfdl_n0g/grp_positive))
    result=0.5*np.exp(np.log(gfdl_gammag /6.0)/gfdl_alphag)/lambdag*1.0e6
    result=xr.where(result<gfdl_regmin, gfdl_regmin, result)
    result=xr.where(result>gfdl_regmax, gfdl_regmax, result)
    result=xr.where(snow>min_qg, result, 0)
    return result

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
            boston=Topos(str(lat[i])+' N', str(lon[j])+' E',elevation_m=0.0)

            # the angle between the two vectors: earth center to satellite and earth center to observer
            theta = satellite.at(tm).separation_from(boston.at(tm))

            # geometry 
            difference = satellite - boston
            geometry = difference.at(tm).altaz()

            # angles involving satellite 
            scan = np.round(180 - (90 + geometry[0].degrees + theta.degrees), 2)
            zenith = 90-np.round(geometry[0].degrees, 2)
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

def create_profile2d(f,fobs):
    # f is the GRiB2 file
    # fo is the file with the target grid
    # constants 
    Rd = 287.0
    Rv = 461.0
    fv = Rv / Rd - 1

    # open GRiB2 file 
    gfs = xr.open_dataset(f1, engine='pynio')

    y0=np.arange(20.0,45.0,0.25)
    x0=np.arange(255,290,0.25)
    n=len(y0)*len(x0)

    # p levels  
    pint=gfs.lv_ISBL0


    # heights 
    hgt=gfs[find_var(gfs,"HGT_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')

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
    qq=gfs[find_var(gfs,"RH_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    #  interploate q to all levels (upper levels has no moisture)
    #  qpint=qq.interp(lv_ISBL5=pint,kwargs={"fill_value": 0.0})
    qpint=qq # in new input file qq are on all p levels

    # temp 
    t = gfs[find_var(gfs,"TMP_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')

    # mixing ratio from specific humidity

    sphd=gfs[find_var(gfs,"SPFH_P0_L100_")].sel(lat_0=y0,lon_0=x0,method='nearest')
    mix=sphd/(1-sphd)*1000
    
    # the mixing ration from RH is saved here but not used.
    mix_metpy = mixing_ratio_from_relative_humidity(pint.broadcast_like(t),t,qpint)

    mix=xr.DataArray(mix,coords=t.coords,dims=t.dims)
    mix.attrs['units']='kg/kg'

    # test another way to average mix
    mixavg=mix.rolling(lv_ISBL0=2,center=True).mean()[1:]

    # get the t and p for the layers (instead of levels)
    tavg=t.rolling(lv_ISBL0=2,center=True).mean()[1:]
    dp=xr.DataArray(np.diff(pint),dims=pint.dims).broadcast_like(tavg)
    dz=xr.DataArray(np.diff(hgt,axis=0),dims=hgt.dims).broadcast_like(tavg)
    R=287.05
    g=9.8
    pavg=-R*tavg*(1+mixavg*fv)*dp/dz/g
    rho=-(1+mixavg*fv)*dp/dz/g

    mixavg=mixavg #2000.0 #rho*-dz
    print(mixavg.max())
    # ozone and the five types of "clouds"
    o3=gfs[find_var(gfs,"O3MR_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    #o3_pavg=o3.interp(lv_ISBL12=pavg.coords['lv_ISBL0'],kwargs={"fill_value": 0.0})

    o3_pavg=o3*rho*-dz

#    cld=gfs.CLWMR_P0_L100_GLL0.sel(lat_0=y0, lon_0=x0,method='nearest')
    cld=gfs[find_var(gfs,"CLWMR_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    cld=cld.rename({cld.dims[0]: 'lv_ISBL_HYD'})
    cld_pavg=cld.interp(lv_ISBL_HYD=pavg.coords['lv_ISBL0'],kwargs={"fill_value": 0.0})
    cld_effr=effr_cld(cld_pavg)
    cld_wc=cld_pavg*rho*-dz



    ice_cld=gfs[find_var(gfs,"ICMR_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    ice_cld=ice_cld.rename({ice_cld.dims[0]: 'lv_ISBL_HYD'})
    ice_cld_pavg = ice_cld.interp(lv_ISBL_HYD=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    ice_effr=effr_ice(ice_cld_pavg,t)
    ice_wc=ice_cld_pavg*rho*-dz
    
    
    rain_cld=gfs[find_var(gfs,"RWMR_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    rain_cld=rain_cld.rename({rain_cld.dims[0]: 'lv_ISBL_HYD'})
    rain_cld_pavg = rain_cld.interp(lv_ISBL_HYD=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    rain_effr=effr_rain(rain_cld_pavg)
    rain_wc=rain_cld_pavg*rho*-dz

    
    snow_cld=gfs[find_var(gfs,"SNMR_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    snow_cld=snow_cld.rename({snow_cld.dims[0]: 'lv_ISBL_HYD'})
    snow_cld_pavg = snow_cld.interp(lv_ISBL_HYD=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    snow_effr=effr_snow(snow_cld_pavg)
    snow_wc=snow_cld_pavg*rho*-dz
    
    grp_cld=gfs[find_var(gfs,"GRLE_P0_L100_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    grp_cld=grp_cld.rename({grp_cld.dims[0]: 'lv_ISBL_HYD'})
    grp_cld_pavg=grp_cld.interp(lv_ISBL_HYD=pavg.coords['lv_ISBL0'], kwargs={"fill_value": 0.0})
    grp_effr=effr_snow(grp_cld_pavg)
    grp_wc=grp_cld_pavg*rho*-dz
    
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
    u10=gfs[find_var(gfs,"UGRD_P0_L103_")]
    u10=u10.rename({u10.dims[0]:'lv_HTGL_wind'})
    u10=u10.sel(lv_HTGL_wind=10.0, lat_0=y0, lon_0=x0,method='nearest')
    v10=gfs[find_var(gfs,"VGRD_P0_L103_")]
    v10=v10.rename({v10.dims[0]:'lv_HTGL_wind'})
    v10=v10.sel(lv_HTGL_wind=10.0, lat_0=y0, lon_0=x0,method='nearest')
    
    print(u10.max(),u10.min())

    speed10m = np.sqrt(u10 * u10 + v10 * v10)
    dir10m=np.arctan2(v10, u10)/np.pi*180

    # land mask 
    lm=gfs[find_var(gfs,"LAND_P0_L1_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    sm=1-lm
    snow=gfs[find_var(gfs,"CSNOW_P0_L1_")].sel(lat_0=y0, lon_0=x0,method='nearest')
    ice=gfs[find_var(gfs,"ICEC_P0_L1_")].sel(lat_0=y0, lon_0=x0,method='nearest')

    print(lm.max(),lm.min())
    print(sm.max(),sm.min())
    print(snow.max(),snow.min())
    print(ice.max(),ice.min())

    # surface temp 
    sfctemp=gfs[find_var(gfs,"TMP_P0_L1_")].sel(lat_0=y0, lon_0=x0,method='nearest')

    print(sfctemp)

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
    cld_effr.name="water_cloud_effr"
    ice_wc.name='ice_cloud'
    ice_effr.name="ice_cloud_effr"
    snow_wc.name = 'snow_cloud'
    snow_effr.name='snow_cloud_effr'
    rain_wc.name = 'rain_cloud'
    rain_effr.name='rain_cloud_effr'
    grp_wc.name = 'graupel_cloud'
    grp_effr.name='graupel_cloud_effr'
    lm.name="landmask"
    sm.name="seamask"
    sfctemp.name="sfctemp"
    snow.name="snow_cover"
    ice.name="ice_cover"
    sfctype.name="land_type"
    speed10m.name='wind_speed'
    dir10m.name="wind_dir"
    all_data=xr.merge([xangles,datetimes,pint,pavg,tavg,mixavg,o3_pavg,cld_wc,cld_effr,ice_wc,ice_effr,snow_wc,snow_effr,rain_wc,rain_effr,grp_wc,grp_effr,lm,sfctemp,snow,ice,sfctype,speed10m,dir10m])
    print("combined")
    all_data.to_netcdf("profile2d4_2021_gfs_EMCUPP_qv1000.nc","w")


if __name__ == "__main__":
    create_profile2d(f1,fobs)

