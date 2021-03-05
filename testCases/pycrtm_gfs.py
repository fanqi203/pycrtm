#!/usr/bin/env python3
#import ipdb
#ipdb.set_trace()
import configparser
import os, sys 
import numpy as np
from matplotlib import pyplot as plt
import xarray as xr
thisDir = os.path.dirname(os.path.abspath(__file__))
parentDir = os.path.dirname(thisDir)
sys.path.insert(0,parentDir)
from pyCRTM import pyCRTM, profilesCreate
from metpy.calc import *
import metpy
import xesmf as xe

fin="profiles/profile2d4_2019_dorain_gfs.nc"
fgrid="obs_grid.nc"


def cntrl(profiles):
    return profiles

def test_angle0_0(profiles):
    profiles.Angles[:,0]=0
    return profiles
def test_angle0_30(profiles):
    profiles.Angles[:,0]=30
    return profiles
def test_angle0_60(profiles):
    profiles.Angles[:,0]=60
    return profiles
def test_angle0_90(profiles):
    profiles.Angles[:,0]=90
    return profiles


def test_angle1_0(profiles):
    profiles.Angles[:,1]=0
    return profiles
def test_angle1_30(profiles):
    profiles.Angles[:,1]=30
    return profiles
def test_angle1_60(profiles):
    profiles.Angles[:,1]=60
    return profiles
def test_angle1_90(profiles):
    profiles.Angles[:,1]=90
    return profiles

def test_angle2_0(profiles):
    profiles.Angles[:,2]=0
    return profiles
def test_angle2_30(profiles):
    profiles.Angles[:,2]=30
    return profiles
def test_angle2_60(profiles):
    profiles.Angles[:,2]=60
    return profiles
def test_angle2_90(profiles):
    profiles.Angles[:,2]=90
    return profiles

def test_angle3_0(profiles):
    profiles.Angles[:,3]=0
    return profiles
def test_angle3_30(profiles):
    profiles.Angles[:,3]=30
    return profiles
def test_angle3_60(profiles):
    profiles.Angles[:,3]=60
    return profiles
def test_angle3_90(profiles):
    profiles.Angles[:,3]=90
    return profiles


def test_angle4_0(profiles):
    profiles.Angles[:,4]=0
    return profiles
def test_angle4_30(profiles):
    profiles.Angles[:,4]=30
    return profiles
def test_angle4_60(profiles):
    profiles.Angles[:,4]=60
    return profiles
def test_angle4_90(profiles):
    profiles.Angles[:,4]=90
    return profiles

def test_O3_0(profiles):
    profiles.O3[:,:]=0.0
    return profiles
def test_ice_0(profiles):
    profiles.clouds[:,:,1,0]=0.0
    return profiles
def test_ice_10pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*1.1
    return profiles
def test_ice_20pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*1.2
    return profiles
def test_ice_30pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*1.3
    return profiles
def test_ice_n10pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*0.9
    return profiles
def test_ice_n20pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*0.8
    return profiles
def test_ice_n30pct(profiles):
    profiles.clouds[:,:,1,0]=profiles.clouds[:,:,1,0]*0.7
    return profiles



def test_cld_0(profiles):
    profiles.clouds[:,:,0,0]=0.0
    return profiles
def test_cld_10pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*1.1
    return profiles
def test_cld_20pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*1.2
    return profiles
def test_cld_30pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*1.3
    return profiles
def test_cld_n10pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*0.9
    return profiles
def test_cld_n20pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*0.8
    return profiles
def test_cld_n30pct(profiles):
    profiles.clouds[:,:,0,0]=profiles.clouds[:,:,0,0]*0.7
    return profiles

def test_rain_0(profiles):
    profiles.clouds[:,:,2,0]=0.0
    return profiles
def test_rain_10pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*1.1
    return profiles
def test_rain_20pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*1.2
    return profiles
def test_rain_30pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*1.3
    return profiles
def test_rain_n10pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*0.9
    return profiles
def test_rain_n20pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*0.8
    return profiles
def test_rain_n30pct(profiles):
    profiles.clouds[:,:,2,0]=profiles.clouds[:,:,2,0]*0.7
    return profiles

def test_snow_0(profiles):
    profiles.clouds[:,:,3,0]=0.0
    return profiles
def test_snow_10pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*1.1
    return profiles
def test_snow_20pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*1.2
    return profiles
def test_snow_30pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*1.3
    return profiles
def test_snow_n10pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*0.9
    return profiles
def test_snow_n20pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*0.8
    return profiles
def test_snow_n30pct(profiles):
    profiles.clouds[:,:,3,0]=profiles.clouds[:,:,3,0]*0.7
    return profiles


def test_grpl_0(profiles):
    profiles.clouds[:,:,4,0]=0.0
    return profiles
def test_grpl_10pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*1.1
    return profiles
def test_grpl_20pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*1.2
    return profiles
def test_grpl_30pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*1.3
    return profiles
def test_grpl_n10pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*0.9
    return profiles
def test_grpl_n20pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*0.8
    return profiles
def test_grpl_n30pct(profiles):
    profiles.clouds[:,:,4,0]=profiles.clouds[:,:,4,0]*0.7
    return profiles


def test_R_ice_10pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*1.1
    return profiles
def test_R_ice_20pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*1.2
    return profiles
def test_R_ice_30pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*1.3
    return profiles
def test_R_ice_n10pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*0.9
    return profiles
def test_R_ice_n20pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*0.8
    return profiles
def test_R_ice_n30pct(profiles):
    profiles.clouds[:,:,1,1]=profiles.clouds[:,:,1,1]*0.7
    return profiles


def test_R_cld_10pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*1.1
    return profiles
def test_R_cld_20pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*1.2
    return profiles
def test_R_cld_30pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*1.3
    return profiles
def test_R_cld_n10pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*0.9
    return profiles
def test_R_cld_n20pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*0.8
    return profiles
def test_R_cld_n30pct(profiles):
    profiles.clouds[:,:,0,1]=profiles.clouds[:,:,0,1]*0.7
    return profiles

def test_R_rain_10pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*1.1
    return profiles
def test_R_rain_20pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*1.2
    return profiles
def test_R_rain_30pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*1.3
    return profiles
def test_R_rain_n10pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*0.9
    return profiles
def test_R_rain_n20pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*0.8
    return profiles
def test_R_rain_n30pct(profiles):
    profiles.clouds[:,:,2,1]=profiles.clouds[:,:,2,1]*0.7
    return profiles


def test_R_snow_10pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*1.1
    return profiles
def test_R_snow_20pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*1.2
    return profiles
def test_R_snow_30pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*1.3
    return profiles
def test_R_snow_n10pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*0.9
    return profiles
def test_R_snow_n20pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*0.8
    return profiles
def test_R_snow_n30pct(profiles):
    profiles.clouds[:,:,3,1]=profiles.clouds[:,:,3,1]*0.7
    return profiles

def test_R_grpl_10pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*1.1
    return profiles
def test_R_grpl_20pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*1.2
    return profiles
def test_R_grpl_30pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*1.3
    return profiles
def test_R_grpl_n10pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*0.9
    return profiles
def test_R_grpl_n20pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*0.8
    return profiles
def test_R_grpl_n30pct(profiles):
    profiles.clouds[:,:,4,1]=profiles.clouds[:,:,4,1]*0.7
    return profiles

def test_T_sfc_plus1deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]+1.0
    return profiles

def test_T_sfc_plus2deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]+2.0
    return profiles

def test_T_sfc_min1deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-1.0
    return profiles

def test_T_sfc_min2deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-2.0
    return profiles

def test_T_sfc_min3deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-3.0
    return profiles

def test_T_sfc_min4deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-4.0
    return profiles

def test_T_sfc_min5deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-5.0
    return profiles

def test_T_sfc_min6deg(profiles):
    profiles.surfaceTemperatures[:,:]=profiles.surfaceTemperatures[:,:]-6.0
    return profiles


def test_wind_10pct(profiles):
    profiles.windSpeed10m[:]=profiles.windSpeed10m[:]*1.1
    return profiles

def test_wind_20pct(profiles):
    profiles.windSpeed10m[:]=profiles.windSpeed10m[:]*1.2
    return profiles

def test_wind_m10pct(profiles):
    profiles.windSpeed10m[:]=profiles.windSpeed10m[:]*0.9
    return profiles

def test_wind_m20pct(profiles):
    profiles.windSpeed10m[:]=profiles.windSpeed10m[:]*0.8
    return profiles

def test_lai_10pct(profiles):
    profiles.LAI[:]=profiles.LAI[:]*1.1
    return profiles

def test_lai_20pct(profiles):
    profiles.LAI[:]=profiles.LAI[:]*1.2
    return profiles

def test_lai_m10pct(profiles):
    profiles.LAI[:]=profiles.LAI[:]*0.9
    return profiles

def test_lai_m20pct(profiles):
    profiles.LAI[:]=profiles.LAI[:]*0.8
    return profiles

def test_aero_100pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002,0.0)
    profiles.aerosols[:,:,0,1]=2.0
    profiles.aerosolType[:]=1
    # 
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_10pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002*1.1,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_20pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002*1.2,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_m10pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002*0.9,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_m20pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002*0.8,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_R_10pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0*1.1
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_R_20pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0*1.2
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def test_aero_R_m10pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0*0.9
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles


def test_aero_R_m20pct(profiles):
    __tmp__=profiles.clouds[:,:,0,0]
    profiles.aerosols[:,:,0,0]=np.where(__tmp__>0.0,0.002,0.0) #0.002*1.1
    profiles.aerosols[:,:,0,1]=2.0*0.8
    profiles.aerosolType[:]=1
    #
    profiles.clouds[:,:,0,0]=0.0
    return profiles

def main(coefficientPath, sensor_id, fin, experiment):
    # the profile from grb2 has 33 levels 
    gfs=xr.open_dataset(fin)
    lat=gfs.lat_0
    lon=gfs.lon_0

    profiles = profilesCreate(len(lat)*len(lon),33, nAerosols=0, nClouds = 5)
    angles=gfs.xangles.stack(z=("lat_0","lon_0")).transpose()

#    profiles.Angles[:,:] = angles[:,:]
    profiles.Angles[:,0] =  angles[:,0] #h5['zenithAngle'][()]
    profiles.Angles[:,1] =  angles[:,1]
    profiles.Angles[:,2] =  angles[:,3]  # 100 degrees zenith below horizon.
    profiles.Angles[:,3] =  angles[:,4] # zero solar azimuth 
    profiles.Angles[:,4] =  angles[:,2] # h5['scanAngle'][()]

    # date time not used. 
    datetimes=gfs.valid_time.stack(z=("lat_0","lon_0")).transpose()
    profiles.DateTimes[:,:]=datetimes[:,:]
    
    profiles.Pi[:,:]=gfs.plevel.broadcast_like(gfs.player).stack(z=("lat_0","lon_0")).transpose()/100.0
    profiles.P[:,:]=gfs.player[1:,:,:].stack(z=("lat_0","lon_0")).transpose()/100.0

    profiles.T[:,:] = gfs.temp[1:,:,:].stack(z=("lat_0","lon_0")).transpose()

    Q=gfs.moisture[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    
    profiles.Q[:,:] = xr.where(Q>0, Q, 0) 
    
    profiles.O3[:,:]=gfs.o3[1:,:,:].stack(z=("lat_0","lon_0")).transpose()

    cld=gfs.water_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    cld=np.where(cld<0,0, cld)
    profiles.clouds[:,:,0,0]=cld 
    profiles.clouds[:,:,0,1]=10.0

    ice=gfs.ice_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    ice=np.where(ice<0,0,ice)
    profiles.clouds[:,:,1,0]=ice 
    profiles.clouds[:,:,1,1]=75

    rain=gfs.rain_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    rain=np.where(rain<0,0,rain)
    profiles.clouds[:,:,2,0]=rain 
    profiles.clouds[:,:,2,1]=50 

    snow=gfs.snow_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    snow=np.where(snow<0,0,snow)
    profiles.clouds[:,:,3,0]=snow 
    profiles.clouds[:,:,3,1]=50 

    graupel=gfs.graupel_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    graupel=np.where(graupel<0,0,graupel)
    profiles.clouds[:,:,4,0]=graupel 
    profiles.clouds[:,:,4,1]=50 


    profiles.cloudType[:,0]=1  
    profiles.cloudType[:,1]=2  
    profiles.cloudType[:,2]=3  
    profiles.cloudType[:,3]=4  
    profiles.cloudType[:,4]=5  

    profiles.climatology[:]=6 #US_STANDARD_ATMOSPHERE #6 # place holder don't know 


    lands=gfs.landmask.stack(z=("lat_0","lon_0")).transpose()
    snow=gfs.snow_cover.stack(z=("lat_0","lon_0")).transpose()
    ice=gfs.ice_cover.stack(z=("lat_0","lon_0")).transpose()

    profiles.surfaceFractions[:,0]=np.where(lands==1,int(1),int(0))
    profiles.surfaceFractions[:,1]=np.where(lands==0,int(1),int(0))
    profiles.surfaceFractions[:,2]=0 
    profiles.surfaceFractions[:,3]=0 

    profiles.surfaceTemperatures[:,0] =gfs.sfctemp.stack(z=("lat_0","lon_0"))
    profiles.surfaceTemperatures[:,1] =gfs.sfctemp.stack(z=("lat_0","lon_0"))
    profiles.surfaceTemperatures[:,2] =gfs.sfctemp.stack(z=("lat_0","lon_0"))
    profiles.surfaceTemperatures[:,3] =gfs.sfctemp.stack(z=("lat_0","lon_0"))
 
    profiles.S2m[:,1] = 33.0 # just use salinity out of S2m for the moment.
    wind_speed=gfs.wind_speed.stack(z=("lat_0","lon_0")).transpose()
    
    profiles.windSpeed10m[:] = wind_speed[:] # wind_speed

    # lai not affecting 
    profiles.LAI[:] = 0.57 # h5['LAI'][()] 
    wind_dir=gfs.wind_dir.stack(z=("lat_0","lon_0")).transpose()
    wind_dir=np.where(wind_dir<0,wind_dir+360,wind_dir)
    
    profiles.windDirection10m[:] = wind_dir[:] # h5['windDirection10m'][()]

    landtype=gfs.land_type[:,:,0].stack(z=("lat_0","lon_0")).transpose()
    landtype=np.where(landtype > 14, 16.0, landtype)
    landtype=np.where(landtype < 1.0, 1.0, landtype)
    
    profiles.surfaceTypes[:,0] = landtype 
    profiles.surfaceTypes[:,1] = 6 
    profiles.surfaceTypes[:,2] = 2 
    profiles.surfaceTypes[:,3] = 1 
    profiles.surfaceTypes[:,4] = 3 
    profiles.surfaceTypes[:,5] = 1 

    profiles=experiment(profiles)

    try:
        crtmOb = pyCRTM()
        crtmOb.profiles = profiles
        crtmOb.coefficientPath = pathInfo['CRTM']['coeffs_dir']
        crtmOb.sensor_id = sensor_id
        crtmOb.nThreads = 1
        
        crtmOb.loadInst()
        crtmOb.runDirect()
        forwardTb = crtmOb.Bt
    except:
        print("error in ", experiment.__name__)
        #    forwardEmissivity = crtmOb.surfEmisRefl[0,:]
        #    #make K matrix run surfEmisRefl again.
        #    crtmOb.surfEmisRefl = []
        #    crtmOb.runK()
        #    kTb = crtmOb.Bt
        #    kEmissivity = crtmOb.surfEmisRefl[0,:]
    finally:
        result=xr.DataArray(np.reshape(forwardTb.T,(10,len(lat),len(lon))),dims=("channel","lat_0","lon_0"), coords=(np.arange(0,10),gfs.lat_0,gfs.lon_0))


        #fig,ax=plt.subplots(5,2)
        #    result[:,30,30].plot()
        #    for i in range(5):
        #        for j in range(2):
        #            ax[i,j].pcolor(result[i*2+j,:,:].where(result[i*2+j,:,:]>-100), vmin=201, vmax=300)
        #


        # prepare for interpolation 
        ds_in=result.rename({"lat_0": 'lat', "lon_0": 'lon'})    
        
        fo=xr.open_dataset(fgrid)
        lato=fo.lat.data
        lono=fo.lon.data

        ds_out = xr.Dataset({'lat': (['lat'], lato),'lon': (['lon'], lono),})
        regridder = xe.Regridder(ds_in, ds_out, 'bilinear',reuse_weights=True)

        dr_out = regridder(result)
        out = dr_out.where(fo.obs.data>0)
        
        fout=fin[0:len(fin)-3]+"_output_"+experiment.__name__+".nc"

        if os.path.exists(fout):
            os.remove(fout)
        result.to_netcdf(fout,"w")
    
    
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(parentDir,'crtm.cfg') ) 
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'abi_g16'
    main(coefficientPath, sensor_id, fin, EXP)
