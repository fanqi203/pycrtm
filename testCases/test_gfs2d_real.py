#!/usr/bin/env python3
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


f1="profile2d4_2020.nc"



def main(coefficientPath, sensor_id):
    # the profile from grb2 has 33 levels 
    profiles = profilesCreate(4000,33, nAerosols=0, nClouds = 5)
    #gfs=xr.open_dataset(f1,engine='pynio')
    gfs=xr.open_dataset(f1)

    angles=gfs.xangles.stack(z=("lat_0","lon_0")).transpose()

#    profiles.Angles[:,:] = angles[:,:]
    profiles.Angles[:,0] = angles[:,0] #h5['zenithAngle'][()]
    profiles.Angles[:,1] = angles[:,1]
    profiles.Angles[:,2] = angles[:,3]  # 100 degrees zenith below horizon.
    profiles.Angles[:,3] = angles[:,4] # zero solar azimuth 
    profiles.Angles[:,4] = angles[:,2] # h5['scanAngle'][()]

    # date time not used. 
    datetimes=gfs.valid_time.stack(z=("lat_0","lon_0")).transpose()
    profiles.DateTimes[:,:]=datetimes[:,:]
    
    profiles.Pi[:,:]=gfs.plevel.broadcast_like(gfs.player).stack(z=("lat_0","lon_0")).transpose()/100.0
    profiles.P[:,:]=gfs.player[1:,:,:].stack(z=("lat_0","lon_0")).transpose()/100.0

    profiles.T[:,:] = gfs.temp[1:,:,:].stack(z=("lat_0","lon_0")).transpose()




    profiles.Q[:,:] = gfs.moisture[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.O3[:,:]=gfs.o3[1:,:,:].stack(z=("lat_0","lon_0")).transpose()

    cld=gfs.water_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    cld=np.where(cld>0,0, -cld)
    print(cld.shape)
    print(cld.max())
    print(cld.min())

#    fig, ax = plt.subplots()
#    ax.pcolormesh(cld.sum(dim=')
#    #plt.plot(profiles.T)
#    plt.show()

    profiles.clouds[:,:,0,0]=cld #fs.cldavg[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.clouds[:,:,0,1]=10.0

    ice=gfs.ice_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    ice=np.where(ice>0,0,-ice)
    profiles.clouds[:,:,1,0]=ice #fs.cldavg[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.clouds[:,:,1,1]=20.0

    rain=gfs.rain_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    rain=np.where(rain>0,0,-rain)
    profiles.clouds[:,:,2,0]=rain #fs.cldavg[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.clouds[:,:,2,1]=50.0

    snow=gfs.snow_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    snow=np.where(snow>0,0,-snow)
    profiles.clouds[:,:,3,0]=snow #fs.cldavg[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.clouds[:,:,3,1]=50.0

    graupel=gfs.graupel_cloud[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    graupel=np.where(graupel>0,0,-graupel)
    profiles.clouds[:,:,4,0]=graupel #fs.cldavg[1:,:,:].stack(z=("lat_0","lon_0")).transpose()
    profiles.clouds[:,:,4,1]=50.0


    profiles.cloudType[:,0]=1 #place holder 
    profiles.cloudType[:,1]=2 #place holder 
    profiles.cloudType[:,2]=3 #place holder 
    profiles.cloudType[:,3]=4 #place holder 
    profiles.cloudType[:,4]=5 #place holder 

    #print(profiles.cloudType.shape)

    profiles.climatology[:]=6 #US_STANDARD_ATMOSPHERE #6 # place holder don't know 


    lands=gfs.landmask.stack(z=("lat_0","lon_0")).transpose()
    snow=gfs.snow_cover.stack(z=("lat_0","lon_0")).transpose()
    ice=gfs.ice_cover.stack(z=("lat_0","lon_0")).transpose()

    #gfs.ice_cover.plot()
    #plt.show()
    print(np.sum(snow))
    print(np.sum(ice))
  
    #for i in np.arange(0,4000):
#        land_frac=1 
#        water_frac=0 
#        snow_frac=0 
#        ice_frac=0 
#    profiles.surfaceFractions[:,:]=np.array([1, 0, 0, 0]) 


    #land_frac=np.where((lands==1) & (snow == 0) & (ice == 0) ,int(1),int(0))
    #water_frac=np.where(lands==0,int(1),int(0))
    #snow_frac=np.where((lands==1) & (snow > 0),int(1),int(0))
    #ice_frac=np.where((lands==1) & (ice > 0),int(1),int(0))
    #profiles.surfaceFractions[i,:]=np.array([0,1,0,0])

    profiles.surfaceFractions[:,0]=np.where(lands==1,int(1),int(0))
    profiles.surfaceFractions[:,1]=np.where(lands==0,int(1),int(0))
    profiles.surfaceFractions[:,2]=0 #np.where((lands==1) & (snow > 0),int(1),int(0))
    profiles.surfaceFractions[:,3]=0 #np.where((lands==1) & (ice > 0),int(1),int(0))

#    profiles.surfaceFractions[:,0]=1 #np.where(lands>0, int(1),int(0))
#
#    profiles.surfaceFractions[:,1]=0 #np.where(lands==0,int(1),int(0))
#
#    profiles.surfaceFractions[:,2]=0 #np.where((lands==1) & (snow > 0),int(1),int(0))
#
#    profiles.surfaceFractions[:,3]=0 # np.where((lands==1) & (ice > 0),int(1),int(0))


    
#    print(type(profiles.surfaceFractions))
#    print(np.sum(profiles.surfaceFractions))
#
#    quit()
#    print(np.sum(profiles.surfaceFractions))
#
#    print(profiles.surfaceFractions)
#    
#    print(profiles.surfaceFractions.shape)
#
#
#    plt.contourf(profiles.surfaceFractions[:,:])
#    plt.show()



 
 #   print(profiles.surfaceFractions[:,1])


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
    
    print(landtype.shape)
#    quit()
    #plt.plot(landtype)
    #plt.show()
    #plt.show()
    print(landtype.max())
    print(landtype.min())
#    landtype[0:2001]=16.0
#    landtype[2001:3001]=14.0
#    landtype[3001:4000]=13.0
    landtype=np.where(landtype > 14, 16.0, landtype)
    landtype=np.where(landtype<1.0, 1.0, landtype)
    print(any(landtype[:]==15))
    print(landtype.max())
    print(landtype.min())
    print(landtype.astype(int))
    
#    quit()
#    landtype=np.where(landtype==15.0, 1, landtype)
#    landtype=np.where(landtype==0.0, 1, landtype)

#    print(any(landtype[:,0]==0.0))
#    print(any(landtype[:,0]==15.0))

    #plt.plot(landtype,"+")
    #plt.show()
    
    profiles.surfaceTypes[:,0] = landtype #h5['landType'][()]
    profiles.surfaceTypes[:,1] = 6 #landtype[:,1] # h5['soilType'][()]
    profiles.surfaceTypes[:,2] = 2 #h5['vegType'][()]
    profiles.surfaceTypes[:,3] = 1 #h5['waterType'][()]
    profiles.surfaceTypes[:,4] = 3 # ['snowType'][()]
    profiles.surfaceTypes[:,5] = 1 # h5['iceType'][()]


    crtmOb = pyCRTM()
    crtmOb.profiles = profiles
    crtmOb.coefficientPath = pathInfo['CRTM']['coeffs_dir']
    crtmOb.sensor_id = sensor_id
    crtmOb.nThreads = 1

    crtmOb.loadInst()
    
    crtmOb.runDirect()
    forwardTb = crtmOb.Bt
    forwardEmissivity = crtmOb.surfEmisRefl[0,:]
    #make K matrix run surfEmisRefl again.
    crtmOb.surfEmisRefl = []
    crtmOb.runK()
    kTb = crtmOb.Bt
    kEmissivity = crtmOb.surfEmisRefl[0,:]

    print(forwardTb.T.shape)
    
    result=xr.DataArray(np.reshape(forwardTb.T,(10,40,100)),dims=("channel","lat_0","lon_0"), coords=(np.arange(0,10),gfs.lat_0,gfs.lon_0))
    print(result)
    #fig,ax=plt.subplots(5,2)
    #    result[:,30,30].plot()
#    for i in range(5):
#        for j in range(2):
#            ax[i,j].pcolor(result[i*2+j,:,:].where(result[i*2+j,:,:]>-100), vmin=201, vmax=300)
#

    plt.pcolor(result[8]) #,vmin=230,vmax=290)
    
    plt.show()
    print(result.max())
    print(result.min())
#    plt.plot(forwardTb[:,:])
#    plt.show()

            
    

    
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(parentDir,'crtm.cfg') ) 
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    sensor_id = 'abi_g16'
    main(coefficientPath, sensor_id)
