import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import os
import os.path

class aero:
    def __init__(self,ymdh_ini,fcst):
        cef=xr.open_dataset("gfs_ctrl.nc")
        print("hello")
        self.hybrid_coefs=[cef.vcoord[0,1:65],cef.vcoord[1,1:65]]
        ymdh_str=str(ymdh_ini)
        y=ymdh_str[0:4]
        m=ymdh_str[4:6]
        d=ymdh_str[6:8]
        h=ymdh_str[8:10]
        url="https://noaa-gefs-pds.s3.amazonaws.com/gefs."+y+m+d+"/"+h+"/chem/pgrb2ap5/gefs.chem.t"+h+"z.a3d_0p50.f"+str(fcst).zfill(3)+".grib2"
        self.url=url
    def download(self):
        self.f=self.url.split('/')[-1]
        if not os.path.isfile(self.f):
            os.system("wget "+url)
            
    def sig2p(self):
        p0=1000
        hyam=self.hybrid_coefs[0]
        hybm=self.hybrid_coefs[1]
        psfc=1013 #xr.full_like(dust[0,:,:],1013)
        p=hyam/100.0+hybm*psfc
        self.p=p

    def extract_dust(self):
        gefs=xr.open_dataset(self.f,engine='pynio')
        tmp1=gefs.PMTC_P48_L105_GLL0_A62001
        tmp2=gefs.PMTC_P48_L105_GLL0_A62001_1
        tmp3=gefs.PMTC_P48_L105_GLL0_A62001_2
        tmp4=gefs.PMTF_P48_L105_GLL0_A62001
        tmp5=gefs.PMTF_P48_L105_GLL0_A62001_1
        
        dust=sum([tmp1,tmp2,tmp3,tmp4,tmp5])
        dust_size=(tmp1*5+tmp2*9+tmp3*11+tmp4*1+tmp5*3)/dust        

        dust=dust.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        dust_size=dust_size.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        pnew=np.arange(1000,90,-10)
        dust.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        dust_size.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        return (dust,dust_size)

    def extract_salt(self):
        gefs=xr.open_dataset(self.f,engine='pynio')
        tmp1=gefs.PMTC_P48_L105_GLL0_A62008
        tmp2=gefs.PMTC_P48_L105_GLL0_A62008_1
        tmp3=gefs.PMTC_P48_L105_GLL0_A62008_2
        tmp4=gefs.PMTF_P48_L105_GLL0_A62008 
        tmp5=gefs.PMTF_P48_L105_GLL0_A62008_1
        
        salt=sum([tmp1,tmp2,tmp3,tmp4,tmp5])
        salt_size=(tmp1*2+tmp2*6+tmp3*15+tmp4*0.6+tmp5*0.1)/salt        

        salt=salt.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        salt_size=salt_size.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        pnew=np.arange(1000,90,-10)
        salt.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        salt_size.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        return (salt,salt_size)

    def extract_Sulphate(self):
        gefs=xr.open_dataset(self.f,engine='pynio')
        tmp1=gefs.PMTF_P48_L105_GLL0_A62006
        
        sulp=tmp1
        sulp_size=xr.full_like(sulp,0.139)

        sulp=sulp.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        sulp_size=sulp_size.assign_coords({"lv_HYBL0":("lv_HYBL0",self.p[::-1].values)})
        pnew=np.arange(1000,90,-10)
        sulp.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        sulp_size.sel(lv_HYBL0=pnew,method="nearest").sel(lon_0=np.arange(250,300,1)).sel(lat_0=np.arange(20,60,1))[20].plot()
        plt.show()
        return (sulp,sulp_size)



a=aero(2021050412,96)
a.download()
a.sig2p()
#dust=a.extract_dust()
#sulp=a.extract_Sulphate()
salt=a.extract_salt()
