# run ipdb debug:

#import ipdb
#ipdb.set_trace()

import xarray as xr
import matplotlib
#matplotlib.use('AGG')# or PDF, SVG or PS
 
import matplotlib.pyplot as plt

import numpy as np

import os 
from pathlib import Path

def preprocess(ds):
    import os
    ds=ds.expand_dims("cases")
    filename=ds.encoding["source"]
    fn=os.path.basename(filename)
    ds=ds.assign_coords(cases=("cases",[fn[39:len(fn)-3]]))
    ds=ds.rename({"__xarray_dataarray_variable__":"bt"})
    print(ds.bt[0,7,:,:].mean().values)
    return ds

def preprocessobs(ds):
    import os
    ds=ds.expand_dims("cases")
    filename=ds.encoding["source"]
    fn=os.path.basename(filename)
    ds=ds.assign_coords(cases=("cases",[fn[39:len(fn)-3]]))
    #ds=ds.rename({"__xarray_dataarray_variable__":"bt"})
    return ds


def hist(ch7):
    ch7flat=ch7.where(ch7>0).values.flatten()
    [obshist, obsbin] = np.histogram(ch7flat, density=True, bins=np.arange(190, 320,3))
    avg_BT=ch7flat[np.where(ch7flat>0)].mean()
    labelstr=str(np.round(avg_BT,0))+' '+str(ch7.cases.values)
    return [obsbin, obshist,  labelstr]

def plot(cases, figname):
    
# comment out this for now
#    f2=cases #"profiles/profile2d4_2019_dorain_gfs_output_test_angle*.nc"
#    xall=xr.open_mfdataset(f2,preprocess=preprocess,concat_dim="cases")
#    ch7=xall.bt

    f0="profiles/profile2d4_2019_dorain_gfs_output_test_cld_30pct.nc"
    f0="profiles/profile2d4_2019_dorain_gfs_output_test_R_ice_30pct.nc"
    f0="profiles/profile2d4_2019_dorain_gfs_output_test_T_sfc_min6deg"
    f0="profiles/profile2d4_2019_dorain_gfs_output_hydrozero.nc"
    f0=cases # "profiles/profile2d4_2019_dorain_gfs_output_cntrl.nc"
    x=xr.open_mfdataset(f0,preprocess=preprocess,concat_dim="cases")
    ch7_cntrl=0.5*(x.bt[0,6]+x.bt[0,7])


    fobs="profiles/obs_dorian*.nc"
    obs=xr.open_mfdataset(fobs,preprocess=preprocessobs,concat_dim="cases")
    ch7_obs=obs.obs


    fsbt="profiles/sbt_dorian_regrid.nc"
    sbt=xr.open_mfdataset(fsbt,preprocess=preprocessobs,concat_dim="cases")
    ch7_sbt=sbt.SBT124_P0_L8_GLC0


    fig,(ax0,ax1,ax2) = plt.subplots(3, 1)

#    for i in np.arange(ch7.shape[0]):
#        ax.plot(hist(ch7[i,7])[0][1:],hist(ch7[i,7])[1][0:],label=hist(ch7[i,7])[2]) 
    
    def cond1(var,st,ed):
        return np.where((var > st) & (var < ed), var, np.nan)

    st=300
    ed=320
    ax0.pcolor(cond1(ch7_cntrl,st,ed), vmin=220, vmax=320)
    ax0.set_title('cntrl'+str(ch7_cntrl.mean().values))

    ax1.pcolor(cond1(ch7_obs[0],st,ed), vmin=220, vmax=320)
    ax1.set_title('obs')

    ax2.pcolor(cond1(ch7_sbt[0],st,ed), vmin=220, vmax=320)
    ax2.set_title('sbt')

    profile=xr.open_dataset("profiles/profile2d4_2019_dorain_gfs.nc")

#    ax[1].plot(hist(ch7_obs)[0][1:],hist(ch7_obs)[1][0:],label="obs")
#    ax[2].plot(hist(ch7_sbt)[0][1:],hist(ch7_sbt)[1][0:],label="upp")
#
#    legend=ax.legend(loc='best',shadow=False,fontsize='x-small')
#    legend.get_frame().set_facecolor('C0')
#    plt.legend()
#    fig.tight_layout()
#    plt.savefig(figname)
    plt.show()




from pathlib import Path

cases="profile2d4_2019_dorain_gfs_output_cntrl.nc"
figname=cases+'.newq.png' #"test.png"
plot(cases,figname)
