# run ipdb debug:

#import ipdb
#ipdb.set_trace()

import xarray as xr
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('legend', fontsize=9, loc='upper left')
import numpy as np

def preprocess(ds):
    import os
    ds=ds.expand_dims("cases")
    filename=ds.encoding["source"]
    fn=os.path.basename(filename)
    ds=ds.assign_coords(cases=("cases",[fn[39:len(fn)-3]]))
    ds=ds.rename({"__xarray_dataarray_variable__":"bt"})
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
    [obshist, obsbin] = np.histogram(ch7flat, density=True, bins=np.arange(180, 320, 3))
    avg_BT=ch7flat[np.where(ch7flat>0)].mean()
    mean1=ch7flat[np.where(ch7flat<240)].mean()
    mean1=np.round(mean1,1)
    mean2=ch7flat[np.where((ch7flat>240) & (ch7flat<280))][0].mean()
#    mean3=ch7flat[np.where(280<ch7flat<320)].mean()
    mean2=np.round(mean2,1)
    mean3=ch7flat[np.where(ch7flat>280)].mean()
    mean3=np.round(mean3,1)
    
    labelstr=str(np.round(avg_BT,0))+' '+str(ch7.cases.values)+' T1='+str(mean1)+' mean2='+str(mean2)+' mean3='+str(mean3)
    
    return [obsbin, obshist,  labelstr]




def plot(cases, figname):
    
    f2=cases #"profiles/profile2d4_2019_dorain_gfs_output_test_angle*.nc"
    xall=xr.open_mfdataset(f2,preprocess=preprocess,concat_dim="cases")
    ch7=xall.bt

    f0="profile2d4_2019_dorain_gfs_output_test_control.nc"
    x=xr.open_mfdataset(f0,preprocess=preprocess,concat_dim="cases")
    ch7_cntrl=0.5*(x.bt[0,7]+x.bt[0,7])

    fobs="profiles/obs_dorian*.nc"
    obs=xr.open_mfdataset(fobs,preprocess=preprocessobs,concat_dim="cases")
    ch7_obs=obs.obs


    fsbt="profiles/sbt_dorian_regrid.nc"
    sbt=xr.open_mfdataset(fsbt,preprocess=preprocessobs,concat_dim="cases")
    ch7_sbt=sbt.SBT124_P0_L8_GLC0


    fig, ax = plt.subplots(1, 1)

    for i in np.arange(ch7.shape[0]):
        ax.plot(hist(ch7[i,7])[0][1:],hist(ch7[i,7])[1][0:],label=hist(ch7[i,7])[2]) 
        


    ax.plot(hist(ch7_cntrl)[0][1:],hist(ch7_cntrl)[1][0:],label=hist(ch7_cntrl)[2],linewidth=3.0,c='black')

#    ax.plot(hist(ch7_obs)[0][1:],hist(ch7_obs)[1][0:],label="obs")
#    ax.plot(hist(ch7_sbt)[0][1:],hist(ch7_sbt)[1][0:],label="upp")

    legend=ax.legend(loc='best',shadow=False)
#fontsize='x-small',
    legend.get_frame().set_facecolor('C0')
    plt.legend()
    #plt.savefig(figname)
    plt.show()



cases="profile2d4_2019_dorain_gfs_output_test_R_ice*.nc" 
figname="test_t_sfc_min.png"
plot(cases,figname)
