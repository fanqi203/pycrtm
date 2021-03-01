#import ipdb
#ipdb.set_trace()
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np

def preprocess(ds):
    import os
    ds=ds.expand_dims("cases")
    filename=ds.encoding["source"]
    fn=os.path.basename(filename)
    ds=ds.assign_coords(cases=("cases",[fn[39:len(fn)-3]]))
    ds=ds.rename({"__xarray_dataarray_variable__":"bt"})
    return ds

def hist(ch7):
    ch7flat=ch7.where(ch7>0).values.flatten()
    [obshist, obsbin] = np.histogram(ch7flat, density=True, bins=np.arange(190, 320,1))
    avg_BT=ch7flat[np.where(ch7flat>0)].mean()
    labelstr="mean bt "+str(np.round(avg_BT,0))+" for "+str(ch7.cases.values)
    return [obsbin, obshist,  labelstr]

def plot(cases, figname):
    
    f2=cases #"profiles/profile2d4_2019_dorain_gfs_output_test_angle*.nc"
    xall=xr.open_mfdataset(f2,preprocess=preprocess,concat_dim="cases")
    ch7=xall.bt

    f0="profiles/profile2d4_2019_dorain_gfs_output_test_cntrl*.nc"
    x=xr.open_mfdataset(f0,preprocess=preprocess,concat_dim="cases")
    ch7_cntrl=x.bt

    fig, ax = plt.subplots(1, 1)

    for i in np.arange(ch7.shape[0]):
        ax.plot(hist(ch7[i])[0][1:],hist(ch7[i])[1][0:],label=hist(ch7[i])[2]) 
    
        ax.plot(hist(ch7_cntrl)[0][1:],hist(ch7_cntrl)[1][0:],label=hist(ch7_cntrl)[2])

        plt.legend()
        plt.savefig(figname)
        #show()


