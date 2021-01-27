import ipdb
ipdb.set_trace()
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
x=xr.open_dataset("res.nc")
print(x)
data=x.__xarray_dataarray_variable__
ch7=data[7]
ch7flat=ch7.where(ch7>0).values.flatten()
print(ch7flat.shape)
[obshist, obsbin] = np.histogram(ch7flat, density=True, bins=np.arange(190, 320,1))
print(obshist)
print(obsbin)

fig, ax = plt.subplots(1, 1)
avg_BT=ch7flat[np.where(ch7flat>0)].mean()
labelstr="ch7 and the mean is "+str(np.round(avg_BT,0))

ax.plot(obsbin[1:],obshist[0:],label=labelstr) #  and the mean is "+str(np.round(avg_BT,0)))
#_=plt.hist(ch7flat,bins='auto',density=True) # legend()
plt.legend()
plt.show()

#plt.savefig("pdf_"+str(thresh)+".png")
