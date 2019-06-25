#!/usr/bin/env python3
import configparser
import h5py 
from pycrtm import pycrtm 
import numpy as np
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as mplcm
from cycler import cycler 
from matplotlib import pyplot as plt

"""
  INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
  INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
  INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
  INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
  INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
  INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics



  INTEGER, PARAMETER :: INVALID_MODEL          = 0
  INTEGER, PARAMETER :: TROPICAL               = 1
  INTEGER, PARAMETER :: MIDLATITUDE_SUMMER     = 2
  INTEGER, PARAMETER :: MIDLATITUDE_WINTER     = 3
  INTEGER, PARAMETER :: SUBARCTIC_SUMMER       = 4
  INTEGER, PARAMETER :: SUBARCTIC_WINTER       = 5
  INTEGER, PARAMETER :: US_STANDARD_ATMOSPHERE = 6

"""
def plotJacobians(chan_list, p, t, q, jacobians, instrument, wavenumbers, ofWhat):
    matplotlib.rc('xtick', labelsize=10) 
    plt.figure()
    NUM_COLORS = len(chan_list)
    cm = plt.get_cmap('brg')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    plt.gca().set_prop_cycle(cycler('color', [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]))
    jacSub = []
    if(ofWhat.lower()=='ozone'):
        for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*q
    elif(ofWhat.lower()=='water'):
         for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*q
    elif(ofWhat.lower()=='temperature'):
          for i in range(jacobians.shape[0]):
            jacobians[i,:] = jacobians[i,:]*t

    for i in range(jacobians.shape[0]):
        if(i+1 in chan_list):
            plt.plot( jacobians[i,:], p )
            jacSub.append(jacobians[i,:])
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    legendList = []
    for i,c in enumerate(chan_list):
        legendList.append('{} {:4.3f} {:4.3f}'.format(c, wavenumbers[c-1], 10000.0/wavenumbers[c-1]))
    plt.legend(legendList, fontsize=6,  ncol=3)
    plt.ylabel('Pressure [hPa]')
    if(ofWhat.lower() =='temperature'): plt.xlabel('Jacobian [K]')
    else: plt.xlabel('Jacobian [K]')
    plt.savefig(instrument+'_jacobians_'+ofWhat+'.png')
    
    if(len(chan_list)>30): matplotlib.rc('xtick', labelsize=6) 
    else: matplotlib.rc('xtick', labelsize=10) 
    plt.close()

    plt.figure()
    if(ofWhat.lower() == 'temperature'): plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel')
    else: plt.title(ofWhat.capitalize()+' Jacobian by Instrument Channel')
    jacSub = np.asarray(jacSub)
    #plt.pcolor(np.arange(wf.shape[0]), myProfiles.P[0,0:99], wf[:,0:99].T, norm = LogNorm( vmin = wf[:,0:99].min(), vmax = 0.3 ))#vmax = wf[:,0:99].max() ) )
    plt.pcolor(np.arange(len(chan_list)+1), p, jacSub[:,:].T ) 
    plt.colorbar()
    plt.ylabel('Pressure [hPa]')
    plt.xlabel('Instrument Channel')
    plt.gca().set_yscale('log')
    plt.gca().invert_yaxis()
    plt.xticks(np.arange(len(chan_list)), chan_list, rotation='vertical')
    plt.yticks(np.array([1000.0, 100.0, 10.0, 1.0, 0.1]),['1000.0','100.0','10.0','1.0','0.1'])
    plt.tight_layout() 
    plt.savefig(instrument+'_jacobians_'+ofWhat+'_pcolor.png')
    plt.close()

def main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m):
    p = {}
    p['US_Std'] ={}
    """    
    with open('Temperature_CRTM1.bin') as f:
        p['US_Std']['temperatureLayers'] = np.fromfile(f,dtype='<f8')
    with open('aerosolConcentration_CRTM1.bin') as f:
        p['US_Std']['aerosolConcentration'] = np.fromfile(f,dtype='<f8')
    with open('aerosolEffectiveRadius_CRTM1.bin') as f:
        p['US_Std']['aerosolEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('cloudEffectiveRadius_CRTM1.bin') as f:
        p['US_Std']['cloudEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('co2_CRTM1.bin') as f:
        p['US_Std']['co2ConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('ozone_CRTM1.bin') as f:
        p['US_Std']['ozoneConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('pressureLevels_CRTM1.bin') as f:
        p['US_Std']['pressureLevels'] = np.fromfile(f,dtype='<f8')
    with open('pressure_CRTM1.bin') as f:
        p['US_Std']['pressureLayers'] = np.fromfile(f,dtype='<f8')
    with open('waterContent_CRTM1.bin') as f:
        p['US_Std']['cloudConcentration'] = np.fromfile(f,dtype='<f8')
    with open('waterVapor_CRTM1.bin') as f:
        p['US_Std']['humidityLayers'] = np.fromfile(f,dtype='<f8')
    """
    with open('Temperature_CRTM2.bin') as f:
        p['US_Std']['temperatureLayers'] = np.fromfile(f,dtype='<f8')
    with open('aerosolConcentration_CRTM2.bin') as f:
        p['US_Std']['aerosolConcentration'] = np.fromfile(f,dtype='<f8')
    with open('aerosolEffectiveRadius_CRTM2.bin') as f:
        p['US_Std']['aerosolEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('cloudEffectiveRadius_CRTM2.bin') as f:
        p['US_Std']['cloudEffectiveRadius'] = np.fromfile(f,dtype='<f8')
    with open('co2_CRTM2.bin') as f:
        p['US_Std']['co2ConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('ozone_CRTM2.bin') as f:
        p['US_Std']['ozoneConcLayers'] = np.fromfile(f,dtype='<f8')
    with open('pressureLevels_CRTM2.bin') as f:
        p['US_Std']['pressureLevels'] = np.fromfile(f,dtype='<f8')
    with open('pressure_CRTM2.bin') as f:
        p['US_Std']['pressureLayers'] = np.fromfile(f,dtype='<f8')
    with open('waterContent_CRTM2.bin') as f:
        p['US_Std']['cloudConcentration'] = np.fromfile(f,dtype='<f8')
    with open('waterVapor_CRTM2.bin') as f:
        p['US_Std']['humidityLayers'] = np.fromfile(f,dtype='<f8')

    p['US_Std']['aerosolType'] = 2 # 1 Dust., 2 Sea salt
    p['US_Std']['cloudType'] = 3 # 1 Water ( I think ) 3, rain
    # Land, water, snow, ice 
    #surfaceTemperatures = np.asarray([272.0, 275.0, 265.0, 269.0])
    surfaceTemperatures = np.asarray([318.0, 0.0, 0.0, 0.0])
    #surfaceFractions = np.asarray([0.1, 0.5, 0.25, 0.15]) 
    surfaceFractions = np.asarray([1.0, 0.0, 0.0, 0.0]) 
    n_absorbers = 2
    climatology = 1
    #LAI = 0.17
    LAI = 0.65
    TUNDRA_SURFACE_TYPE         = 10  # NPOESS Land surface type for IR/VIS Land SfcOptics
    SCRUB_SURFACE_TYPE          =  7  # NPOESS Land surface type for IR/VIS Land SfcOptics
    COARSE_SOIL_TYPE            =  1  # Soil type                for MW land SfcOptics
    GROUNDCOVER_VEGETATION_TYPE =  7  # Vegetation type          for MW Land SfcOptics
    BARE_SOIL_VEGETATION_TYPE   = 11  # Vegetation type          for MW Land SfcOptics
    SEA_WATER_TYPE              =  1  # Water type               for all SfcOptics
    FRESH_SNOW_TYPE             =  2  # NPOESS Snow type         for IR/VIS SfcOptics
    FRESH_ICE_TYPE              =  1  # NPOESS Ice type          for IR/VIS SfcOptics


    landType = SCRUB_SURFACE_TYPE
    soilType = COARSE_SOIL_TYPE
    vegType = BARE_SOIL_VEGETATION_TYPE
    waterType = SEA_WATER_TYPE
    snowType = FRESH_SNOW_TYPE
    iceType = FRESH_ICE_TYPE
    chan_list = [577, 607, 626, 650, 667]
    for k in list(p.keys()):
        N_LAYERS = p[k]['pressureLayers'].shape[0]
        print(surfaceTemperatures,surfaceFractions,LAI,windSpeed10m,windDirection10m, n_absorbers)
        Tb1, Transmission1,\
        emissivity1 = pycrtm.wrap_forward( coefficientPath, sensor_id,\
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan,\
                        p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                        p[k]['co2ConcLayers'],\
                        p[k]['aerosolEffectiveRadius'], p[k]['aerosolConcentration'], p[k]['aerosolType'], \
                        p[k]['cloudEffectiveRadius'], p[k]['cloudConcentration'], p[k]['cloudType'], climatology, \
                        surfaceTemperatures, surfaceFractions, LAI, windSpeed10m, windDirection10m, n_absorbers,\
                        landType, soilType, vegType, waterType, snowType, iceType)



        Tb, Transmission,\
        temperatureJacobian,\
        humidityJacobian,\
        ozoneJacobian, emissivity = pycrtm.wrap_k_matrix( coefficientPath, sensor_id,\
                        zenithAngle, scanAngle, azimuthAngle, solarAngle, nChan,\
                        p[k]['pressureLevels'], p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], p[k]['ozoneConcLayers'],\
                        p[k]['co2ConcLayers'],\
                        p[k]['aerosolEffectiveRadius'], p[k]['aerosolConcentration'], p[k]['aerosolType'],\
                        p[k]['cloudEffectiveRadius'], p[k]['cloudConcentration'], p[k]['cloudType'], climatology, \
                        surfaceTemperatures, surfaceFractions, LAI, windSpeed10m, windDirection10m, n_absorbers,\
                        landType, soilType, vegType, waterType, snowType, iceType)
                     
        print('done crtm bit.')
        h5 = h5py.File('cris_wavenumbers.h5','r')
        wavenumbers = np.asarray(h5['wavenumbers'])
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['ozoneConcLayers'], ozoneJacobian, 'CrIS_'+k+'_', wavenumbers, 'ozone')
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['ozoneConcLayers'], temperatureJacobian, 'CrIS_'+k+'_', wavenumbers, 'temperature')
        plotJacobians(chan_list, p[k]['pressureLayers'], p[k]['temperatureLayers'], p[k]['humidityLayers'], humidityJacobian, 'CrIS_'+k+'_', wavenumbers, 'water')
        with open('bt.bin') as f:
            dataTb = np.fromfile(f, dtype='<f8' )
        plt.figure()
        plt.plot(wavenumbers,Tb1-dataTb)
        plt.savefig('spectrum.png')

        with open('emissivity.bin') as f:
            emissivitySaved = np.fromfile(f, dtype='<f8' )
        plt.figure()
        plt.plot(wavenumbers,emissivity1-emissivitySaved)
        plt.savefig('emissivity.png') 
if __name__ == "__main__":
    pathInfo = configparser.ConfigParser()
    pathInfo.read('crtm.cfg')
    coefficientPath = pathInfo['CRTM']['coeffs_dir']
    #sensor_id = 'atms_n20'
    sensor_id = 'cris_npp'
    zenithAngle = float(30.0)
    scanAngle = float(26.37293341421)
    azimuthAngle = float(0.0)
    solarAngle = float(0.0)
    nChan = int(1305)
    surfaceType = int(1)
    windSpeed10m = float(0.0)
    windDirection10m = float(0.0)
    surfaceTemperature = 0
    # assuming this is zenith? probably need Az too?
    main(coefficientPath, sensor_id,\
        zenithAngle, scanAngle, azimuthAngle, solarAngle,\
        nChan, surfaceType, surfaceTemperature, windSpeed10m, windDirection10m)
 
