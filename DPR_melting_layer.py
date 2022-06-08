# -*- coding: utf-8 -*-
"""
@author: Arnaud Reboud (IGE)
"""

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import proplot
import numpy as np
import numpy.ma as ma
import h5py
import pandas as pd
import re

import funcDPR

Short_name_sat = "Alps"  #Define folder station name
name_station = '_'.join(Short_name_sat.split(' '))
path_dir = "./Data/"+name_station+"/"
hdf_list = [path_dir + x for x in os.listdir(path_dir)]
path_plot = "./Plots/"
path_dir_out = "./Outputs/"

#Stations data
Stations = pd.read_csv('./Data/Stations_metadata.txt',
                       sep=';', index_col=0)

hdffile='C:/Users/rebouda/Documents/02-Data/HDF_data/Alps/2A/DPR/2016-2022/2A-CS-4E47N7E43N-39.GPM.DPR.V9-20211125.20211031-S232814-E232924.043609.V07A.HDF5'
print(hdf_list)

#Some constants and parameters
lat_loc = 45.2 #Grenoble 45.2  #Huaraz -9.52 #Stuttgart 48.78 #Clermont 45.78
lon_loc = 5.76 #Grenoble 5.76 #Huaraz -77.53 #Stuttgart 9.17 #Clermont 3.16
thr_dist = 0.5 # the radius of the circle where the value are considered in degree

geoid_Offs = -51 #m difference between Ellipsoid elevation and geoid elevation, -51m at Grenoble

#The dataframe where to save the melting layer elevations and statistics
BB_GPM = pd.DataFrame(columns=['Time_sat','heightBB_med','topBB_med','botBB_med','#pxl','heightBB_std','heightClutterFree_mean','heightClutterFree_med','heightClutterFree_std'])

### Loop over the whole list of HDF files
for hdffile in hdf_list:
    print(hdffile)
  
    with h5py.File(hdffile, mode='r') as f:
        list_keys = list(f.keys())
        
        cat='FS'
        data_elev = f['/'+cat+'/PRE/elevation'][:] + geoid_Offs
        height_sat = f['/'+cat+'/PRE/height'][:] + geoid_Offs
        heightBB = f['/'+cat+'/CSF/heightBB'][:]
        binBBpeak = f['/'+cat+'/CSF/binBBPeak'][:]
        topBB = f['/'+cat+'/CSF/binBBTop'][:]
        botBB = f['/'+cat+'/CSF/binBBBottom'][:]
        #binDFRmMLBottom = f['/'+cat+'/CSF/binDFRmMLBottom'][:]
        #binDFRmMLTop = f['/'+cat+'/CSF/binDFRmMLTop'][:]
        widthBB = f['/'+cat+'/CSF/widthBB'][:]
        flagBB = f['/'+cat+'/CSF/flagBB'][:]
        flagSnow = f['/'+cat+'/Experimental/flagSurfaceSnowfall'][:]
        heightZeroDeg = f['/'+cat+'/VER/heightZeroDeg'][:] + geoid_Offs
        binZeroDeg = f['/'+cat+'/VER/binZeroDeg'][:]
        airTemperature = f['/'+cat+'/VER/airTemperature'][:]
        phaseNS = f['/'+cat+'/SLV/phaseNearSurface'][:]
        phase = f['/'+cat+'/DSD/phase'][:]
        localZA = f['/FS/PRE/localZenithAngle'][:]
        ellipsoidBinOffset = f['/FS/PRE/ellipsoidBinOffset'][:]
        binEchoBottom = f['/FS/SLV/binEchoBottom'][:]
        binClutterFreeBottom = f['/FS/PRE/binClutterFreeBottom'][:]
        binRealSurface = f['/FS/PRE/binRealSurface'][:]
        zFactorFinal = f['/FS/SLV/zFactorFinal'][:]
        zFactorMeasured = f['/FS/PRE/zFactorMeasured'][:]
        
        data_elev = funcDPR.setnan(data_elev)
        phaseNS = funcDPR.setnan(phaseNS)
        heightBB = funcDPR.setnan(heightBB)
        topBB = funcDPR.setnan(topBB)
        botBB = funcDPR.setnan(botBB)
        binBBpeak = funcDPR.setnan(binBBpeak)
        #binDFRmMLBottom = funcDPR.setnan(binDFRmMLBottom)
        #binDFRmMLTop = funcDPR.setnan(binDFRmMLTop)
        widthBB = funcDPR.setnan(widthBB)
        flagBB= funcDPR.setnan(flagBB)
        heightZeroDeg = funcDPR.setnan(heightZeroDeg)
        airTemperature = funcDPR.setnan(airTemperature)
        phase = funcDPR.setnan(phase)
        localZA = funcDPR.setnan(localZA)
        ellipsoidBinOffset = funcDPR.setnan(ellipsoidBinOffset)
        binEchoBottom = funcDPR.setnan(binEchoBottom)
        binClutterFreeBottom = funcDPR.setnan(binClutterFreeBottom)
        height_sat = funcDPR.setnan(height_sat)
        zFactorFinal = funcDPR.setnan(zFactorFinal)
        binRealSurface = funcDPR.setnan(binRealSurface)
        zFactorMeasured = funcDPR.setnan(zFactorMeasured)
        heightBB.name= 'BB height'
        botBB.name = 'BB Bottom'
        topBB.name = 'BB Top'
        
        #mask the bright band (BB) values where the BB is not detected (flagBB<=0)
        heightBB = np.ma.masked_where(flagBB==0, heightBB) + geoid_Offs
        botBB = np.ma.masked_where(flagBB==0, botBB)
        topBB = np.ma.masked_where(flagBB==0, topBB)
        
        # Get the geolocation data.
        latitude = f[cat+'/Latitude'][:]
        longitude = f[cat+'/Longitude'][:]
        #compute the distance between the city of Grenoble and the radius extremities
        dist = np.sqrt((latitude - lat_loc)**2+(longitude - lon_loc)**2)
        
        #retrieve info about the overpass
        date = re.search(pattern = r"-((\d+).(\d+))-", 
                         string=os.path.basename(hdffile)).group(3)
        hour = re.search(pattern = r"-S(\d+)-",
                         string=os.path.basename(hdffile)).group(1)
        band = re.search(pattern = r"GPM.(\w+).", 
                         string=os.path.basename(hdffile)).group(1)
        
        if (np.sum(phaseNS.mask==False) >0):
            #plot the BB height and the Phase Near Surface
            funcDPR.proplot_map_var(data1= funcDPR.bin_to_height(topBB,localZA,ellipsoidBinOffset,
                                                         geoid_Offset=geoid_Offs,band=band),
                            data2=phaseNS,data_elev=data_elev,
                            longitude=longitude, latitude=latitude,
                            lon_loc=lon_loc,lat_loc=lat_loc,radius=thr_dist,
                            Stations=Stations, name_station=name_station, 
                            date=date, hour=hour, band=band, cat=cat, 
                            path_dir_out=path_plot)
            #plot the BB height, top, and bottom
            funcDPR.proplot_map_BB(data1= heightBB,
                            data2=funcDPR.bin_to_height(topBB,localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band),
                            data3=funcDPR.bin_to_height(botBB,localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band),
                            data_elev=data_elev,
                            longitude=longitude, latitude=latitude,
                            lon_loc=lon_loc,lat_loc=lat_loc, radius=thr_dist,
                            Stations=Stations, name_station=name_station, 
                            date=date, hour=hour, band=band, cat=cat, 
                            path_dir_out=path_plot)
        else:
            print('no phaseNS values to plot in file: '+date+band+cat)
        
        #plot the detection of snowfall by the DPR in the vicinity of the stations
        funcDPR.proplot_map_snow(data= flagSnow, data_elev=data_elev, longitude=longitude, 
                        latitude=latitude, Stations=Stations, name_station=name_station, 
                        date=date, hour=hour, band=band, cat=cat, 
                        path_dir_out=path_plot)
        
        
        #%Compute the median and statistics of variables inside the radius (if there is at least 1 pxl of BB inside the circle)
        flagBB_mask = ma.masked_where(flagBB<=0,flagBB)
        if funcDPR.data_circle(flagBB_mask,latitude,longitude,lon_loc,lat_loc,thr_dist).count()>0 :
            heightBB_stats = {'counts': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightBB),latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                              'mean': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightBB),latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                              'variance': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightBB),latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                              'std': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightBB),latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                              'median': [ma.median(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightBB),latitude,longitude,lon_loc,lat_loc,thr_dist))]
                              }
            heightBB_stats = pd.DataFrame(heightBB_stats)
            
            topBB_stats = {'counts': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, topBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).count()],
                           'mean': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, topBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).mean()],
                           'variance': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, topBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).var()],
                           'std': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, topBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).std()],
                           'median': [ma.median(funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, topBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band))]
                           }
            topBB_stats = pd.DataFrame(topBB_stats)
            
            botBB_stats = {'counts': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, botBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).count()],
                              'mean': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, botBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).mean()],
                              'variance': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, botBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).var()],
                              'std': [funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, botBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).std()],
                              'median': [ma.median(funcDPR.bin_to_height(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, botBB),latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band))]
                              }
            botBB_stats = pd.DataFrame(botBB_stats)
            
            heightClutterFree_stats = {'counts': [funcDPR.bin_to_height(funcDPR.data_circle(binClutterFreeBottom,latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).count()],
                           'mean': [funcDPR.bin_to_height(funcDPR.data_circle(binClutterFreeBottom,latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).mean()],
                           'variance': [funcDPR.bin_to_height(funcDPR.data_circle(binClutterFreeBottom,latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).var()],
                           'std': [funcDPR.bin_to_height(funcDPR.data_circle(binClutterFreeBottom,latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band).std()],
                           'median': [ma.median(funcDPR.bin_to_height(funcDPR.data_circle(binClutterFreeBottom,latitude,longitude,lon_loc,lat_loc,thr_dist),localZA,ellipsoidBinOffset,geoid_Offset=geoid_Offs,band=band))]
                           }
            heightClutterFree_stats = pd.DataFrame(heightClutterFree_stats)
            
            #height of the isotherm 0 (only where the BB is detected)
            heightZeroDeg_stats = {'counts': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightZeroDeg),latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                                   'mean': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightZeroDeg),latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                   'variance': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightZeroDeg),latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                                   'std': [funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightZeroDeg),latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                                   'median': [ma.median(funcDPR.data_circle(ma.masked_where(flagBB_mask.mask, heightZeroDeg),latitude,longitude,lon_loc,lat_loc,thr_dist))]
                                   }
            heightZeroDeg_stats = pd.DataFrame(heightZeroDeg_stats)
            
            BB_GPM = BB_GPM.append(pd.DataFrame([[ date+' '+hour, heightBB_stats['median'][0], topBB_stats['median'][0], 
                                                  botBB_stats['median'][0], heightBB_stats['counts'][0], heightBB_stats['std'][0], 
                                                  heightClutterFree_stats['mean'][0], heightClutterFree_stats['median'][0], heightClutterFree_stats['std'][0] ]],
                                                columns=['Time_sat','heightBB_med','topBB_med','botBB_med','#pxl','heightBB_std',
                                                         'heightClutterFree_mean','heightClutterFree_med','heightClutterFree_std']), 
                                   ignore_index=True)
        else:
            print('No data of BB falling into the circled area')
            heightBB_stats = np.nan
        
        #calculate the phase and temperature (3-dimension variables) inside the rdius
        if funcDPR.data_circle(phase[:,:,0],latitude, longitude, lon_loc, lat_loc, thr_dist).count() > 0:
            height_phase = ma.masked_where(phase.mask, height_sat)
            ground_elev_phase = funcDPR.data_circle(ma.masked_where(phase[:,:,0].mask,data_elev),latitude,longitude,lon_loc, lat_loc, thr_dist)
            
            bins=0
            phase_dict_stats = {'counts': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                                'mean': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                'variance': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                                'std': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                                'median': [ma.median(funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                'mean_height': [funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                'median_height': [ma.median(funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                'std_height': [funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()]
                                }
            phase_stats = pd.DataFrame(phase_dict_stats)
            

            for bins in range(1,phase.shape[2]):
                phase_dict_stats = {'counts': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                                       'mean': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                       'variance': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                                       'std': [funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                                       'median': [ma.median(funcDPR.data_circle(phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                       'mean_height': [funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                       'median_height': [ma.median(funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                       'std_height': [funcDPR.data_circle(height_phase[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()]
                                       }
                phase_stats = phase_stats.append(pd.DataFrame(phase_dict_stats),ignore_index=True)
            
            phase_stats.index.name = 'binRange'
        else:
            print('No phase data falling into the defined area')
            phase_stats = np.nan

        bins=0
        airTemperature_dict_stats = {'counts': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                            'mean': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                            'variance': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                            'std': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                            'median': [ma.median(funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                            'mean_height': [funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                            'median_height': [ma.median(funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                            'std_height': [funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()]
                            }
        airTemperature_stats = pd.DataFrame(airTemperature_dict_stats)
        
        for bins in range(1,airTemperature.shape[2]):
            airTemperature_dict_stats = {'counts': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).count()],
                                   'mean': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                   'variance': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).var()],
                                   'std': [funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()],
                                   'median': [ma.median(funcDPR.data_circle(airTemperature[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                   'mean_height': [funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).mean()],
                                   'median_height': [ma.median(funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist))],
                                   'std_height': [funcDPR.data_circle(height_sat[:,:,bins],latitude,longitude,lon_loc,lat_loc,thr_dist).std()]
                                   }
            airTemperature_stats = airTemperature_stats.append(pd.DataFrame(airTemperature_dict_stats),ignore_index=True)
        
        airTemperature_stats.index.name = 'binRange'
        
        #plot verical profile of temperature with the BB elevations
        plt.figure(figsize=(3,4))
        plt.axvline(x=100, c='k')
        plt.axvline(x=200, c='k')
        plt.axvline(x=0, c='k', linestyle='--')
        
        if np.isscalar(heightBB_stats) == False :
            plt.axhline(y=topBB_stats['median'][0], label='BB_Top', c='orange')
            plt.fill_between(x=np.array([-50,300]), y1=topBB_stats['median'][0] - topBB_stats['std'],
                              y2=topBB_stats['median'] + topBB_stats['std'], alpha=0.1,
                              color='orange')
            plt.axhline(botBB_stats['median'][0], label='BB_Botom', c='green')
            plt.fill_between(x=np.array([-50,300]), y1=botBB_stats['median'][0] - botBB_stats['std'],
                              y2=botBB_stats['median'] + botBB_stats['std'], alpha=0.1,
                              color='green')
            plt.axhline(heightBB_stats['median'][0], label='BB_height', c='tab:blue')
            plt.fill_between(x=np.array([-50,300]), y1=heightBB_stats['median'][0] - heightBB_stats['std'],
                              y2=heightBB_stats['median'] + heightBB_stats['std'], alpha=0.1,
                              color='tab:blue')
            plt.axhline(heightZeroDeg_stats['median'][0], label='heightZeroDeg', c='pink')
            plt.fill_between(x=np.array([-50,300]), y1=heightZeroDeg_stats['median'][0] - heightZeroDeg_stats['std'],
                              y2=heightZeroDeg_stats['median'] + heightZeroDeg_stats['std'], alpha=0.1,
                              color='pink')
           
            plt.plot(airTemperature_stats['median']-273.15,airTemperature_stats.mean_height, 
                     'x-',ms=3, c='r', label='Air Temperature')
           
            plt.ylim(0,3000)
            plt.xlim(-15,15)
            plt.title('Vertical profile '+band+'_'+cat+
                      '\nDate= '+date+' '+hour+'\n Loc= Grenoble (+'+str(thr_dist)+
                      '° radius)')
            plt.xlabel('Temperature [°C]')
            plt.ylabel('height [m.a.s.l]')
            plt.grid(which='both',ls='--')
            plt.legend()
        
            if (np.isscalar(heightBB_stats) | np.isscalar(phase_stats)) == False:
                plt.text(0, 0, '#values= '+str(heightBB_stats['counts'][0]),
                         #+'\n#phase= '+str(phase_stats['counts'][0]),
                         verticalalignment='bottom', horizontalalignment='right',
                         color='black', fontsize=10)
                plt.legend(loc='center left', bbox_to_anchor=(0, 0.4))
                plt.tight_layout()
            elif np.isscalar(heightBB_stats) == False:
                plt.text(0, 0, '#heightBB= '+str(heightBB_stats['counts'][0]),
                         verticalalignment='bottom', horizontalalignment='center',
                         color='black', fontsize=6)
                plt.legend()
                plt.tight_layout()
            elif np.isscalar(phase_stats) == False:
                plt.text(0, 0, '#phase= '+str(phase_stats['counts'][0]),
                         verticalalignment='bottom', horizontalalignment='center',
                         color='black', fontsize=6)
                plt.legend()
                plt.tight_layout()

            plt.savefig(path_plot+'profile_'+date+'_'+hour+'_rad'+str(thr_dist)+'.png', format='png')
        
        #mpl.rcParams.update(mpl.rcParamsDefault)



#export the results of the melting layer elevations to a csv file
BB_GPM.Time_sat = pd.to_datetime(BB_GPM.Time_sat)
BB_GPM.to_csv(path_dir_out+name_station+'_BB_GPM_'+band+'_rad'+str(int(thr_dist*100))+'km.csv',
              sep=';', index=False)
print(BB_GPM)
