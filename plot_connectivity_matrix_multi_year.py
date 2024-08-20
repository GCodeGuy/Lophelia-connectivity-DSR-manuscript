# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:04:00 2022

@author: Graem
"""

from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4_3D, ParticleFile
from parcels import plotTrajectoriesFile, ErrorCode
from glob import glob
import numpy as np
from datetime import timedelta as delta
from matplotlib import path
from parcels import rng as random
import math
from numpy import random as nr
from netCDF4 import Dataset
import  matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
from copy import copy
from timeit import default_timer as timer
#%%
'''Set Parameters'''
year_arr=['2019'] #['2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019']  #water depth
dep_arr=['Seabed'] #['Surface','Mid-water','Seabed']
#tot_arr=[10143,15548,26588]
season_arr= ['January', 'July']#,'February','March','April','May','June','July','August','September','October','November','December']
#season_arr=  ['Spring','Summer','Winter', 'Fall']# different seasons
#season='Winter'# different seasons
cell_size= 0.1
particle_space = '0.01' #unit:degree
dt='10'  #time step minutes
day1='1' #output results every day
day2='60' #final duration
durations = np.array([61])
comp = 30
repeatdt = 0
#kh_arr = np.array([0,50,100,150,200,250,300])
kh = 0 #52.1 for scotian slope
areas = 29
trial = 'CorLocCADoff2001500'

#%%
'''Set Parameters'''

for idx_y in range(len(year_arr)):
    start = timer()
    year=year_arr[idx_y]
    print("Year is %s"%year)
    for idx_d in range(len(dep_arr)):
        dep=dep_arr[idx_d]    
        for idx_s in range(len(season_arr)):
            start = timer()
            season = season_arr[idx_s]
            for idx_1 in range(len(durations)):
                idx_2 = durations[idx_1]
                num_index = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\"+trial+' number for each site ' +str(cell_size)+"cell size "+particle_space+"space.txt")
                num_index = num_index.astype(int)
                count_crossing = np.loadtxt("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\"+trial+dep+"_result_count_crossingno_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+str(repeatdt)+"dt"+"_B"+str(kh)+".txt")
                count_retention = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\Graduate Project\\Parcels Connectivity Matrices\\"+trial+dep+"_result_count_retentionno"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+str(repeatdt)+"dt"+"_B"+str(kh)+".txt")
                for i in range(len(count_crossing)):
                    count_crossing[i,i]=count_retention[i]
                
                    
                for i in range(areas):
                    if i==0:
                        tot=num_index[i]
                    else:
                        tot=(num_index[i]-num_index[i-1])
                    for j in range(areas):
                            #count_crossing[i,j]=count_crossing[i,j]/tot_arr[i]*1
                            count_crossing[i,j]=count_crossing[i,j]/tot*1
                            
                #np.savetxt("D:/Output data/Scotian Shelf/backward3DSS_"+position+"_result_connevtivity_proportion_"+season+day2+"days_"+day1+"hours_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+"_B"+str(kh)+".txt", count_crossing)
                def getnums(s, e,i):
                    return list(range(s, e,i))


                #start, end, interval = 1, 214,1
                
                #names = getnums(start, end, interval)
                #tnames = names.reverse()
                names = ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29')
                tnames = ('29','28','27','26','25','24','23','22','21','20','19','18','17','16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1')
                
                font1 = {'family' : 'Arial',
                         'weight' : 'normal',
                         'size'   : 20} #weight: normal
                
                z015=np.zeros([areas,areas])
                
                for i in range(areas):
                    for j in range(areas):
                        z015[areas-1-j,i]=count_crossing[i,j]
                #z_min, z_max = -np.abs(z).max(), np.abs(z).max()
                fig, ax1 = plt.subplots(figsize = [12,12])
                
                if idx_d==1:
                    dep1='mid-water'
                else:
                    dep1=dep   
                ax1.set_title('%s %s at %s for %s days with %s comp %s cell %s space released at %s dt %s kh'%(season,year,dep1,idx_2-1,comp,cell_size,particle_space,repeatdt, kh),font1)
                plt.tick_params(labelsize=16)
                labels = ax1.get_xticklabels() + ax1.get_yticklabels()
                [label.set_fontname('Arial') for label in labels]
                
                
                   
                #c=ax1.pcolor(x, y, z015, cmap='Blues', vmin=0, vmax=1,edgecolors='k', linestyles = '-', linewidths=0.2)
                z015 = np.ma.masked_where(z015 < 0.0001, z015) #create masked array
                
                palette = copy(plt.cm.YlGnBu)
                
                palette.set_under('w', 0.0001)
                palette.set_bad('w')
                
                im = ax1.imshow(z015, interpolation='none', cmap=palette, vmin=0.0, vmax=1)
                cbar = plt.colorbar(im, extend='both', shrink=0.7, ax=ax1)
                cbar.ax.tick_params(labelsize=14) 
                cbar.set_label('Proportion',family = 'Arial',weight = 'normal',size =20)
                
                
                z015_text = np.round(z015, 2)
                for i in range(areas):
                    for j in range(areas):
                        text = ax1.text(j,i, z015_text[i,j], ha='center', va='center', color='k', size=8)
                   
                
                for i in range(areas):
                    ax1.axhline(i+0.5, lw=1, color='k', zorder=5)
                    ax1.axvline(i+0.5, lw=1, color='k', zorder=5)
                x = np.arange(0,areas,1)
                y = np.arange(0,areas,1)
                
                plt.xticks(x, names)
                plt.yticks(y, tnames)
                
                plt.xlabel("Source Area",font1)
    
                plt.ylabel("Receiving Area",font1)
         
                
                
               # plt.savefig('%s %s at %s for %s days %s competency %d cell %d kh stationary particles removed'%(season,year,dep,idx_2-1,comp,cell_size,kh),bbox_inches='tight',dpi=900)
    
                
                plt.show()
                end = timer()
                print("Time to plot 1 connectivity matrix is "+str(end-start)+' seconds')
