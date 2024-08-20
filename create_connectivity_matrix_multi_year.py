# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 17:19:23 2022

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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from shapely.geometry import Polygon, Point
import shapely
from timeit import default_timer as timer
import os
import seaborn as sns
import pandas as pd

#%%
'''Set Parameters'''
dep_arr= ["Seabed"] #['Surface','Mid-water','Seabed']
year_arr=['2005']#['2005','2006','2007','2008','2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019']  
#tot_arr=[3509] number of particles per release - line 68 does automatically
#tot_par=[147620] #total number of particles over all releases
#season_arr= ['Spring','Summer', 'Winter', 'Fall']# different seasons
season_arr= ['April'] #'July']#'February','March','April','May','June','July','August','September','October','November','December']

cell_size= 0.1

particle_space = '0.01' #unit:degree
dt='10'  #time step minutes
day1='1' #output results every day
day2='60' #final duration
durations = np.array([61, 81, 101]) #PLD duration
comp = np.array([30, 50, 70]) #competency period - need real estimate
trial = 'Prim seed MT Prelim'
repeatdt = '1'
releases = int(durations)/int(repeatdt)
vertvel = 'Vel.0001'
#kh_arr = np.array([0,50,100,150,200,250,300])
kh = 52.1 #52.1 for scotian slope
'''Set Parameters'''

#%%

""" Identify the particles in/out polygon  """     
def  inpolygon(xq, yq, xv, yv):
    shape = xq.shape
    xq = xq.reshape(-1)
    yq = yq.reshape(-1)
    xv = xv.reshape(-1)
    yv = yv.reshape(-1)
    q = [(xq[i], yq[i]) for i in range(xq.shape[0])]
    p = path.Path([(xv[i], yv[i]) for i in range(xv.shape[0])])
    return p.contains_points(q).reshape(shape)                                           
""" Identify the particles in/out polygon  """

''''import grouped occurrences'''''
grouped_occur = np.loadtxt('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Particle_Info/'+trial+'.txt')
num_index = np.loadtxt('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Particle_Info/'+trial+'.txt')
areas = len(num_index)
tot = int(num_index[-1]) #*(int(durations)/int(repeatdt))
''''import grouped occurrences'''''
#%%

for idx_y in range(len(year_arr)):
    year=year_arr[idx_y]       
    print("Year is %s"%year)
    for id_s in range(len(season_arr)):    
        season = season_arr[id_s]
        start=timer()
        print("Season is %s"%season)
        for idx_p in range(len(dep_arr)): 
            dep=dep_arr[idx_p]
    
            ''''import trajectory resutls'''''
            dataset1 = Dataset('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Particle_Info/'+trial+'.txt')
            lat= dataset1.variables['lat'][:,0:61]
            lon= dataset1.variables['lon'][:,0:61]
            depth= dataset1.variables['z'][:,0:61]
            distance = dataset1.variables['distance'][:,0:61]
            time_var = dataset1.variables['time'][:,0:61]
            

            #dataset2 = Dataset("C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Output netCDF\\04082022\\"+trial+' '+dep+"_"+year+season+day2+"days_"+day1+"days_"+repeatdt+"repeatdt_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(tot)+"_B"+str(kh)+".nc")
            #lat_novert = dataset2.variables['lat'][:,0:61]
            #lon_novert = dataset2.variables['lon'][:,0:61]
            #depth_novert = dataset2.variables['z'][:,0:61]
            #distance_novert = dataset2.variables['distance'][:,0:61]

            count_crossing = np.zeros([areas,areas],dtype=np.int32)
            count_retention = np.zeros([areas,1],dtype=np.int32)
            ending=[] #ending area for particle i
            source=[] #source area for particle i
            particle_index = [] # particle index for particle i
            time=[] #time step particle i is in area M
            time_arrival=[]
            time_out=[]
            lat_arr = [] #lat when particle enters a source location
            lon_arr = [] #lon when particle enters a source location
            depth_dist = [] #depth when particle enters a source location

            for idx_1 in range(int(releases)):
                idx_2 = durations[0]       
                for M in range(areas): 
                    if M==0:
                        begin=0 + idx_1*num_index[areas-1] 
                        final=num_index[0] + idx_1*num_index[areas-1] 
    
                    elif M>0:
                        begin = num_index[M-1] + idx_1*num_index[areas-1]
                        final=num_index[M] + idx_1*num_index[areas-1]
                    for i in range(int(begin),int(final)):
                        for t in range(int(comp),idx_2): #durations - estimate with competency period of 'comp' days
                            xq = lon[i,t]
                            yq = lat[i,t]
                            dp = depth[i,t]
                            if str(xq) == '--':
                                break
                            for j in range(len(grouped_occur)):
                                x_occur = grouped_occur[j,1]
                                y_occur = grouped_occur[j,0]
                                xv = np.array([x_occur-cell_size/2,x_occur+cell_size/2,x_occur+cell_size/2,x_occur-cell_size/2])
                                yv = np.array([y_occur-cell_size/2,y_occur-cell_size/2,y_occur+cell_size/2,y_occur+cell_size/2])
                                results = inpolygon(xq, yq, xv, yv)
                                if results==True:
                                   N=j
                                   if M==N:
                                       if lat[i,t] == lat[i,t-1]: #gets rid of particles which don't move
                                           count_retention[M]+=1
                                           time.append(t)
                                           source.append(M)
                                           ending.append(N)
                                           particle_index.append(i)
                                           lat_arr.append(yq)
                                           depth_dist.append(dp)
                                           lon_arr.append(xq)
                                       break
                                   else:
                                       time.append(t)
                                       source.append(M)
                                       ending.append(N)
                                       particle_index.append(i)
                                       lat_arr.append(yq)
                                       depth_dist.append(dp)
                                       lon_arr.append(xq)
                                   break
            lat_arr_in =[]
            lon_arr_in = []
            depth_dist_in = []
            particle_index_in = []
            source_in=[]
            ending_out=[]
            for k in range(len(ending)): ##this block of code is for getting ride of particles which stay in ending locations for more than 1 day
                if k==0: #i.e only counting a particle that reaches an ending location once
                    time_arrival.append(time[k])
                    count_crossing[source[k],ending[k]]+=1 #add 1 to count_crossing array - sources as rows, endings as columns
                    source_in.append(source[k])
                    ending_out.append(ending[k])
                    particle_index_in.append(particle_index[k])
                    lat_arr_in.append(lat_arr[k])
                    lon_arr_in.append(lon_arr[k])
                    depth_dist_in.append(depth_dist[k])
                elif k>0 and k<len(ending)-1:
                    if source[k]==source[k-1]:
                        if particle_index[k]==particle_index[k-1]:
                            if ending[k]==ending[k-1]:
                                continue;
                            else:
                                count_crossing[source[k],ending[k]]+=1
                                time_out.append(time[k-1])
                                time_arrival.append(time[k])
                                source_in.append(source[k])
                                ending_out.append(ending[k])
                                particle_index_in.append(particle_index[k])
                                lat_arr_in.append(lat_arr[k])
                                lon_arr_in.append(lon_arr[k])
                                depth_dist_in.append(depth_dist[k])
                        else:
                            time_out.append(time[k-1])
                            time_arrival.append(time[k])
                            count_crossing[source[k],ending[k]]+=1
                            source_in.append(source[k])
                            ending_out.append(ending[k])
                            particle_index_in.append(particle_index[k])
                            lat_arr_in.append(lat_arr[k])
                            lon_arr_in.append(lon_arr[k])
                            depth_dist_in.append(depth_dist[k])
                    else:
                        time_out.append(time[k-1])
                        time_arrival.append(time[k])
                        count_crossing[source[k],ending[k]]+=1
                        source_in.append(source[k])
                        ending_out.append(ending[k])
                        particle_index_in.append(particle_index[k])
                        lat_arr_in.append(lat_arr[k])
                        lon_arr_in.append(lon_arr[k])
                        depth_dist_in.append(depth_dist[k])
                else: #not consider no particles of areas
                    time_out.append(time[k])    
            num = np.zeros([areas,1],dtype=np.int32) 
            time_arrival_new = []
            time_out_new = [] 
            source_in = np.array(source_in) 
            ending_out = np.array(ending_out) 
            num[0]+=np.sum(count_crossing[0,:]) #adds up first row of count_crossing
            #make index correspond to the count_crossing
            
# =============================================================================
#                 for M in range(areas):
#                     if M>0:
#                         num[M]+=np.sum(count_crossing[M,:])+num[M-1] #adds up all rows of count_crossing, adds rows together in continuous sum
#                     for N in range(areas):
#                         if count_crossing[M,N]!=0:
#                             if M==0:
#                                 for i in range(int(num[0])): #somehow reorders the values in time_arrival and time_out
#                                     if source_in[i]==M and ending_out[i]==N:
#                                         time_arrival_new.append(time_arrival[i])
#                                         time_out_new.append(time_out[i])
#                             else:
#                                 for i in range(int(num[M-1]),int(num[M])):
#                                     if source_in[i]==M and ending_out[i]==N:
#                                         time_arrival_new.append(time_arrival[i])
#                                         time_out_new.append(time_out[i])
#                 #make index correspond to the count_crossing  
#                 time_arrival_new=np.array(time_arrival_new) #unsure what the _new values are referring too
#                 time_out_new=np.array(time_out_new)      
# =============================================================================
            
            
            np.savetxt('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Output/'+trial+dep+"_result_count_retentionno"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+vertvel+".txt", count_retention)
            np.savetxt('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Output/'+trial+dep+"_result_count_crossingno"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+vertvel+".txt", count_crossing)
# =============================================================================
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_lat_arr_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", lat_arr_in)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_lon_arr_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", lon_arr_in)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_time_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", time_arrival)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_depth_dist_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", depth_dist_in)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_particle_id_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", particle_index_in)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_source_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", source_in)
#                 np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_ending_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", ending_out)
# =============================================================================

#%%


                
                
            #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_result_arrival_timeno"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", time_arrival_new)
            #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Connectivity Matrices\\'+trial+dep+"_result_out_timeno_"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", time_out_new)
                
            depth_dist_in = np.array(depth_dist_in)
            lon_arr_in = np.array(lon_arr_in)
            lat_arr_in = np.array(lat_arr_in)
            time_arrival = np.array(time_arrival)
            particle_index_in= np.array(particle_index_in)
            source_in = np.array(source_in)
            ending_out = np.array(ending_out)
            number = int(len(source_in))
            particle_arr_data = np.zeros([number,7])#[particle_index,lat_arr,lon_arr,depth_dist,source,ending])
            for i in range(7):
                for j in range(number):
                    if i == 0:
                        particle_arr_data[j,i] = particle_index_in[j]
                    if i == 1:
                        particle_arr_data[j,i] = lat_arr_in[j]
                    if i == 2:
                        particle_arr_data[j,i] = lon_arr_in[j]
                    if i == 3:
                        particle_arr_data[j,i] = depth_dist_in[j]
                    if i == 4:
                        particle_arr_data[j,i] = source_in[j]
                    if i == 5:
                        particle_arr_data[j,i] = ending_out[j]
                    if i == 6:
                        particle_arr_data[j,i] = time_arrival[j]
                        
            
            np.savetxt('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Output/'+trial+'particle data'+year+season+dep+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeatdt+"dt"+"_B"+str(kh)+".txt", particle_arr_data)    
            end = timer()

            #print("Time to run is "+str(end-start)+' seconds')
