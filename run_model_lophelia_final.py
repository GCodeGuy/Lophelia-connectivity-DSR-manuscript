# -*- coding: utf-8 -*-


"""
Created on Thu Mar 24 15:39:40 2022

@author: Graem
"""

from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4_3D, ParticleFile, Variable, ParcelsRandom
from parcels import DiffusionUniformKh, ErrorCode
from glob import glob
import numpy as np
from datetime import timedelta as delta
import math
from netCDF4 import Dataset
from operator import attrgetter

'''Set Parameters'''
year_arr=['2005'] #, '2006', '2007', '2008', '2009','2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019']  #trial years
dep_arr= ['Seabed']
season_arr=['January','March','May','July','September','November'] #['Spring','Summer','Fall','Winter']# different seasons
cell_size= '0.1'
particle_space = ['0.01', '0.005', '0.003', '0.0001'] #unit:degree
dt='10'  #time step minutes
day1='1' #output results days
day2='60' #final duration
repeat = '0' #release particles every x days
kh='VARkh'
Area = 'Lophelia_final_run' #'FCBB' 
vertvel='VARvel'


class DistParticle(JITParticle):  # Define a new particle class that contains three extra variables
    distance = Variable('distance', initial=0., dtype=np.float32)  # the distance travelled
    prev_lon = Variable('prev_lon', dtype=np.float32, to_write=False,
                        initial=attrgetter('lon'))  # the previous longitude
    prev_lat = Variable('prev_lat', dtype=np.float32, to_write=False,
                        initial=attrgetter('lat'))  # the previous latitude.
    particle_swim = Variable('particle_swim', initial=0., dtype=np.float32)
    particle_sink = Variable('particle_sink', initial=0., dtype=np.float32)
    
def TotalDistance(particle, fieldset, time):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
    
    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
    particle.prev_lat = particle.lat  
    
def DeleteParticle(particle, fieldset, time):
    particle.delete()


def VertVel(particle, fieldset, time): #add sink and swim speeds
    swim_speed_mean = -0.00072 #m/s #Stromberg and Larssson 2017
    swim_speed_sd = 0.00018 #m/s
    sink_speed_mean = 0.00072 #m/s #Stromberg and Larssson 2017
    sink_speed_sd = 0.00018 #m/s
    
    if time > 0 and time < 86400*30 and particle.depth > 5:
        swim = ParcelsRandom.normalvariate(swim_speed_mean, swim_speed_sd)
        particle.depth += swim * particle.dt
        particle.particle_swim = swim
        
    if time >= 86400*30:
        sink = ParcelsRandom.normalvariate(sink_speed_mean, sink_speed_sd)
        particle.depth += sink * particle.dt
        particle.particle_sink = sink

        
data_path = '/home/gguy/projects/def-metaxas/gguy/ComputeCanada/' 
input_netCDF = 'BNAM/'
domain_resolution = np.loadtxt(data_path +"Particle_Info/Domain_resolution_in_x.txt") #in metres
domain_lat = np.loadtxt(data_path +"Particle_Info/Domain_latp.txt")
for idx_y in range(len(year_arr)):
    year=year_arr[idx_y]
    print("Year is %s"%year)
    for idx_d in range(len(particle_space)):
        dep=dep_arr[0]
        particle_space=particle_space[idx_y]
        print("Particle space is %s"%particle_space)
        for idx_s in range(len(season_arr)):
            season = season_arr[idx_s]
            print("Season is %s"%season)
            ufiles = sorted(glob(data_path+input_netCDF+year+'/'+season+'/'+'3DCurrents_U*.nc'))
            vfiles = sorted(glob(data_path+input_netCDF+year+'/'+season+'/'+'3DCurrents_V*.nc'))
            wfiles = sorted(glob(data_path+input_netCDF+year+'/'+season+'/'+'3DCurrents_W*.nc'))
            mesh_mask =  data_path+input_netCDF+year+'/'+'3DCurrents_mask_'+year+'_BNAM.nc'
            order_1= [2,3,0,1,4,5] #order for Jan,Mar,July,Oct
            order_2 = [2,3,4,5,0,1] #order for Feb,Aug
            order_3 = [0,1,4,5,2,3] #order for April, August
            order_4 = [4,5,2,3,0,1] #order for May, June, Sept, Oct
            order_5 = [2,3,0,1] #order for November
            order_6 = [0,1] #order for Dec
        

            if season=="January" or season=="March" or season=="July":
                ufiles = [ufiles[i] for i in order_1]
                vfiles = [vfiles[i] for i in order_1]
                wfiles = [wfiles[i] for i in order_1]
            elif season=="February":
                ufiles = [ufiles[i] for i in order_2]
                vfiles = [vfiles[i] for i in order_2]
                wfiles = [wfiles[i] for i in order_2]
            elif season=='April' or season=="August":
                ufiles = [ufiles[i] for i in order_3]
                vfiles = [vfiles[i] for i in order_3]
                wfiles = [wfiles[i] for i in order_3]
            elif season=='May'or season=='June' or season=="September"or season=="October":
                ufiles = [ufiles[i] for i in order_4]
                vfiles = [vfiles[i] for i in order_4]
                wfiles = [wfiles[i] for i in order_4]
            elif season == 'November':
                ufiles = [ufiles[i] for i in order_5]
                vfiles = [vfiles[i] for i in order_5]
                wfiles = [wfiles[i] for i in order_5]

     
            
            filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                         'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                         'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles}}
            
            variables = {'U': 'U',
                         'V': 'V',
                         'W': 'W'}
            dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'deptht', 'time': 'time_counter'},
                          'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'deptht', 'time': 'time_counter'},
                          'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'deptht', 'time': 'time_counter'}}

            fieldset = FieldSet.from_nemo(filenames, variables, dimensions, allow_time_extrapolation=True)
        
       
            dataset = Dataset(mesh_mask)
            kh_zonal = dataset.variables['fmask'][0,:,:,:]  
            kh_meridional = dataset.variables['fmask'][0,:,:,:]
                        
            kh_z = (0.0103*(domain_resolution*100)**(4/3))/10000
            kh_m = (0.0103*(domain_resolution*100)**(4/3))/10000
            
            fieldset.add_constant("dres", 0.01)
            fieldset.add_field(Field('Kh_zonal', kh_zonal*kh_z, depth=fieldset.U.grid.depth,lat=fieldset.U.grid.lat, 
                                                  lon=fieldset.U.grid.lon,  mesh='spherical'))
            fieldset.add_field(Field('Kh_meridional', kh_meridional*kh_m, depth=fieldset.U.grid.depth,lat=fieldset.U.grid.lat, 
                                                  lon=fieldset.U.grid.lon,  mesh='spherical')) 

            lonp = np.loadtxt(data_path +"Particle_Info/"+Area+' lonp '+cell_size+"cell size "+particle_space+"space.txt")
            latp = np.loadtxt(data_path +"Particle_Info/"+Area+' latp '+cell_size+"cell size "+particle_space+"space.txt")
            depth = np.loadtxt(data_path +"Particle_Info/"+Area+' depth '+cell_size+"cell size "+particle_space+"space.txt")
             
            depth=-depth  ##CHANGE DEPENDING ON RELEASE DEPTH = /2 for mid-water
            
            
            tot=len(lonp)
            print ('Total particles: ',tot)
            #depth = -np.ones(len(lonp))*3
            repeatdt=delta(days=20)
            pset = ParticleSet.from_list(fieldset=fieldset, pclass=DistParticle, 
                                         lon=lonp,
                                         lat=latp,
                                         depth=depth,
                                        repeatdt=repeatdt)
                                         
            pfile = ParticleFile('/home/gguy/projects/def-metaxas/gguy/ComputeCanada/Output/'+Area+dep+"_"+year+season+day2+"days_"+day1+"days_"+repeat+"repeatdt_"+dt+"minutes_"+cell_size+"cell_"+particle_space+"degree_"+str(tot)+"_B"+kh+vertvel,
                                     pset, outputdt=delta(days=int(day1)))
     
    
            kernels = pset.Kernel(TotalDistance) + DiffusionUniformKh + AdvectionRK4_3D + pset.Kernel(VertVel) ##+pset.Kernel(variable_kh)
            pset.execute(kernels, runtime=delta(days=int(day2)), dt=delta(minutes=int(dt)), recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle},
                         output_file=pfile) 
            pfile.export()
     


