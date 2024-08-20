# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:43:22 2022

@author: Graem
"""

import numpy as np
from matplotlib import path
import math
from netCDF4 import Dataset
import  matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import csv
import itertools
import copy
import statistics
import igraph as ig
import statistics as stats
import openpyxl
#%%
'''Set Parameters'''
seasons_arr = ["Winter","Spring","Summer","Fall"]
#season = "Summer"
year_arr=['B2005','B2006','B2007','B2008','B2009','B2010','B2011','B2012','B2013','B2014','B2015', 'B2016','B2017','B2018'] 
dep_arr=['Seabed']
#month_arr=['February','April'] #,'May','July','September','November']
#month = 'January'# different months
#year= 'B2005'
dep = 'Seabed'
#month_num = 1
particle_space = '0.003'
month_num_arr = [2,4] #,5,7,9,11]
season_num = [1,2,3,4]
cell_size= '0.1'
comp_arr = ['30'] #,'20']
comp = '30'
particle_spaces = '0.003' #, '0.005'] #unit:degree
dt='10'  #time step minutes
day1='1' #output results days
day2='60' #final duration
repeat = '0' #release particles every x days
durations = np.array([61])
idx_2=durations[0]
kh='VARkh'
trial = 'Bdom_Lophelia_final_run' #'FCBB' 
Run = 'Bdom_Loph_linearVELincrease'
vertvel='VARvel'

#%%
seasonal_dict = {
    "Winter": ["January",'February','March'],
    "Spring": ['April','May','June'],
    "Summer": ['July','August','September'],
    "Fall": ['October','November', 'December']}

#%%
year = 'B2005'
year2 = 'B2006'
year3 = 'B2007'
year4 = 'B2008'

season = 'Winter'
season_number = 1
particle_space = '0.003'
comp = '30'

#%%            
"SEASONAL MULTI YEAR AVERAGES"

'''Graph theory analysis for seasonal averages over entire time period'''


coords = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\Final Runs\\"+trial+'.txt')
#coords = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\Final Runs\\Lophelia_GOM_sources_SummarizeWithin.txt", delimiter = ',', skiprows=1)
coords[:,2] = range(len(coords))
coords = pd.DataFrame(coords)
coords = coords.sort_values(by=[1])
coords = np.array(coords)
#np.savetxt('BdomLF_SourceLocations_byLongitude.csv',coords, delimiter=',')

for idx_c in range(len(comp_arr)):
    comp = comp_arr[idx_c]
    print("Comp is %s"%comp)
    for idx_s in range(len(seasons_arr)):
        season = seasons_arr[idx_s]
        #month_num = month_num_arr[idx_s]
        print("Season is %s"%season)
        #for i in range(len(particle_spaces)):
        #particle_space = particle_spaces[i]
        num_index = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\Final Runs\\"+trial+" number for each site "+str(cell_size)+"cell size "+particle_space+"space.txt")
        num_index = num_index.astype(int)
        #print("Particle space is %s" %particle_space)
        ave_count_crossing = np.loadtxt("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\"+Run+'nts_epd'+dep+"_averagecount"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+comp+'comp'+repeat+"dt_B"+kh+vertvel+".txt")
        
        cc_rearrange = np.zeros([len(ave_count_crossing),len(ave_count_crossing)])
        
        for i in range(len(ave_count_crossing)):  ## REarranges node in matrix by longitude from West to East
            for j in range(len(ave_count_crossing)):
                cc_rearrange[i,j] = ave_count_crossing[int(coords[i,2]), int(coords[j,2])]
        
        G = nx.DiGraph(cc_rearrange)
        
        G2 = ig.Graph.Weighted_Adjacency(cc_rearrange, mode = 'directed')
        pagerank = G2.pagerank(directed = True, weights = 'weight')
        
        betweenness = G2.betweenness(directed = True) #use unweighted betweeness


        edgelist = nx.to_edgelist(G)
        
        edgelist = list(edgelist)  #turns edgelist to list of tuples
                  
        probability_matrix = cc_rearrange/((num_index[1]-num_index[0])*3)*100  #probability of arrival - proportion of individuals originating from a donor planning unit which arrive into a recipient planning unit - fration units
        migration_matrix = cc_rearrange/cc_rearrange.sum(axis=0)*100 #percentage of individuals arriving at a recipient unit that originated from a donor planning unit. elements of matrix are relative to the destination
        
        retention=np.diag(cc_rearrange)
        self_recruitment = retention/cc_rearrange.sum(axis=0)*100 #particles released at source location which also recruit at source location/ particles released from all origins which recruit at source location ##correct as of 14/03/2023 
        local_retention = retention/((num_index[1]-num_index[0])*3)*100 #multiplied by 3 because 3 releases per season - Correct as of 23/01/2023
        
        in_degree=np.count_nonzero(cc_rearrange,axis=0) #Correct as of 23/01/2023 #sum all source areas for a destination
        out_degree=np.count_nonzero(cc_rearrange,axis=1) #Correct as of 23/01/2023 #sum all destination areas for a source
        
        node_table = np.zeros([len(G._node),9])
        
        for i in range(len(G._node)):
            node_table[i,0] = int(i)
            node_table[i,1] = coords[i,1]
            node_table[i,2] = coords[i,0]
            node_table[i,3] = int(in_degree[i])
            node_table[i,4] = int(out_degree[i])
            node_table[i,5] = pagerank[i]
            node_table[i,6] = betweenness[i]
            node_table[i,7] = self_recruitment[i]
            node_table[i,8] = local_retention[i]
        

        
        node_table = pd.DataFrame(data = node_table,
                                  columns = ['id','x','y','in_degree','outdegree','pagerank','betweenessC','self_rectruitment','local_retention'])
        
        node_header = ['id','x','y','in_degree','outdegree','pagerank','betweennessC','self_rectruitment','local_retention']
        
        for i, nlrow in node_table.iterrows():
            G._node[nlrow['id']] = nlrow[1:3].to_dict()
        
        
        node_positions = {node[0]: (node[1]['x'], node[1]['y']) for node in G.nodes(data=True)}
        
        connection_edge_table = np.zeros([len(G.edges),9]) #node1, x1, y1, node2,x2, y2 weight
        
        for i in range(len(G.edges)):
            connection_edge_table[i,0]= edgelist[i][0]
            connection_edge_table[i,1]= node_positions[int(connection_edge_table[i][0])][0]
            connection_edge_table[i,2] = node_positions[int(connection_edge_table[i][0])][1]
            connection_edge_table[i,3] = edgelist[i][1]
            connection_edge_table[i,4] = node_positions[int(connection_edge_table[i][3])][0]
            connection_edge_table[i,5] = node_positions[int(connection_edge_table[i][3])][1]
            connection_edge_table[i,6] = edgelist[i][2]['weight']
            connection_edge_table[i,7] = probability_matrix[int(connection_edge_table[i,0]),int(connection_edge_table[i,3])]
            connection_edge_table[i,8] = migration_matrix[int(connection_edge_table[i,0]),int(connection_edge_table[i,3])] 
            
        connection_edge_table = pd.DataFrame(data = connection_edge_table,
                                             columns = ['node1','x1','y1','node2','x2','y2','edge_weight','probability_matrix','migration_matrix'])
        
        
        edge_header = ['node1','x1','y1','node2','x2','y2','edge_weight','probability_matrix','migration_matrix']
        
        #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+trial+dep+"_edge_table_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".txt",connection_edge_table, header = str(edge_header)) 
        #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+trial+dep+"_node_table_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".txt",node_table, header = str(node_header))
        
        
        '''save dataframes with headers'''
        
        connection_edge_table.to_excel('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+'LVELinc'+"_edge_table_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+'.xlsx', header=edge_header)
        node_table.to_excel('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+'LVELinc'+"_node_table_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".xlsx", header=node_header)
        
        #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+trial+dep+"Betweenness_Long_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".txt", betweenness_arr)
        #np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+trial+dep+"BetweennessW_Long_nts_epd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".txt", betweennessWeight_arr)
#%%
"INDIVIDUAL SEASONS"
''' Graph theory analysis for individual seasons multi years'''

spacing_list = []
tot_con_list = []
tot_edge_w_list = []
ave_pro_con_list = []
comp_per_list = []
year_list = []
season_list = []
connection_summary = []
stdev_list = []


coords = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\Final Runs\\"+trial+'.txt')
coords[:,2] = range(len(coords))
coords = pd.DataFrame(coords)
coords = coords.sort_values(by=[1])
coords = np.array(coords)

for idx_c in range(len(comp_arr)):
    comp = comp_arr[idx_c]
    print("Comp is %s"%comp)
    
    for idx_s in range(len(seasons_arr)):
        season = seasons_arr[idx_s]
        season_number = season_num[idx_s]
        print("Season is %s"%season)
        
        for idx_y in range(len(year_arr)):
            year = year_arr[idx_y]
            print("Year is %s"%year)
            
            #for i in range(len(particle_spaces)):
                #particle_space = particle_spaces[i]
            num_index = np.loadtxt("C:\\Users\\Graem\\OneDrive - Dalhousie University\Documents\Dal 2020-2021\\Graduate School\\Graduate Project\\Parcels Particle Seeding\\Final Runs\\"+trial+" number for each site "+str(cell_size)+"cell size "+particle_space+"space.txt")
            num_index = num_index.astype(int)
            print("Particle space is %s" %particle_space)
            count_crossing = np.loadtxt("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\"+Run+'nts_epd'+dep+"_Count_total"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+comp+'comp'+repeat+"dt_B"+kh+vertvel+".txt")
           
            cc_rearrange = np.zeros([len(count_crossing),len(count_crossing)])
            
            for i in range(len(count_crossing)):  ## REarranges node in matrix by longitude from West to East
                for j in range(len(count_crossing)):
                    cc_rearrange[i,j] = count_crossing[int(coords[i,2]), int(coords[j,2])]
            
            G = nx.DiGraph(cc_rearrange)
            
            G2 = ig.Graph.Weighted_Adjacency(cc_rearrange, mode = 'directed')
            pagerank = G2.pagerank(directed = True, weights = 'weight')
            betweenness = G2.betweenness(directed = True) #use unweighted betweeness
        
            edgelist = nx.to_edgelist(G)

            edgelist = list(edgelist)  #turns edgelist to list of tuples
            
            retention=np.diag(cc_rearrange)
            
            
           
            probability_matrix = cc_rearrange/((num_index[1]-num_index[0])*3)*100  #probability of arrival - proportion of individuals originating from a donor planning unit which arrive into a recipient planning unit
            migration_matrix = cc_rearrange/cc_rearrange.sum(axis=0)*100 #proportion of individuals arriving at a recipient unit that originated from a donor planning unit. elements of matrix are relative to the destination
            
            in_degree=np.count_nonzero(cc_rearrange,axis=0)  #number of incomming connecting nodes (irrespective of weight)    
            out_degree=np.count_nonzero(cc_rearrange,axis=1)  # number of outgoing connecting nodes (irrespective of weight)
            retention=np.diag(cc_rearrange)
  
            self_recruitment = retention/cc_rearrange.sum(axis=0)*100 #particles released at source location which also recruit at source location/ particles released from all origins which recruit at source location ##correct as of 14/03/2023 
            local_retention = retention/((num_index[1]-num_index[0])*3)*100 #multiplied by 3 because 3 releases per season - Correct as of 23/01/2023
            
            
            node_table = np.zeros([len(G._node),9])
            
            for i in range(len(G._node)):
                node_table[i,0] = int(i)
                node_table[i,1] = coords[i,1]
                node_table[i,2] = coords[i,0]
                node_table[i,3] = int(in_degree[i])
                node_table[i,4] = int(out_degree[i])
                node_table[i,5] = pagerank[i]
                node_table[i,6] = betweenness[i]
                node_table[i,7] = self_recruitment[i]
                node_table[i,8] = local_retention[i]
            
            
            
            node_table = pd.DataFrame(data = node_table,
                                      columns = ['id','x','y','in_degree','outdegree','pagerank','betweennessC','self_recruitment','local_retention'])
            
            node_header = ['id','x','y','in_degree','outdegree','pagerank','betweennessC','self_recruitment','local_retention']
            
            for i, nlrow in node_table.iterrows():
                G._node[nlrow['id']] = nlrow[1:3].to_dict()
            
        
            node_positions = {node[0]: (node[1]['x'], node[1]['y']) for node in G.nodes(data=True)}
            
            
            connection_edge_table = np.zeros([len(G.edges),9]) #node1, x1, y1, node2,x2, y2 weight
            
            for i in range(len(G.edges)):
                connection_edge_table[i,0]= edgelist[i][0]
                connection_edge_table[i,1]= node_positions[int(connection_edge_table[i][0])][0]
                connection_edge_table[i,2] = node_positions[int(connection_edge_table[i][0])][1]
                connection_edge_table[i,3] = edgelist[i][1]
                connection_edge_table[i,4] = node_positions[int(connection_edge_table[i][3])][0]
                connection_edge_table[i,5] = node_positions[int(connection_edge_table[i][3])][1]
                connection_edge_table[i,6] = edgelist[i][2]['weight']
                connection_edge_table[i,7] = probability_matrix[int(connection_edge_table[i,0]),int(connection_edge_table[i,3])]
                connection_edge_table[i,8] = migration_matrix[int(connection_edge_table[i,0]),int(connection_edge_table[i,3])] 
            
                
            connection_edge_table = pd.DataFrame(data = connection_edge_table,
                                                 columns = ['node1','x1','y1','node2','x2','y2','edge_weight','probability_matrix','migration_matrix'])
            
            edge_header = ['node1','x1','y1','node2','x2','y2','edge_weight','probability_matrix','migration_matrix']
            
            
            probability_average = (sum(connection_edge_table['probability_matrix'])/len(connection_edge_table)) #in percent
            stdev = stats.stdev(connection_edge_table['probability_matrix'])
            total_edge_weight = sum(connection_edge_table['edge_weight'])
            total_connections = len(connection_edge_table)
            year_count = year.replace("B","")
            year_count = int(year_count)
            season_count = season_number
            comp_count = comp
            space_count = particle_space
            
            
            spacing_list.append(space_count)
            tot_con_list.append(total_connections)
            tot_edge_w_list.append(total_edge_weight)
            ave_pro_con_list.append(probability_average)
            comp_per_list.append(comp_count)
            year_list.append(year_count)
            season_list.append(season_count)
            stdev_list.append(stdev)
            
            
            connection_edge_table.to_excel('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+"LVELinc_edge_table_nts_epd"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".xlsx", header=edge_header) 
            node_table.to_excel('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+"LVELinc_node_table_nts_epd"+year+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".xlsx", header=node_header)
#%%                
connection_summary = np.zeros([112,8])#[years,season, particle_spaces,comp_period,tot_conn,tot_edge_weight,ave_probability_connections])
for i in range(8):
    for j in range(112):
        if i == 0:
            connection_summary[j,i] = year_list[j]
        if i == 1:
            connection_summary[j,i] = season_list[j]
        if i == 2:
            connection_summary[j,i] = spacing_list[j]
        if i == 3:
            connection_summary[j,i] = comp_per_list[j]
        if i == 4:
            connection_summary[j,i] = tot_con_list[j]
        if i == 5:
            connection_summary[j,i] = tot_edge_w_list[j]
        if i == 6:
            connection_summary[j,i] = ave_pro_con_list[j]
        if i == 7:
            connection_summary[j,i] = stdev_list[j]

connection_summary = pd.DataFrame(data = connection_summary,
                                     columns = ['year','season','particle_spacing','competency_period','total_connections','total_edge_weight','average_connection_probability','stdev_connection_probability'])            
summary_header = ['year','season','particle_spacing','competency_period','total_connections','total_edge_weight','average_connection_probability','stdev_connection_probability']
   
connection_summary.to_excel('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+'LVELinc'+"_seasonal_connection_data"+year_arr[0]+'-'+year_arr[-1]+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+repeat+"dt"+"_B"+str(kh)+vertvel+".xlsx", header = summary_header) 
#%%
connection_stdev = np.zeros([112,5])
for i in range(5):
    for j in range(112):
        if i == 0:
            connection_stdev[j,i] = year_list[j]
        if i == 1:
            connection_stdev[j,i] = season_list[j]
        if i == 2:
            connection_stdev[j,i] = comp_per_list[j]
        if i == 3:
            connection_stdev[j,i] = ave_pro_con_list[j]
        if i == 4:
            connection_stdev[j,i] = stdev_list[j]
        

connection_stdev = pd.DataFrame(data = connection_stdev,
                                columns = ['year', 'season', 'competency_period', 'ave_con_prob', 'ave_con_stdev'])

headers = ['year', 'season', 'competency_period', 'ave_con_prob', 'ave_con_stdev']

connection_stdev.to_csv('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\'+trial+dep+"_seasonal_ave_conn_stdev"+year_arr[0]+'-'+year_arr[-1]+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+repeat+"dt"+"_B"+str(kh)+vertvel+".csv",sep=',',header = headers) 
#%%
winter_connections_30 = []
summer_connections_30 = []
fall_connections_30 = []
spring_connections_30 = []

for i in range(len(connection_summary)):
    if connection_summary['season'][i]==1 and connection_summary['competency_period'][i] ==30:
        winter_connections_30.append(connection_summary['total_connections'][i])
        
    if connection_summary['season'][i]==2 and connection_summary['competency_period'][i] ==30:
        spring_connections_30.append(connection_summary['total_connections'][i])
    
    if connection_summary['season'][i]==3 and connection_summary['competency_period'][i] ==30:
        summer_connections_30.append(connection_summary['total_connections'][i])
    
    if connection_summary['season'][i]==4 and connection_summary['competency_period'][i] ==30:
        fall_connections_30.append(connection_summary['total_connections'][i])


wi_con_30_stdev = statistics.pstdev(winter_connections_30)
sp_con_30_stdev = statistics.pstdev(spring_connections_30)
su_con_30_stdev = statistics.pstdev(summer_connections_30)
fa_con_30_stdev = statistics.pstdev(fall_connections_30)

wi_con_30_mean = statistics.mean(winter_connections_30)
sp_con_30_mean = statistics.mean(spring_connections_30)
su_con_30_mean = statistics.mean(summer_connections_30)
fa_con_30_mean = statistics.mean(fall_connections_30)



#%%
            '''
            for i in range(len(edge_list)):
                for j in range(3):
                    if j < 2:
                        edge_list[i,j] = edgelist[i][j]
                    else:
                        edge_list[i,j] = edgelist[i][j]['weight']
                        '''
                    
           
            
            '''
            # Add edges and edge attributes
            for i, elrow in edgelist.iterrows():
                G.add_edge(elrow[0], elrow[1], attr_dict=elrow[2:].to_dict())
            '''
            
            
            
         

            # Preview of node_positions with a bit of hack (there is no head/slice method for dictionaries).
            
            plt.figure(figsize=(12, 12))
            nx.draw_networkx_nodes(G, pos=node_positions, node_size=10, node_color='red',alpha = 0.9)
            nx.draw_networkx_edges(G, pos=node_positions, alpha = 0.5, edge_color = 'blue')
            plt.title('NWA Lophelia connections', size=15)
            plt.show()
            
            
          
            

            
''''            np.savetxt('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+'edge_weight_test.txt',edge_weight )
            for i in range(5):
                for j in range(areas):
                    if i == 0:
                        connection_edge_table[j,i] = 
                    if i == 1:
                        connection_edge_table[j,i] = 
                    if i == 2:
                        connection_edge_table[j,i] = 
                    if i == 3:
                        connection_edge_table[j,i] = 
                    if i == 4:
                        connection_edge_table[j,i] = 
                    if i == 5:
                        connection_edge_table[j,i] =
                    if i == 6:
                        connection_edge_table[j,i] = 
                        
                        
                        
            ## len(G.edges)
            ## create new table with 7 columns - source lat, source lat, sink long, sink long, source #, sink #, weight
            ## number of rows = len(G.edges)'''
#%%
# Grab edge list data hosted on Gist
edgelist = pd.read_csv('https://gist.githubusercontent.com/brooksandrew/e570c38bcc72a8d102422f2af836513b/raw/89c76b2563dbc0e88384719a35cba0dfc04cd522/edgelist_sleeping_giant.csv')


