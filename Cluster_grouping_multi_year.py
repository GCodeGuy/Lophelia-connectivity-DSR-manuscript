# -*- coding: utf-8 -*-
"""
Created on Sun Feb 26 12:58:23 2023

@author: Graem
"""

import numpy as np
import pandas as pd
import leidenalg as la
import igraph as ig

#%%
'''Set Parameters'''
seasons_arr = ["Winter","Spring","Summer","Fall"]
year_arr=['B2005','B2006','B2007','B2008','B2009','B2010','B2011','B2012','B2013','B2014','B2015', 'B2016','B2017','B2018'] 
dep_arr=['Seabed']
dep = 'Seabed'
month_num_arr = [2,4] #,5,7,9,11]
season_nums = [1,2,3,4]
cell_size= '0.1'
particle_space = '0.003' #, '0.005'] #unit:degree
comp_arr=['30'] #,'20']
dt='10'  #time step minutes
day1='1' #output results days
day2='60' #final duration
repeat = '0' #release particles every x days
durations = np.array([61])
idx_2=durations[0]
kh='VARkh'
trial = 'Bdom_Lophelia_final_run'
trial2 = 'Bdom_Loph_linearVELincrease'
shtr = "LVELinc"
vertvel='VARvel'
#%%

year = 'B2005'
season = 'Winter'
particle_space = '0.003'
comp = '30'
dep = 'Seabed'
season_num = '1'
#%%
coords = np.loadtxt("Path to nodal coordinates.txt")
coords[:,2] = range(len(coords))
coords = pd.DataFrame(coords)
coords = coords.sort_values(by=[1])
coords = np.array(coords)
#%%
'''Partition nodes based on connectivity matrix'''

year_list = []
season_list = []
comp_list = []
q_list = []
partition_count = []

#for idc in range(len(comp_arr)):
    #comp = comp_arr[idc]
for idy in range(len(year_arr)):
    year = year_arr[idy]
    for ids in range(len(seasons_arr)):
        season = seasons_arr[ids]
        season_num = season_nums[ids]
        #for idy in range(len(year_arr)):
        #year = year_arr[idy]  
        count_crossing = np.loadtxt("Path to connectivity matrix.txt")
        count_crossing = count_crossing.round(decimals = 0)  
         
        cc_rearrange = np.zeros([len(count_crossing),len(count_crossing)])
        
        for i in range(len(count_crossing)):  ## Rearranges node in matrix by longitude from West to East
            for j in range(len(count_crossing)):
                cc_rearrange[i,j] = count_crossing[int(coords[i,2]), int(coords[j,2])]
        
        G = ig.Graph.Weighted_Adjacency(cc_rearrange, mode = 'directed')
        G_dict = G.to_dict_dict(edge_attrs='weight')
        
        
        # Set initial partition
        partition = la.ModularityVertexPartition(G , weights = G.es['weight'])
        refined_partition = la.ModularityVertexPartition(G, weights = G.es['weight'])
        partition_agg = refined_partition.aggregate_partition()
        optimiser = la.Optimiser()
        optimiser.set_rng_seed(0)
        
        while optimiser.move_nodes(partition_agg):
        
          # Get individual membership for partition
          partition.from_coarse_partition(partition_agg, refined_partition.membership)
        
          # Refine partition
          refined_partition = la.ModularityVertexPartition(G, weights = G.es['weight'])
          optimiser.merge_nodes_constrained(refined_partition, partition)
        
          # Define aggregate partition on refined partition
          partition_agg = refined_partition.aggregate_partition()
        
        
          # Use membership of actual partition
          aggregate_membership = [None] * len(refined_partition)
          for i in range(G.vcount()):
            aggregate_membership[refined_partition.membership[i]] = partition.membership[i]
          partition_agg.set_membership(aggregate_membership)
        
        partition_list = list(partition)
        refined_partition_list = list(refined_partition)
        partition_agg_list = list(partition_agg)
        
        partition_df = pd.DataFrame(partition_list).T
        refined_partition_df = pd.DataFrame(refined_partition_list).T
        partition_agg_df = pd.DataFrame(partition_agg_list).T

        season_list.append(season)
        comp_list.append(comp)
        q_list.append(partition.q)
        partition_count.append(len(partition_list))
        year_list.append(year)
        
    
        '''Read in node table and attach cluster groupings and refined cluster groupings'''
    
        
        node_table = pd.read_excel("Path to node table.txt")
        node_table['cluster_group']=np.zeros(len(node_table))
        node_table['refined_partition_group']=np.zeros(len(node_table))
        
        for i in range(len(partition_df.columns)):
            rows = np.array(partition_df[i])
            rows = rows[np.isfinite(rows)]
            
            for l in range(len(rows)):       
                node_table['cluster_group'][partition_df[i][l]] = partition.membership[int(partition_df[i][l])]
                
        for i in range(len(refined_partition_df.columns)):
            rows = np.array(refined_partition_df[i])
            rows = rows[np.isfinite(rows)]
            
            for l in range(len(rows)):
                node_table['refined_partition_group'][refined_partition_df[i][l]] = refined_partition.membership[int(refined_partition_df[i][l])]
                  
            
                 
        node_table.to_excel("Path to updated node table.csv")
                
 #%%
## Create data table containing only clustering information
data_table = np.zeros([len(q_list),5])
data_table =  pd.DataFrame(data = data_table,
                          columns = ['year','season','comp','clusters','q'])

data_table['year'] = year_list
data_table['season'] = season_list
data_table['comp'] = comp_list
data_table['clusters'] = partition_count
data_table['q'] = q_list
dt_header =  ['year','season','comp','clusters','q']

data_table.to_excel('Path to data table.csv')
   

   
    
        


#%%

'''Read in node tables and calculate Google Pagerank for each INDIVIDUAL cluster'''

node_table = pd.read_csv("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\"+trial+dep+'Long_node_table_nts_epd'+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+comp+'comp'+repeat+"dt_B"+kh+vertvel+".csv", sep=',')
node_table['cluster_group']=np.zeros(len(node_table))
node_table['pagerank']=np.zeros(len(node_table))

for i in range(len(partition_df.columns)):
    rows = np.array(partition_df[i])
    rows = rows[np.isfinite(rows)]
    columns = np.array(partition_df[i])
    columns = columns[np.isfinite(columns)]
    
    grouped_cc = np.zeros([len(rows), len(columns)])
    for r in range(len(rows)):
        for c in range(len(columns)):
            grouped_cc[r,c] = cc_rearrange[int(rows[r]), int(columns[c])]
    
    
    G_cluster = ig.Graph.Weighted_Adjacency(grouped_cc, mode = 'directed')
    
    print(G_cluster.is_connected())
    #G_centrality = G_cluster.eigenvector_centrality(directed = True, weights='weight') 
    G_cluster_pagerank = G_cluster.pagerank(directed = True, weights = 'weight')
    for l in range(len(rows)):
        node_table['pagerank'][partition_df[i][l]] = G_cluster_pagerank[l]
        node_table['cluster_group'][partition_df[i][l]] = partition.membership[int(partition_df[i][l])]
            
        
        node_table.to_csv('C:\\Users\\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\ComputeCanada Output\\'+trial+dep+"PageRL_node_table_ntsepd"+year_arr[0]+'-'+year_arr[-1]+season+str(idx_2)+"days_"+dt+"minutes_"+str(cell_size)+"cell_"+particle_space+"degree_"+str(comp)+"comp"+repeat+"dt"+"_B"+str(kh)+vertvel+".csv", sep=',')
        
        partition_df.to_csv("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\" +season+"_all_partitions"+comp+"comp"+".csv", index=False, header=False)
        refined_partition_df.to_csv("C:\\Users\Graem\\OneDrive - Dalhousie University\\Documents\\Dal 2020-2021\\Graduate School\\Graduate Project\\Chapter 1 Analysis\\" +season+"_refined_partitions"+comp+"comp"+".csv", index=False, header=False)





