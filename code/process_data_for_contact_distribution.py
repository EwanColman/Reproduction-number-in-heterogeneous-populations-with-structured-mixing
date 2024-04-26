# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from time import time
from datetime import datetime
from numpy.random import choice
import math
import pickle as pk

# data from https://zenodo.org/communities/social_contact_data
country='Belgium'

contacts_data_path='../data/'+country+'/CoMix_BE_contact_common.csv'
participants_data_path='../data/'+country+'/CoMix_BE_participant_common.csv'
sday_data_path='../data/'+country+'/CoMix_BE_sday.csv'

contacts_df=pd.read_csv(contacts_data_path).dropna(subset=['cnt_age_est_max','cnt_age_est_min'],inplace=False)
participants_df=pd.read_csv(participants_data_path).dropna(subset=['part_age'],inplace=False)
sday_df=pd.read_csv(sday_data_path)

# isolate the wave here to avoide unneccessary computing

age_groups={'[0,1)':(0,1),
            '[1,6)':(1,6),
            '[6,12)':(6,12),
            '[12,18)':(12,18),
            '[18,30)':(18,30),
            '[30,40)':(30,40),
            '[40,50)':(40,50),
            '[50,60)':(50,60),
            '[60,70)':(60,70),
            '[70,120)':(70,120)}

# give each group a number
index_of={}
i=0
for age_group in age_groups:
    index_of[age_group]=i
    i=i+1

####### GET AGES ###########
# get age and degree of every node
ID_list=participants_df['part_id'].tolist()
age_group_list=participants_df['part_age'].tolist()
# store the age group of every participant in a dictionary
age_group_of={}
while ID_list:
    age=age_group_list.pop()
    ID=ID_list.pop()
    age_group_of[ID]=age

###### CONTACT DATA IN A DICTIONARY ######
# get the degree of every node categorized by age of contact

# get the degree of every node
degree={}
for ID in age_group_of:
    # start at zero then later start counting
    degree[ID]=0
# count number of times each ID appears in this list
part_list=contacts_df['part_id'].tolist()
while part_list:
    ID=part_list.pop()
    if ID in degree:
        degree[ID]=degree[ID]+1

print('Surveys:',len(sday_df))
print('Participants IDs:',len(list(set(sday_df['part_id']))))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
sday_df['days_since_t0']=[(datetime.strptime(str(d), '%Y.%m.%d')-time_zero).days for d in sday_df['sday_id'].tolist()]

plot_data={}

# choose a wave
for wave in [20,21]:
#wave=10
    period_sday_df=sday_df[sday_df['wave']==wave]
    # find the range of times
    print('Wave '+str(wave)+', from '+min(period_sday_df['sday_id'])+' to '+max(period_sday_df['sday_id'])+', '+str(max(period_sday_df['days_since_t0'])-min(period_sday_df['days_since_t0']))+' days')
    # make a list
    participant_list=period_sday_df['part_id'].tolist()
    # a bit of data cleaning...
    participant_list=[p for p in participant_list if p in age_group_of]
    
    
    # make n - the array of the number of individuals with 1, 2, etc. contacts in each age band
    
    
    dist={}
    for a in age_groups:
        dist[a]=[0 for i in range(10)]
    
    for participant in participant_list:
    
        a=age_group_of[participant]
        
        index=int(math.log2(degree[participant]+1))
        dist[a][index]=dist[a][index]+1
    
    while sum(dist[a][-1] for a in age_groups)==0:
        for a in age_groups:
            dist[a].pop()

    plot_data[wave]={'from':min(period_sday_df['sday_id']),
                     'to':max(period_sday_df['sday_id']),
                     'dist':dist
                     }

pk.dump(plot_data,open('../pickles/data_for_contact_distribution_plots.p','wb'))


