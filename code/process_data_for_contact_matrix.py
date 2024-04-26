# -*- coding: utf-8 -*-
import pandas as pd
import pickle as pk
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from time import time
from datetime import datetime
from numpy.random import choice
np.random.seed(0)
# data from https://zenodo.org/communities/social_contact_data

country='Belgium'

contacts_data_path='../data/'+country+'/CoMix_BE_contact_common.csv'
participants_data_path='../data/'+country+'/CoMix_BE_participant_common.csv'
sday_data_path='../data/'+country+'/CoMix_BE_sday.csv'

contacts_df=pd.read_csv(contacts_data_path).dropna(subset=['cnt_age_est_max','cnt_age_est_min'],inplace=False)
participants_df=pd.read_csv(participants_data_path).dropna(subset=['part_age'],inplace=False)
sday_df=pd.read_csv(sday_data_path)


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


############
# population data
pop_df=pd.read_csv('../data/Populations/Euro_pops.csv')
# slice out the country and year
pop_df=pop_df[(pop_df['Region, subregion, country or area *']==country)&(pop_df['Year']==2020)]
# make a list of weights for sampling (over 100s spread out uniformly)
age_dist=[pop_df[str(i)].iloc[0] for i in range(100)]+[pop_df['100+'].iloc[0]/20  for i in range(21)]

# store population counts in a dictionary
population={}
for age in age_groups:
    # sum over all the ages in the group
    population[age]=sum(age_dist[age_groups[age][0]:age_groups[age][1]])

####### IMPUTE AGES #########
# ages are given as ranges, so use population data to guess the exact age
contact_age=[]
min_list=contacts_df['cnt_age_est_min'].tolist()
max_list=contacts_df['cnt_age_est_max'].tolist()

for i,row in contacts_df.iterrows():
#while min_list:
#    min_age=int(min_list.pop())
#    max_age=int(max_list.pop())
    #print(len(min_list))
#    print(i)
    min_age=int(row['cnt_age_est_min'])#int(min_age_list.pop(0))
    max_age=int(row['cnt_age_est_max'])#int(max_age_list.pop(0))
    # get fraction of population in the age range that are in each exact age
    weights=[age_dist[a]/sum(age_dist[min_age:max_age]) for a in range(min_age,max_age)]
    # # choose one randomly using these weights
    random_age=choice(range(min_age,max_age),p=weights)
    contact_age.append(random_age)
# add this as a new column to the dataframe
contacts_df['cnt_age']=contact_age

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
degree={}
for ID in age_group_of:
    # start at zero then later start counting
    degree[ID]={}
    for age_group in age_groups:
        # number of contacts of ID in age group 
        degree[ID][age_group]=0

# count number of times each ID appears in this list
part_list=contacts_df['part_id'].tolist()
age_list=contacts_df['cnt_age'].tolist()
while part_list:
    ID=part_list.pop()
    age=age_list.pop()
    # get the age group

    for a in age_groups:
        if age>=age_groups[a][0] and age<age_groups[a][1]:
            age_group=a
    
    if ID in degree:
        degree[ID][age_group]=degree[ID][age_group]+1

print('Surveys:',len(sday_df))
print('Participants IDs:',len(list(set(sday_df['part_id']))))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
sday_df['days_since_t0']=[(datetime.strptime(str(d), '%Y.%m.%d')-time_zero).days for d in sday_df['sday_id'].tolist()]

# inuitialise the output
plot_data={}

# choose a wave
#for wave in range(1,44):
for wave in [20,21]:
    period_sday_df=sday_df[sday_df['wave']==wave]
    # find the range of times
    print('Wave '+str(wave)+', from '+min(period_sday_df['sday_id'])+' to '+max(period_sday_df['sday_id'])+', '+str(max(period_sday_df['days_since_t0'])-min(period_sday_df['days_since_t0']))+' days')
    # make a list
    participant_list=period_sday_df['part_id'].tolist()
    # a bit of data cleaning...
    participant_list=[p for p in participant_list if p in age_group_of]
    
    ### Weigh by day of week #############
    weekend_participants=sday_df[(~sday_df['dayofweek'].isin([1,2,3,4,5]))]['part_id'].tolist()
    weekday_participants=sday_df[(sday_df['dayofweek'].isin([1,2,3,4,5]))]['part_id'].tolist()

    # reduce to only those in the participant list
    weekend_participants=[p for p in weekend_participants if p in participant_list]
    weekday_participants=[p for p in weekday_participants if p in participant_list]

    print('Weekend:',len(weekend_participants))
    print('Weekday:',len(weekday_participants))
    day_weight={}
    for participant in weekend_participants:
        day_weight[participant]=(2/7)/(len(weekend_participants)/(len(weekend_participants)+len(weekday_participants)))
    for participant in weekday_participants:
        day_weight[participant]=(5/7)/(len(weekday_participants)/(len(weekend_participants)+len(weekday_participants)))

    ###############
    
    
    
    # calculate the mixing matrix from the degree dictionary
    # count the number of participants and contacts by age group
    
    ####### CONTACT MATRIX #########
    # initialise dictionary
    number_of_contacts={}
    # and the dictionary for number in each age group
    number_of_participants={}
    
    for a in age_groups:
        number_of_participants[a]=0
        for b in age_groups:
            number_of_contacts[(a,b)]=0
    
    for participant in participant_list:    
        a=age_group_of[participant]
        for b in age_groups:
            #print((a,b))
            number_of_contacts[(a,b)]=number_of_contacts[(a,b)]+degree[participant][b]*day_weight[participant]
            #mean_number_of_contacts[(a,b)]=mean_number_of_contacts[(a,b)]+(degree[participant][b]/nboots)
        
        number_of_participants[a]=number_of_participants[a]+1
    
    
    total_contacts=sum(number_of_contacts.values())
    # make the edges matrix with
    C=[]
    number_of_edges={}
    for a in age_groups:
        c=[]
        #c=[]
        for b in age_groups:
            # using the correction formula        
            #c.append(100*number_of_contacts[(a,b)]/total_contacts)
           
            mean_contacts=(1/(2*population[a]))*((number_of_contacts[(a,b)]*population[a]/number_of_participants[a])+(number_of_contacts[(b,a)]*population[b]/number_of_participants[b]))
            
            #print('People in',a,'have a mean of',mean_contacts,'contacts in',b)
            c.append(np.log(1+mean_contacts))
           
        C.append(c)
    C=np.array(C)
    
    plot_data[wave]=C
 
pk.dump(plot_data,open('../pickles/data_for_matrix_plots.p','wb'))


