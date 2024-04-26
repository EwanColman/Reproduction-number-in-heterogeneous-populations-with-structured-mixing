# -*- coding: utf-8 -*-
import pandas as pd
import pickle as pk
import numpy as np
from numpy import linalg as LA
from datetime import datetime
from numpy.random import choice
np.random.seed(0)
# data from https://zenodo.org/communities/social_contact_data
#country='Slovenia'
beta=0.2
z_score_threshold=3.5
nboots=10000
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

    min_age=int(row['cnt_age_est_min'])#int(min_age_list.pop(0))
    max_age=int(row['cnt_age_est_max'])#int(max_age_list.pop(0))
    # get fraction of population in the age range that are in each exact age
    weights=[age_dist[a]/sum(age_dist[min_age:max_age]) for a in range(min_age,max_age)]
    # # choose one randomly using these weights
    random_age=choice(range(min_age,max_age),p=weights)
    contact_age.append(random_age)
# add this as a new column to the dataframe
pk.dump(contact_age,open('../pickles/contact_age_temp.p','wb'))
contacts_df['cnt_age']=contact_age

#contacts_df['cnt_age']=pk.load(open('../pickles/contact_age_temp.p','rb'))

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
    #participant_list=[p for p in participant_list if p in age_group_of and sum(degree[p][age_group] for age_group in age_groups)<50]

    participant_list=[p for p in participant_list if p in age_group_of]# and sum(degree[p][age_group] for age_group in age_groups)<cut_off]    
    degree_list=[(p,sum(degree[p][age_group] for age_group in age_groups)) for p in participant_list]
    #degree_list=sorted(degree_list, key = lambda item: item[1])
    #participant_list=[d[0] for d in degree_list]
    #participant_list=participant_list[0:len(participant_list)-removals]
   
    # find the degree corresponding to the z-score
    mu=np.mean([np.log(1+d[1]) for d in degree_list])
    sd=np.std([np.log(1+d[1]) for d in degree_list])
   
    # cutoff is the degree that is z sds above mean 
    #use threold of 4
    cut_off=np.exp(mu+z_score_threshold*sd)-1
    print('cut off:',cut_off)
    print('Removed:',len([p for p in participant_list if sum(degree[p][age_group] for age_group in age_groups)>=cut_off]))
    participant_list=[p for p in participant_list if sum(degree[p][age_group] for age_group in age_groups)<cut_off]
    
    largest_degree=max([sum(degree[p][age_group] for age_group in age_groups) for p in participant_list])
    
    # make the samples
    bootstrap_samples=[]
    boot=0
    while boot<nboots:
        # bootstrapping: select sample from the participant list
        sample=choice(participant_list,len(participant_list))
        # check that all age groups have at least one participant in the sample
        #if len(set([age_group_of[participant] for participant in sample]))==len(age_groups):
        boot=boot+1
        bootstrap_samples.append(sample)
    
    plot_data[wave]={'Well-mixed & homogeneous':[], 
               'Age-dependent & homogeneous':[], 
               'Well-mixed & heterogeneous':[], 
               'Age-dependent & heterogeneous':[]
               }

    count=0    
    for bootstrap_sample in bootstrap_samples:
        print()
        print(count)
        count=count+1
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
        
        for participant in bootstrap_sample:    
            a=age_group_of[participant]
            for b in age_groups:
                #
                number_of_contacts[(a,b)]=number_of_contacts[(a,b)]+degree[participant][b]

            number_of_participants[a]=number_of_participants[a]+1
    
        # make the edges matrix with population adjustment
        C=[]
        number_of_edges={}
        for a in age_groups:
            c=[]
            #c=[]
            for b in age_groups:
                # using the correction formula        
                population_contacts=(1/2)*((number_of_contacts[(a,b)]*population[a]/max(1,number_of_participants[a]))+(number_of_contacts[(b,a)]*population[b]/max(1,number_of_participants[b])))
                c.append(population_contacts)
            C.append(c)
        C=np.array(C)
    
        # make n - the array of the number of individuals with 1, 2, etc. contacts in each age band
        K=largest_degree+1
        n=[[0 for a in age_groups] for k in range(K)]
    
    
        for participant in bootstrap_sample:
            # get number of contacts
            k=sum(degree[participant][age_group] for age_group in age_groups)
            # get their age
            a=age_group_of[participant]
            # add it to the array (+1 but do the population weighting here too)
            n[k][index_of[a]]=n[k][index_of[a]]+(population[a]/max(1,number_of_participants[a]))
        
        n=np.array(n)
    
        ages=len(C)
        N=sum(n[:])
        c=sum(C[:])
        
        # naive method
        mean_degree=np.mean([sum(degree[participant][age_group] for age_group in age_groups) for participant in bootstrap_sample])
        print('Well-mixed & homogeneous:',mean_degree*beta)
        plot_data[wave]['Well-mixed & homogeneous'].append(mean_degree*beta)
        #### age structure, no heterogeneity
        Z=[]
        for a in range(ages):
            z=[]
            for b in range(ages):
                z.append(beta*C[a,b]/max(1,N[b]))
            Z.append(z)
        Z=np.array(Z)
        #print(Z)
        
        #print('Old approach')
        eigenvalues,eigenvectors=LA.eig(Z)
        #print('R0=')
        R0=max([np.real(z) for z in eigenvalues])
        #t1=t1+time()-start
        print('Age-dependent & homogeneous',R0)
        plot_data[wave]['Age-dependent & homogeneous'].append(R0)
        
        #### Heterogeneity, no age structure    
        R0=beta*(1/sum(c))*sum([sum([(k**2)*n[k,b] for k in range(0,K)]) for b in range(ages)])
        print('Well-mixed & heterogeneous',R0)
        plot_data[wave]['Well-mixed & heterogeneous'].append(R0)
    
        ### Combined method
        M=[]
        for a in range(ages):
            m=[]
            for b in range(ages):
                m.append(C[a,b]/max(1,c[a]*c[b])) 
            M.append(m)
        
        M=np.array(M)
        # use the formula
        G=[]
        for a in range(ages):
            g=[]
            for b in range(ages):
                #print(a,b,beta*M[b,a]*sum([(k**2)*n[k-1,b] for k in range(1,K+1)]))
                g.append(beta*M[b,a]*sum([(k**2)*n[k,b] for k in range(0,K)]))
            G.append(g)
            
        G=np.array(G)
        
        eigenvalues,eigenvectors=LA.eig(G)
        R0=max([np.real(z) for z in eigenvalues])
        print('Age-dependent & heterogeneous',R0)
        plot_data[wave]['Age-dependent & heterogeneous'].append(R0)

pk.dump(plot_data,open('../pickles/data_for_R0_bootsrap_plots','wb'))
    

