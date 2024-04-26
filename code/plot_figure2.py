# -*- coding: utf-8 -*-
import pandas as pd
import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
from time import time
from datetime import datetime

moments={'Primary and secondary schools reopening':'2021-04-20',
         'End of restrictions on higher education & restaurants':'2021-08-25',
         'Restrictions introduced on large gaherings and restaurants':'2021-10-28'}

time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(4,10)]\
    +['01/'+str(i)+'/2020' for i in range(10,13)]\
    +['01/0'+str(i)+'/2021' for i in range(1,10)]\
    +['01/'+str(i)+'/2021' for i in range(10,13)]\
    +['01/'+str(i)+'/2022' for i in range(1,4)]
    
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]


# make the plot and set the gridspec
fig = plt.figure(figsize=(15,5))
gs = fig.add_gridspec(2,2,width_ratios=[3,1])
plt.subplots_adjust(hspace=0,wspace=0.2)


plot_data=pk.load(open('../pickles/data_for_all_wave_plots_unfiltered.p','rb'))

ax=fig.add_subplot(gs[0,0])
for x in dates_numerical:
    plt.axvline(x,linestyle=':',linewidth=0.5,color='k')

plt.text(50,75,'All participants included',backgroundcolor='w')
plt.xticks([])
plt.yticks([0,20,40,60])
plt.xlim([35,740])
plt.ylim([0,90])
plt.ylabel('R0 estimate')

plt.text(-8,80,'A',backgroundcolor='w',size=20)


for wave in plot_data:

    date=plot_data[wave]['Date']

    t=(datetime.strptime(str(date), '%Y.%m.%d')-time_zero).days
    
    R0_bootstraps=sorted(plot_data[wave]['R0_bootstraps'])
    
    R0_lower=R0_bootstraps[int(len(R0_bootstraps)*0.975)]
    R0_median=R0_bootstraps[int(len(R0_bootstraps)*0.5)]
    R0_upper=R0_bootstraps[int(len(R0_bootstraps)*0.025)]
    
    plt.scatter([t],[R0_median],marker='s',c='grey')
    plt.plot([t,t],[R0_lower,R0_upper],c='grey')
    

for moment in moments:
    # convert to numerical
    t=(datetime.strptime(str(moments[moment]), '%Y-%m-%d')-time_zero).days
    plt.axvline(t,linestyle=':',linewidth=2,color='r')

plot_data=pk.load(open('../pickles/data_for_all_wave_plots_filtered.p','rb'))
ax=fig.add_subplot(gs[1,0])
for x in dates_numerical:
    plt.axvline(x,linestyle=':',linewidth=0.5,color='k')

plt.text(50,7,'Outliers removed',backgroundcolor='w')
plt.xticks(dates_numerical,dates_words,rotation=60)
plt.yticks([0,2,4,6])
plt.xlim([35,740])
plt.ylim([0,8.5])
plt.ylabel('R0 estimate')
for wave in plot_data:

    date=plot_data[wave]['Date']
    t=(datetime.strptime(str(date), '%Y.%m.%d')-time_zero).days
    
    R0_bootstraps=sorted(plot_data[wave]['R0_bootstraps'])
    
    R0_lower=R0_bootstraps[int(len(R0_bootstraps)*0.975)]
    R0_median=R0_bootstraps[int(len(R0_bootstraps)*0.5)]
    R0_upper=R0_bootstraps[int(len(R0_bootstraps)*0.025)]
    
    plt.scatter([t],[R0_median],marker='s',c='grey')
    plt.plot([t,t],[R0_lower,R0_upper],c='grey')


    
moment='Primary and secondary schools reopening'
# convert to numerical
t=(datetime.strptime(str(moments[moment]), '%Y-%m-%d')-time_zero).days
plt.axvline(t,linestyle=':',linewidth=2,color='r')
plt.plot([t-30,t],[-3,0],linestyle=':',linewidth=2,color='r',clip_on=False)
plt.text(t-30,-4,moment,horizontalalignment='right')

moment='End of restrictions on higher education & restaurants'
t=(datetime.strptime(str(moments[moment]), '%Y-%m-%d')-time_zero).days
plt.axvline(t,linestyle=':',linewidth=2,color='r')
plt.plot([t-40,t],[-4,0],linestyle=':',linewidth=2,color='r',clip_on=False)
plt.text(t-40,-5,moment,horizontalalignment='right')

moment='Restrictions introduced on large gaherings and restaurants'
t=(datetime.strptime(str(moments[moment]), '%Y-%m-%d')-time_zero).days
plt.axvline(t,linestyle=':',linewidth=2,color='r')
plt.plot([t-50,t],[-5,0],linestyle=':',linewidth=2,color='r',clip_on=False)
plt.text(t-50,-6,moment,horizontalalignment='right')

start=35
end=740
# years on axes
lh= -2#-(m+10)*0.25
th= -3#-(m+10)*0.35
year1=(datetime.strptime(str('01/01/2021'), '%d/%m/%Y')-time_zero).days
year2=(datetime.strptime(str('01/01/2022'), '%d/%m/%Y')-time_zero).days
year3=(datetime.strptime(str('01/01/2023'), '%d/%m/%Y')-time_zero).days
plt.plot([start,year1-15],[lh,lh],c='k',linewidth=1,clip_on=False)
plt.text((start+year1)/2-20,th,'2020')
plt.plot([year1+15,year2-15],[lh,lh],c='k',linewidth=1,clip_on=False)
plt.text((year1+year2)/2-20,th,'2021')
plt.plot([year2+15,end],[lh,lh],c='k',linewidth=1,clip_on=False)
plt.text((year2+end)/2-10,th,'2022')

####### NEXT ONE ########

plot_data=pk.load(open('../pickles/data_for_R0_vs_cutoff_plots.p','rb'))
ax=fig.add_subplot(gs[0,1])
ax.set_xlim([2,6.8])
ax.set_ylim([0,27])
ax.set_yticks([0,5,10,15,20])
ax.set_xticks([])
ax.set_ylabel('Estimated R')

#impact_of_removals=[i for i in plot_data[wave]]

lower,upper,median=[],[],[]
for i in range(25):
    # make list
    R0_list=sorted([plot_data[wave]['R0'][i] for wave in plot_data])
    lower.append(R0_list[0])
    upper.append(R0_list[-1])
    median.append(R0_list[int(len(R0_list)/2)])

plt.plot([0.2*i for i in range(10,35)],median,linewidth=3,c='k')
plt.fill_between([0.2*i for i in range(10,35)],lower,upper,color='k',linewidth=0,alpha=0.1)   

plt.text(1.2,24,'B',backgroundcolor='w',size=20)
#for wave in plot_data:
    
    #R0_with_no_removals=plot_data[wave]['R0'][-1]
    #impact_of_removals=[100*(R0_with_no_removals-i)/R0_with_no_removals for i in plot_data[wave]]
    #impact_of_nth_removal=[(plot_data[wave][i-1]-plot_data[wave][i])/plot_data[wave][i-1] for i in range(1,20)]
    #ax.plot([0.2*i for i in range(10,35)],plot_data[wave]['R0'][:25],c='k',alpha=0.1)
    

#ax2 = ax.twinx()
ax2=fig.add_subplot(gs[1,1])
#ax2.spines['right'].set_color('g')
#ax2.tick_params('y', colors='g')
ax2.set_xlim([2,6.8])
ax2.set_yticks([0,20,40,60])
ax2.set_ylim([0,90])
ax2.set_ylabel('Outliers omitted')
ax2.set_xlabel('Critical z-score for outlier removal')
#for wave in plot_data:
#    ax2.plot([0.2*i for i in range(10,35)],plot_data[wave]['number of outliers'][:25],c='r',alpha=0.1)

lower,upper,median=[],[],[]
for i in range(25):
    # make list
    R0_list=sorted([plot_data[wave]['number of outliers'][i] for wave in plot_data])
    lower.append(R0_list[0])
    upper.append(R0_list[-1])
    median.append(R0_list[int(len(R0_list)/2)])

plt.plot([0.2*i for i in range(10,35)],median,linewidth=3,c='g')
plt.fill_between([0.2*i for i in range(10,35)],lower,upper,color='g',linewidth=0,alpha=0.1)   


plt.savefig('../figures/figure2.png',format='png',dpi=300,bbox_inches='tight')

