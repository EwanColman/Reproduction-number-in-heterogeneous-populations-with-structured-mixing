# -*- coding: utf-8 -*-
import pandas as pd
import pickle as pk
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from time import time
from datetime import datetime
from numpy.random import choice

group_names=['Under 1','1 to 5','6 to 11','12 to 17','18 to 29','30 to 39','40 to 49','50 to 59','60 to 69','70+']

# make the plot and set the gridspec
fig = plt.figure(figsize=(15,9))
gs = fig.add_gridspec(3,3)
plt.subplots_adjust(hspace=0.1,wspace=0.15)



######## PLOT THE CONTACT DISTRIBUTION ###########
plot_data=pk.load(open('../pickles/data_for_contact_distribution_plots.p','rb'))


wave=20
dist=plot_data[wave]['dist']
ax=fig.add_subplot(gs[0,0])  

plt.text(-0.15,1,'A',backgroundcolor='w',size=20, transform=ax.transAxes)


plt.title('$\\bf{Wave~'+str(wave)+'}$, '+plot_data[wave]['from']+' to '+plot_data[wave]['to'])
colors=plt.cm.Paired(np.linspace(0, 1, 12))
# start at the bottom
#bottom=[0 for i in range(10)]
w=0.08
n=0
for age_group in dist:
    print(dist[age_group])
    #plt.plot(range(10),dist[age_group],color=colors[n],linewidth=3,label=age_group)
    plt.bar([i+w*n-w*5 for i in range(len(dist[age_group]))],dist[age_group],w,color=colors[n],edgecolor='k',linewidth=0.1,label=group_names[n])
    #bottom=[bottom[i]+dist[age_group][i] for i in range(10)]
    n=n+1
plt.yscale('log')
tick_vals=[1,2,3,4,5,6,7,8]
plt.xticks([0]+tick_vals,['0']+[str((2**i)-1)+'-'+str((2**(i+1))-2) for i in tick_vals],rotation=60)
plt.yticks([1,10,100],[1,10,100])
plt.xlabel('Number of contacts')
plt.ylabel('Number of participants')
plt.ylim([0.8,200])
plt.xlim([-1,8.5])
    
#### Next plot
wave=21
dist=plot_data[wave]['dist']
ax=fig.add_subplot(gs[2,0])

plt.text(-0.15,1,'C',backgroundcolor='w',size=20, transform=ax.transAxes)

plt.title('$\\bf{Wave~'+str(wave)+'}$, '+plot_data[wave]['from']+' to '+plot_data[wave]['to'])
colors=plt.cm.Paired(np.linspace(0, 1, 12))
# start at the bottom
#bottom=[0 for i in range(10)]
w=0.08
n=0
for age_group in dist:
    
    #plt.plot(range(10),dist[age_group],color=colors[n],linewidth=3,label=age_group)
    plt.bar([i+w*n-w*5 for i in range(len(dist[age_group]))],dist[age_group],w,color=colors[n],edgecolor='k',linewidth=0.1,label=group_names[n])
    #bottom=[bottom[i]+dist[age_group][i] for i in range(10)]
    n=n+1
plt.yscale('log')
tick_vals=[1,2,3,4,5,6,7,8]
plt.xticks([0]+tick_vals,['0']+[str((2**i)-1)+'-'+str((2**(i+1))-2) for i in tick_vals],rotation=60)
plt.yticks([1,10,100],[1,10,100])
plt.xlabel('Number of contacts')
plt.ylabel('Number of participants')
plt.ylim([0.8,200])
plt.xlim([-1,8.5])


plt.legend(bbox_to_anchor=(1,1.7),ncol=3)
    

####### MATRICES #################
plot_data=pk.load(open('../pickles/data_for_matrix_plots.p','rb'))


# wave 20 first
C=plot_data[20]
ax=fig.add_subplot(gs[0,1])   

plt.text(-0.35,1,'B',backgroundcolor='w',size=20, transform=ax.transAxes)
 
plt.xticks(range(10),group_names,rotation=60)
plt.yticks(range(10),group_names)
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
im=plt.imshow(np.transpose(C),vmin=0,vmax=np.log(1+4.26),cmap='Greys',origin='lower')
plt.xlabel('Age of participant')
plt.ylabel('Age of contact')



# wave 21 second
C=plot_data[21]
ax=fig.add_subplot(gs[2,1])

plt.text(-0.35,1,'D',backgroundcolor='w',size=20, transform=ax.transAxes)
 
plt.xticks(range(10),group_names,rotation=60)
plt.yticks(range(10),group_names)
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
im=plt.imshow(np.transpose(C),vmin=0,vmax=np.log(1+4.26),cmap='Greys',origin='lower')
plt.xlabel('Age of participant')
plt.ylabel('Age of contact')

cbar_ax = fig.add_axes([0.45, 0.49, 0.13, 0.01])
cbar=fig.colorbar(im, cax=cbar_ax, shrink=1,ticks=[np.log(1+i) for i in range(5)],orientation='horizontal')
cbar.set_label('Mean number of contacts',rotation=0)
cbar.ax.set_xticklabels([0,1,2,3,4])  # vertically oriented colorbar
#####################################

w=1/20
nbins=int(8/w)
x_vals=[i*w for i in range(nbins)]

lstyle={'Well-mixed & homogeneous':'-', 
       'Age-dependent & homogeneous':'--', 
       'Well-mixed & heterogeneous':'-.', 
       'Age-dependent & heterogeneous':':'}

colour={'Well-mixed & homogeneous':'orange', 
        'Age-dependent & homogeneous':'tomato', 
        'Well-mixed & heterogeneous':'chartreuse', 
        'Age-dependent & heterogeneous':'gold'}

### Botstrap R0 plots
plot_data=pk.load(open('../pickles/data_for_R0_bootsrap_plots','rb'))

ax=fig.add_subplot(gs[0:3,2])

plt.text(-0.1,1,'E',backgroundcolor='w',size=20, transform=ax.transAxes)
plt.xlabel('$R_{0}$ estimate')
data=[]

for calc in lstyle:
    mean={}
    for wave in [20,21]:
        data.append(sorted(plot_data[wave][calc]))

        mean[wave]=np.mean(plot_data[wave][calc])

    print(calc,'has a percent increase of',round(100*(mean[21]-mean[20])/mean[20],3))


z=0.15
parts=ax.violinplot(data,widths=2*z,positions=[4+z,4-z,3+z,3-z,2+z,2-z,1+z,1-z],showmeans=False, showmedians=False,showextrema=False,vert=False)

for pc in [parts['bodies'][i] for i in [0,2,4,6]]:
    
    pc.set_facecolor('blue')
    pc.set_edgecolor('blue')
    pc.set_alpha(0.7)
    
for pc in [parts['bodies'][i] for i in [1,3,5,7]]:
    
    pc.set_facecolor('purple')
    pc.set_edgecolor('purple')
    pc.set_alpha(0.7)

y=0.6
for x in [1+y,2+y,3+y]:
    plt.axhline(x,linestyle=':',linewidth=0.5,color='k')
plt.ylim([y,4+y])
plt.yticks([4+z,4-z,3+z,3-z,2+z,2-z,1+z,1-z],['20','21','20','21','20','21','20','21'])


i=0
for model in lstyle.keys():
    i=i+1
    plt.text(0.5,5-i+3*z,model)
    plt.text(-0.9,5.25-i,'Wave')#,rotation=90)
#plt.xticks([i-0.3 for i in range(1,5)],lstyle.keys(),rotation=60)
#plt.legend(title='$\\bf{Model~assumptions:}$                         ',bbox_to_anchor=(0.7, -0.3))
plt.savefig('../figures/figure3.png',format='png',dpi=300,bbox_inches='tight')
