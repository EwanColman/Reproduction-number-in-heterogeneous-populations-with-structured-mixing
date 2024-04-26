import pickle as pk
import matplotlib.pyplot as plt
import numpy as np

x_values=[i/10 for i in range(10)]
CV_values=[i/10 for i in range(10)]
output=pk.load(open('../pickles/R0_analytical_and_simulated.p','rb'))


R={'Simulation':[]}
calcs=['Well-mixed & homogeneous', 
       'Group-dependent & homogeneous', 
       'Well-mixed & heterogeneous', 
       'Group-dependent & heterogeneous']


### for the contact distribution plots, make three networks
from adjacency_list_generator import get_adjacency_dictionary
x=0
degree_dist={}
for CV in [0,0.4,0.8]:
    neighbours,nodes_in=get_adjacency_dictionary(10000,100000,x,CV)
    degree={}
    for node in neighbours:
        degree[node]=len(neighbours[node])

    
    degree_dist[CV]=[0 for i in range(40)]
    for d in degree.values():
        degree_dist[CV][d]=degree_dist[CV][d]+1/sum(degree.values())

### for the mini heat maps, make three contact matrices
CV=0
CM={}
for x in [0,0.4,0.8]:
    CM[x]=np.array([[1-x,1-x],[1-x,2*x+1]])
###############


for x in x_values:
    r=[]
    for CV in CV_values:
        #print(output[(x,CV)].keys())
        
        # plt.figure()
        # plt.title('x='+str(x)+', y='+str(CV))
        # plt.xlabel('Generation')
        # plt.ylabel('Reproductive ratio')
        # plt.ylim([0,10])

        # for calculating the mean after n generations
        converged_R0_list=[]
        
        for R0_series in output[(x,CV)]['R0 by generation']:
        
            if R0_series[-1]>0:
                converged_R0_list.append(R0_series[-1])
                #plt.plot(R0_series,'k',alpha=0.1)
        
           
        #print(x,CV,np.mean(converged_R0_list))
        r.append(np.mean(converged_R0_list))
    R['Simulation'].append(r)
R['Simulation']=np.array(R['Simulation'])


for calc in calcs:    
    R[calc]=[]
    for x in x_values:
        r=[]
        for CV in CV_values:
            #print(x,CV,np.mean(converged_R0_list))
            r.append(output[(x,CV)][calc])
        R[calc].append(r)    
    
    R[calc]=np.array(R[calc])
    


# make the plot and set the gridspec
fig = plt.figure(figsize=(15,7))
gs = fig.add_gridspec(12, 24)
plt.subplots_adjust(hspace=1,wspace=1)

#create an axis
# min-plots
ax=fig.add_subplot(gs[0:2,0:2])

### HISTOGRAM the degree distribution ######
plt.bar(range(40),degree_dist[0.8],width=1)
plt.xlabel('Contacts')
plt.ylabel('Freq.')
plt.xticks([0,20])
plt.ylim([0,0.02])
plt.yticks([])
plt.text(0,0.015,'h=0.8')

ax=fig.add_subplot(gs[3:5,0:2])
plt.bar(range(40),degree_dist[0.4],width=1)
plt.xlabel('Contacts')
plt.ylabel('Freq.')
plt.xticks([0,20])
plt.ylim([0,0.02])
plt.yticks([])
plt.text(0,0.015,'h=0.4')

ax=fig.add_subplot(gs[6:8,0:2])
plt.bar(range(40),degree_dist[0],width=1)
plt.xlabel('Contacts')
plt.ylabel('Freq.')
plt.xticks([0,20])
plt.ylim([0,0.02])
plt.yticks([])
plt.text(0,0.015,'h=0')

# big heatmap
ax=fig.add_subplot(gs[0:9,3:12])
plt.text(-0.11,1.01,'A',size=20, transform=ax.transAxes)

#### plotting ####
plt.xticks(range(len(x_values)),x_values,rotation=60)
plt.yticks(range(len(CV_values)),CV_values)
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
im=plt.imshow(np.transpose(R[calc]),cmap='Greys',vmax=14,vmin=4,origin='lower') # displays in color
plt.xlabel('Group-dependent mixing, m')
plt.ylabel('Contact heterogeneity, h')
cbar=fig.colorbar(im, ax=ax, shrink=0.8,ticks=[4,6,8,10,12,14])
cbar.set_label('$R_{0}$',rotation=0)
cbar.ax.set_yticklabels([4,6,8,10,12,14])  # vertically oriented colorbar
######################################

# mini contact plots
ax=fig.add_subplot(gs[10:12,3:5])
plt.xticks(range(len(nodes_in)),nodes_in.keys())
plt.yticks(range(len(nodes_in)),nodes_in.keys())
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.imshow(np.transpose(CM[0]),vmax=3,vmin=0,cmap='Greys',origin='lower') # displays in color
plt.text(0.5,-0.3,'m=0')


ax=fig.add_subplot(gs[10:12,6:8])
plt.xticks(range(len(nodes_in)),nodes_in.keys())
plt.yticks(range(len(nodes_in)),nodes_in.keys())
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.imshow(np.transpose(CM[0.4]),vmax=3,vmin=0,cmap='Greys',origin='lower') # displays in color
plt.text(0.2,-0.3,'m=0.4')



ax=fig.add_subplot(gs[10:12,9:11])
plt.xticks(range(len(nodes_in)),nodes_in.keys())
plt.yticks(range(len(nodes_in)),nodes_in.keys())
#ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.imshow(np.transpose(CM[0.8]),vmax=3,vmin=0,cmap='Greys',origin='lower') # displays in color
plt.text(0.2,-0.3,'m=0.8')

###################################
lstyle={'Well-mixed & homogeneous':'-', 
       'Group-dependent & homogeneous':'--', 
       'Well-mixed & heterogeneous':'-.', 
       'Group-dependent & heterogeneous':':'}

# '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'

# calculation vs sim plots
ax=fig.add_subplot(gs[0:4,13:18])
plt.text(-0.23,0.86,'B',backgroundcolor='w',size=20, transform=ax.transAxes)


CV_index=0
CV=CV_values[CV_index]
plt.text(0,7.5,'h='+str(CV))
for calc in calcs:
    plt.plot(x_values,R[calc][:,CV_index],linewidth=2,linestyle=lstyle[calc],label=calc)
plt.scatter(x_values,R['Simulation'][:,CV_index],c='k',label='Simulation')
plt.xlabel('Group-dependent mixing, m')
plt.ylabel('$R_{0}$',labelpad=10,rotation='horizontal')
plt.ylim([3.5,8.3])

ax=fig.add_subplot(gs[5:9,13:18])
plt.text(-0.23,0.86,'D',backgroundcolor='w',size=20, transform=ax.transAxes)


x_index=0
x=x_values[x_index]
plt.text(0,7.5,'m='+str(x))
for calc in calcs:
    plt.plot(CV_values,R[calc][x_index,:],linewidth=2,linestyle=lstyle[calc],label=calc)
plt.scatter(x_values,R['Simulation'][x_index,:],c='k',label='Simulation')
plt.xlabel('Contact heterogeneity, h')
plt.ylabel('$R_{0}$',labelpad=10,rotation='horizontal')
plt.legend(title='$\\bf{Model~assumptions:}$                         ',bbox_to_anchor=(1.5, -0.3))
plt.ylim([3.5,8.3])


ax=fig.add_subplot(gs[0:4,19:24])
plt.text(-0.23,0.86,'C',backgroundcolor='w',size=20, transform=ax.transAxes)


CV_index=5
CV=CV_values[CV_index]
plt.text(0,7.5,'h='+str(CV))
for calc in calcs:
    plt.plot(x_values,R[calc][:,CV_index],linewidth=2,linestyle=lstyle[calc],label=calc)
plt.scatter(x_values,R['Simulation'][:,CV_index],c='k',label='Simulation')
plt.xlabel('Group-dependent mixing, m')
plt.ylabel('$R_{0}$',labelpad=10,rotation='horizontal')
plt.ylim([3.5,8.3])

ax=fig.add_subplot(gs[5:9,19:24])
plt.text(-0.23,0.86,'E',backgroundcolor='w',size=20, transform=ax.transAxes)


x_index=5
x=x_values[x_index]
plt.text(0,7.5,'m='+str(x))
for calc in calcs:
    plt.plot(CV_values,R[calc][x_index,:],linewidth=2,linestyle=lstyle[calc],label=calc)
plt.scatter(x_values,R['Simulation'][x_index,:],c='k',label='Simulation')
plt.xlabel('Contact heterogeneity, h')
plt.ylabel('$R_{0}$',labelpad=10,rotation='horizontal')
plt.ylim([3.5,8.3])


plt.savefig('../figures/figure1.png',format='png',dpi=300,bbox_inches='tight')
    



