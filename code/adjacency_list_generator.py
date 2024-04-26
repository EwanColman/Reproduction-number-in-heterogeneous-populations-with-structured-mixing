import pickle as pk
import numpy as np
from random import shuffle
from numpy.random import choice
import matplotlib.pyplot as plt
from numpy import linalg as LA


def get_adjacency_dictionary(N,E,x,CV):
    #groups=['A','B']
    #coefficient_of_variation={'A':1,'B':2}
    # store the number of nodes in the network in each age group
    
    # same number of nodes in each group
    # number_of_nodes={'A':int(N/2),'B':int(N/2)}
    
    # max value is 1+sqrt(2)
    # x=(eigenvalue*(eigenvalue-2))**(1/2)
    # print('x =',x)
    # print()
    
    
    number_of_edges={('A','A'):round(E*(1-x)/3),
                  ('A','B'):round(E*(1-x)/3),
                  #('B','A'):c,
                  ('B','B'):round(E*(2*x+1)/3)
                  }

    # to make degree equal
    #number_of_nodes={'A':N*(2*number_of_edges[('A','A')]+number_of_edges[('A','B')])/(2*E),
    #                 'B':N*(2*number_of_edges[('B','B')]+number_of_edges[('A','B')])/(2*E)}


    # mean={'A':(2*number_of_edges[('A','A')]+number_of_edges[('A','B')])/number_of_nodes['A'],
    #       'B':(2*number_of_edges[('B','B')]+number_of_edges[('A','B')])/number_of_nodes['B']}
    
    # print('mean A =',mean['A'])
    # print('mean B =',mean['B'])
    
    new_groups=['Ak','Al','Bk','Bl']
    new_number_of_nodes={'Ak':int(N/4),'Al':int(N/4),'Bk':int(N/4),'Bl':int(N/4)}
    
    # new_number_of_nodes={'Ak':round(number_of_nodes['A']/2),
    #                      'Al':round(number_of_nodes['A']/2),
    #                      'Bk':round(number_of_nodes['B']/2),
    #                      'Bl':round(number_of_nodes['B']/2)}
    
    
    
    # new_number_of_edges={('Ak','Ak'):number_of_edges[('A','A')]*(mean['A']+d)*(mean['A']+d)/(4*(mean['A']*mean['A'])),
    #                      ('Ak','Al'):number_of_edges[('A','A')]*(mean['A']+d)*(mean['A']-d)/(2*(mean['A']*mean['A'])),
    #                      #('Al','Ak'):number_of_edges[('A','A')]*(mean['A']+d)*(mean['A']-d)/(4*(mean['A']*mean['A'])),
    #                      ('Al','Al'):number_of_edges[('A','A')]*(mean['A']-d)*(mean['A']-d)/(4*(mean['A']*mean['A'])),
                         
    #                      ('Ak','Bk'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
    #                      ('Ak','Bl'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
    #                      ('Al','Bk'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
    #                      ('Al','Bl'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
                         
    #                      #('Bk','Ak'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
    #                      #('Bk','Al'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
    #                      #('Bl','Al'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
    #                      #('Bl','Ak'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
                         
    #                      ('Bk','Bk'):number_of_edges[('B','B')]*(mean['B']+d)*(mean['B']+d)/(4*(mean['B']*mean['B'])),
    #                      ('Bk','Bl'):number_of_edges[('B','B')]*(mean['B']+d)*(mean['B']-d)/(2*(mean['B']*mean['B'])),
    #                      #('Bl','Bk'):number_of_edges[('B','B')]*(mean['B']-d)*(mean['B']+d)/(4*(mean['B']*mean['B'])),
    #                      ('Bl','Bl'):number_of_edges[('B','B')]*(mean['B']-d)*(mean['B']-d)/(4*(mean['B']*mean['B']))
                         
    #                      }
    
    new_number_of_edges={('Ak','Ak'):number_of_edges[('A','A')]*(1+CV)*(1+CV)/4,
                         ('Ak','Al'):number_of_edges[('A','A')]*(1+CV)*(1-CV)/2,
                         #('Al','Ak'):number_of_edges[('A','A')]*(mean['A']+d)*(mean['A']-d)/(4*(mean['A']*mean['A'])),
                         ('Al','Al'):number_of_edges[('A','A')]*(1-CV)*(1-CV)/4,
                         
                         ('Ak','Bk'):number_of_edges[('A','B')]*(1+CV)*(1+CV)/4,
                         ('Ak','Bl'):number_of_edges[('A','B')]*(1+CV)*(1-CV)/4,
                         ('Al','Bk'):number_of_edges[('A','B')]*(1-CV)*(1+CV)/4,
                         ('Al','Bl'):number_of_edges[('A','B')]*(1-CV)*(1-CV)/4,
                         
                         #('Bk','Ak'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
                         #('Bk','Al'):number_of_edges[('A','B')]*(mean['A']+d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
                         #('Bl','Al'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']-d)/(4*(mean['A']*mean['B'])),
                         #('Bl','Ak'):number_of_edges[('A','B')]*(mean['A']-d)*(mean['B']+d)/(4*(mean['A']*mean['B'])),
                         
                         ('Bk','Bk'):number_of_edges[('B','B')]*(1+CV)*(1+CV)/4,
                         ('Bk','Bl'):number_of_edges[('B','B')]*(1+CV)*(1-CV)/2,
                         #('Bl','Bk'):number_of_edges[('B','B')]*(mean['B']-d)*(mean['B']+d)/(4*(mean['B']*mean['B'])),
                         ('Bl','Bl'):number_of_edges[('B','B')]*(1-CV)*(1-CV)/4
                         
                         }
    
    
    #print()
    for pair in new_number_of_edges:
        new_number_of_edges[pair]=round(new_number_of_edges[pair])
        #print(pair,new_number_of_edges[pair])
     
    #print('Total edges',sum(new_number_of_edges.values()))
    
    # this dictionary will store the node names (integers) of the nodes in each group    
    nodes_in={}
    # name the nodes
    index=0
    for group in new_groups:
        new_index=index+new_number_of_nodes[group]
        nodes_in[group]=[i for i in range(index,new_index)]
        index=new_index
    
    # and a dictionary to stor the neightbours of each node
    neighbours={}
    # make the dic
    degree={}
    for node in range(new_index):
        neighbours[node]=[]
        degree[node]=0
    
    
    # for each pair create a subnetwork that matches the data 
    for pair in new_number_of_edges:
        
        source_nodes=nodes_in[pair[0]].copy()
        target_nodes=nodes_in[pair[1]].copy()
        
        shuffle(source_nodes)
        shuffle(target_nodes)
        
        # spread the edges as evenly as possible
        for i in range(new_number_of_edges[pair]):
            # add edge choosing nodes sequentially
            source=source_nodes[i%len(source_nodes)]
            # use i+1 to avoid self-loops
            target=target_nodes[i%len(target_nodes)]
            
            #edge_list.append((source,target))
            neighbours[source].append(target)
            neighbours[target].append(source)


       
    return neighbours,{'A':nodes_in['Ak']+nodes_in['Al'],'B':nodes_in['Bk']+nodes_in['Bl']}

# neighbours,nodes_in=get_adjacency_dictionary(10000,100000,0.8,0.8)


# degree={}
# for node in neighbours:
#     degree[node]=len(neighbours[node])


# ### HISTOGRAM the degree distribution ######
# plt.figure()
# degree_dist=[0 for i in range(1+max(degree.values()))]
# for d in degree.values():
#     degree_dist[d]=degree_dist[d]+1/sum(degree.values())
# plt.bar(range(1+max(degree.values())),degree_dist,width=1)
# plt.yscale('log')
# plt.ylim([1/100000,1])
# #plt.xlim([0,20])
# plt.xlabel('Number of contacts')

# #######################

# # make the contact matrix C
# C=[]
# for a in nodes_in:
    
#     c=[]
#     for b in nodes_in:
#         # slice the a contacts on age of contact
#         count=0
#         for node in nodes_in[a]:
#             for neighbour in neighbours[node]:
#                 if neighbour in nodes_in[b]:
#                     count=count+1
#         if a==b:
#             count=count/2
#         c.append(count)

#         print(a+' to '+b+',',count,'contacts')        
#     C.append(c)
# C=np.array(C)




# c=sum(C[:])#+np.array([C[i,i] for i in range(len(C))])

# C=C/100000
# print(C)
# eigenvalues,eigenvectors=LA.eig(C)
# print(eigenvalues)
# print('Max eigenvalue of C =',max([np.real(z) for z in eigenvalues]))

# #### plotting ####
# fs=12
# fig=plt.figure()
# ax = fig.add_subplot()
# plt.xticks(range(len(nodes_in)),nodes_in.keys(),size=fs,rotation=60)
# plt.yticks(range(len(nodes_in)),nodes_in.keys(),size=fs)
# #ax.xaxis.tick_top()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.imshow(np.transpose(C),cmap='Greys',origin='lower') # displays in color
# plt.xlabel('Age of participant',size=fs)
# plt.ylabel('Age of contact',size=fs)
