import pickle as pk
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA
from random import choice
from adjacency_list_generator import get_adjacency_dictionary


runs=10
beta=0.2

x_values=[i/10 for i in range(10)]
CV_values=[i/10 for i in range(10)]
output={}
for x in x_values:
    for CV in CV_values:
        print(x,CV)
        # store R0 values in here
        output[(x,CV)]={}
        # generate a network
        neighbours,nodes_in=get_adjacency_dictionary(10000,100000,x,CV)
        
        R0={}
        converged_R0_list=[]
        R0_series_list=[]
        
        for run in range(runs):
            #print(run)
            seed=choice(list(neighbours.keys()))
            
            R0_series=[]
            ### code to make a tree ### move to separate file eventually            
            infected=[seed]                
            # loop over a fixed number of generations
            for generation in range(7):
                #print(run,generation)
                next_infections=[]        
                # for each node decide which neighbours to infect
                for node in infected:
                    # go through neighbours
                    for neighbour in neighbours[node]:
                        r=np.random.random()
                        if r<beta: #and neighbour!=infector[node]:
                            next_infections.append(neighbour)
                            #infector[neighbour]=node
                
                R0_series.append(len(next_infections)/max(1,len(infected)))
                infected=next_infections
            
            # store for output
            R0_series_list.append(R0_series)
            
            # if R0_series[-1]>0:
            #     converged_R0_list.append(R0_series[-1])
        
            #plt.plot(R0_series,'k',alpha=0.1)
        
        output[(x,CV)]['R0 by generation']=R0_series_list
        
        # print()    
        # print('R0=',np.mean(converged_R0_list))
        
        #plt.savefig('meeting_figs_feb20/R_by_generation2.png',format='png',dpi=300,bbox_inches='tight')
        
        # calculate it based on the data
        
        # get n_k,a for and c_ab from the network
        
        
        age_group_of={}
        for age_group in nodes_in:    
            #print(age_group,len(nodes_in[age_group]))
            for node in nodes_in[age_group]:
                age_group_of[node]=age_group
        
        number_of_contacts={}
        for a in nodes_in:
            for b in nodes_in:
                number_of_contacts[(a,b)]=0
        
        
        for node in neighbours:
            a=age_group_of[node]
            c=[]
            for neighbour in neighbours[node]:
                b=age_group_of[neighbour]
                if a==b:
                    number_of_contacts[(a,b)]=number_of_contacts[(a,b)]+1.0
                else:
                    number_of_contacts[(a,b)]=number_of_contacts[(a,b)]+1.0
        
        C=[]
        for a in nodes_in:
            c=[]
            for b in nodes_in:
                c.append(number_of_contacts[(a,b)])
            C.append(c)
        C=np.array(C)
                
        index_of={}
        i=0
        for age_group in nodes_in:
            index_of[age_group]=i
            i=i+1

        # make n - the arrany of the number of part_ids with 1, 2, etc. contacts in each age band
        K=10000
        n=[[0 for a in nodes_in] for k in range(K)]
        
        for node in neighbours:
                
            k=len(neighbours[node])
            # cap it (should we do this?)
            if k<K:
                n[k][index_of[age_group_of[node]]]=n[k][index_of[age_group_of[node]]]+1
            else:
                print('Omitted:',node,'degree =',k)
        n=np.array(n)
        
        #################################################
        #################################################
        #print()
        for group in ['A','B']:
            degree_list=[len(neighbours[node]) for node in nodes_in[group]]
            mean_degree=np.mean(degree_list)
            std_degree=np.std(degree_list)
            print(len(nodes_in[group]),group,CV,std_degree/mean_degree)
        
        degree_list=[len(neighbours[node]) for node in neighbours]
        mean_degree=np.mean(degree_list)
        
        #print('Mean * SAR', beta*mean_degree)
        output[(x,CV)]['Well-mixed & homogeneous']=beta*mean_degree
        
        ages=len(C)
        N=sum(n[:])
        c=sum(C[:])#+np.array([C[i,i] for i in range(len(C))])
        
        #eigenvalues,eigenvectors=LA.eig(C)
        #print('Max eigenvalue of C =',3*max([np.real(z) for z in eigenvalues])/sum(c))
        
        
        #### age structure, no heterogeneity
        
        # using basic method    
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
        
        #print(R0)
        #print('Structured age mixing, no heterogeneity:',R0)
        output[(x,CV)]['Group-dependent & homogeneous']=R0
        #### Heterogeneity, no age structure
        #print('Other old approach')
        
        R0=beta*(1/sum(c))*sum([sum([(k**2)*n[k,b] for k in range(0,K)]) for b in range(ages)])
        #t2=t2+time()-start
        #print('Heterogeneity & proportionate age mixing:',R0)
        output[(x,CV)]['Well-mixed & heterogeneous']=R0
        
        
        ### Combined method
        M=[]
        for a in range(ages):
            m=[]
            for b in range(ages):
                m.append(C[a,b]/(c[a]*c[b]))
        
            
            M.append(m)
            
        M=np.array(M)
        #print('Analytical G')
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
        #print('R0=')
        R0=max([np.real(z) for z in eigenvalues])
        #t3=t3+time()-start
        #print('Combined method:',R0)
        output[(x,CV)]['Group-dependent & heterogeneous']=R0
        
        # for plot
pk.dump(output,open('../pickles/R0_analytical_and_simulated.p','wb'))



