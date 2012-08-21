#!/usr/bin/env python

import pickle
from numpy import *
import copy
import random


def hier_cluster( input = 'titles_and_mcssscorearray_charges.pickle', output= 'clusters.pickle'):

    """Take a pickle file containing a matrix of MCSS scores which is the similarities between each molecules and group molecules to clusters based on the scores.
    

    ARGUMENTS:
        - input: name of input pickle file;  default titles_and_mcssscorearray_charges.pickle
        - output: name of output pickle file; default clusters.pickle 

"""

#parameter
    offline = 0.8 
    small_cluster_size = 5 
    not_large_cluster_size = 1000 
    large_cluster_size = 1000

#    pickle_file = titles_and_mcssscorearray_charges.pickle

#Load MCSS scores
#    file = open( 'titles_and_mcssscorearray_charges.pickle' , 'r')
    file = open( input , 'r')
    (titles, scores, natoms) = pickle.load(file)
#    print titles

    file.close()
#Make a backup copy of the score array to edit destructively in what follows
    modified_scores = copy.copy(scores)
#get the similarity scores from the similar atom numbers by dividing by diagonal number which are the atom numbers of each molecule

#    print "Initial score matrix: ", scores
    Nelem = len(scores)
    max_num = []
    for i in range (Nelem):
        max_num.append(0)
        for k in range (Nelem):
            if modified_scores[i,k] > max_num[i]:
                max_num[i] = modified_scores[i,k]
            
    for i in range (Nelem):
        for j in range (Nelem):
            modified_scores[i,j] = max( modified_scores[i,j], modified_scores[j,i] )
########sliu 04/05/12
            if max_num[i] != 0:
                modified_scores[i,j] = modified_scores[i,j]/max_num[i]
#############
    
    #Once normalized, zero the diagonal elements since we don't want to consider self similarity (always 1)
        modified_scores[i,i] = 0   
#    print "Score matrix after modification: ", modified_scores

    def get_cluster_lists( clusters_by_molecule):
        
        
    #Obtain unique list of cluster numbers
        clusternumbers = unique(clusters_by_molecule.values())
    #Start cluster list
        clusters = {}
    #Build list
        for clusternumber in clusternumbers:
            clusters[clusternumber] = []
            for mol in clusters_by_molecule.keys():
                if clusters_by_molecule[mol] == clusternumber:
                    clusters[ clusternumber ].append( mol )
        return clusters
# get the largest score in each line if there are some same large scores pick one randomly

    clusters_by_molecule = {}
    cluster_num = 1
    largest_cluster_size = 0
#Build some initial clusters put one with the one which most similar to into the initial cluster

#line by line fine the largest number in each line if have several same largest scores randomly pick one
    for i in range (Nelem):
        largest = modified_scores[i,:].max()
#sliu debug 12/19/11#####
        if largest == 0:
            clusters_by_molecule [i] = cluster_num
            cluster_num += 1
#########
        else: 
            indices = where( modified_scores[i,:] == largest )
            index_array = arange (indices[0].size)
            random.shuffle (index_array)
            largest = indices[0][index_array[0]]
            #print modified_scores[i, largest]

            #raw_input()

            #Zero the element we're combining so it won't be looked at again
            modified_scores[i, largest] = 0. 
# put these two molecules into one cluster, if these two molecules do not in the exist build a new cluster
            if clusters_by_molecule.has_key( i ):
                if not clusters_by_molecule.has_key (largest):
                    clusters_by_molecule [largest] = clusters_by_molecule [i]
                #If the largest IS already in a cluster, nothing is done since both are already in clusters
            elif clusters_by_molecule.has_key( largest ):
                if not clusters_by_molecule.has_key (i):
                    clusters_by_molecule [i] = clusters_by_molecule [largest]

            else:
                clusters_by_molecule [i] = cluster_num
                clusters_by_molecule [largest] = cluster_num
                cluster_num +=1
    
    clusters = get_cluster_lists( clusters_by_molecule)
    print "initial cluster ", clusters
    l_cluster_sizea = []
    for j in clusters.values():
        t = len (j)
        l_cluster_sizea.append(t)
    print l_cluster_sizea
        
    #print clusters_by_molecule
    #raw_input()

#######merge clusters
    ct=0
    while True:
        ct+=1
        if ct%1000==0:
            print "Merging clusters; iteration %s..." % ct
        largest = modified_scores.max()

#################sliu 12/21/11        
        if largest < offline:
            cluster_num = len(clusters)                                                                                                                                                     
            cluster_size = []                                                                           
            clusternumbers = unique( clusters_by_molecule.values() )                                    
            clusterkeys = clusters.keys()                                                                                                                                                   
            clusternumbers.sort()                                                                       
            clusterkeys.sort()                                                                          
            for i in clusternumbers:                                                                    
                cluster_size.append ( len(clusters[i]) )                                           
            cluster_size_in_order = sorted (cluster_size)                                          
            largest_cluster_size = cluster_size_in_order[-1]                         
            break                                                                    
                                                                                     
###############
    #Find the best matches in our array
        for i in range (Nelem):
            line_indices = where( modified_scores[i,:] == largest )
        indices = where( modified_scores == largest )
    #Pick a random pair of molecules with the highest similarity to connect
        randentry = random.randint( 0, len(indices[0])-1 )
        largest = modified_scores[indices[0][ randentry],indices[1][randentry]]

        ct=0
        while clusters_by_molecule[ indices[0][ randentry] ] == clusters_by_molecule[ indices[1][randentry] ]: #If the molecules are in the same cluster, pick different ones
            ct+=1
            if ct%1000==0:
                print "Picking alternate cluster, iteration %s..." % ct
            modified_scores [indices[0][ randentry],indices[1][randentry]] = 0 
            largest = modified_scores.max()
            indices = where( modified_scores == largest )
            randentry = random.randint( 0, len(indices[0])-1 )
            largest = modified_scores[indices[0][ randentry],indices[1][randentry]]

    #check the score we get. If the score less than the offline then quit. These scores corresponding to the molecules from different clusters. 
        if largest < offline:
            #TO DO: Why do we have to do anything other than just break, here? (DLM asking...)
            cluster_num = len(clusters)                                                                                                                                                
 
            cluster_size = []
            clusternumbers = unique( clusters_by_molecule.values() )
            clusterkeys = clusters.keys()                                                                                                                                              
 
            clusternumbers.sort()
            clusterkeys.sort()
            for i in clusternumbers:
                cluster_size.append ( len(clusters[i]) )
            cluster_size_in_order = sorted (cluster_size)                                                                                                            
            largest_cluster_size = cluster_size_in_order[-1]                                                                                                         
            break

    #Identify what clusters the molecules are in
        cluster_1 = clusters_by_molecule [indices[0][ randentry]]
        cluster_2 = clusters_by_molecule [indices[1][ randentry]]
    # make sure cluster_1 <= cluster_2
        if cluster_1 > cluster_2:
            a = cluster_1
            cluster_1 = cluster_2
            cluster_2 = a
        if len(clusters[cluster_1])+len(clusters[cluster_2]) > large_cluster_size: 
            cluster_num = len(clusters)
            cluster_size = []
            clusternumbers = unique( clusters_by_molecule.values() )
            clusterkeys = clusters.keys()
            clusternumbers.sort()
            clusterkeys.sort()
            for i in clusternumbers:
                cluster_size.append ( len(clusters[i]) )
            cluster_size_in_order = sorted (cluster_size)
            largest_cluster_size = cluster_size_in_order[-1]    
            break
        for mol in clusters[cluster_2]:
            clusters_by_molecule[mol] = cluster_1
    #Update cluster list by removing cluster_2 and combining it with cluster 1
#        print clusters[cluster_1]
#        print clusters[cluster_2]
        clusters[ cluster_1 ] = clusters[ cluster_1] + clusters[ cluster_2 ]
        clusters.pop( cluster_2 )    
        cluster_num = len(clusters)
#get the largest cluster size
        cluster_size = []
        clusternumbers = unique( clusters_by_molecule.values() )
        clusterkeys = clusters.keys()
        clusternumbers.sort()
        clusterkeys.sort()
        for i in clusternumbers:
            cluster_size.append ( len(clusters[i]) )
        cluster_size_in_order = sorted (cluster_size)
        largest_cluster_size = cluster_size_in_order[-1]
        #TO DO ABOVE: Continue on even after the largest cluster passes the maximum cluster size, merging smaller clusters which are not past the largest cluster size. Could probbaly be done by zeroing the right elements in the score matrix. 


#merge small clusters with others which is not big
#get small clusters and not large clusters
    ct=0
    while True:
        ct+=1
        if ct%1000==0: print "Getting small clusters, iteration %s..." % ct
        small_clusters_name = []
        small_clusters = {}
        for i in clusternumbers:
            if len (clusters[i])< small_cluster_size:
                small_clusters_name.append(i)
                small_clusters [i] = []
                for mol in clusters[i]:
                    small_clusters [i].append (mol)
        
        if len(small_clusters_name) == 0: break

        not_large_clusters_name = []
        not_large_clusters = {}
        for i in clusternumbers:
            if len (clusters[i])< not_large_cluster_size:
                not_large_clusters_name.append(i)
                not_large_clusters [i] = []
                for mol in clusters[i]:
                    not_large_clusters [i].append(mol)

        #print not_larger_clusters
        print "Small cluster", small_clusters
        while True:
# get the possible scores from molecules pairs in clusters
            pairs_scores = []
            for i in small_clusters.keys():
                for j in not_large_clusters.keys():
                    for k in small_clusters[i]:
                        for l in not_large_clusters[j]:
                            pairs_scores.append (modified_scores[k,l])

            max_scores = 0.0

            for num in pairs_scores:
                if num > max_scores:
                    max_scores = num

            print max_scores

#get pairs' position of the max scores

#######sliu 12/22/11
            max_position = []

            for i in small_clusters.keys():
                for j in not_large_clusters.keys():
                    for k in small_clusters[i]:
                        for l in not_large_clusters[j]:
                            if max_scores == modified_scores[k,l]:
                                max_position.append([k,l])

            max_position = random.choice(max_position)
            modified_scores[max_position[0],max_position[1]] = 0
            
            for i in not_large_clusters.keys():
                if max_position[0] in not_large_clusters[i]:
                    first_cluster = i
                if max_position[1] in not_large_clusters[i]:
                    second_cluster = i
            
            if first_cluster != second_cluster:break
#have find a ideal pair and begin to test the score
        if max_scores < offline: break
##if the score is ideal try to merge the cluster

#        print max_position

        cluster_a = first_cluster
        cluster_b = second_cluster

        if cluster_a > cluster_b:
            trans = cluster_a
            cluster_a = cluster_b
            cluster_b = trans

        for mol in clusters[cluster_b]:
            clusters_by_molecule[mol] = cluster_a
        clusters[ cluster_a ] = clusters[ cluster_a] + clusters[ cluster_b ]
        clusters.pop( cluster_b )
        cluster_num = len(clusters)
        clusternumbers = unique( clusters_by_molecule.values() )

    #print clusters
    l_cluster_sizeb = []
    for j in clusters.values():
        t = len (j)
        l_cluster_sizeb.append(t)
    print l_cluster_sizeb
#    print cluster_num
#    print largest_cluster_size
    f = open ( output ,'w')
    pickle.dump(clusters,f)
    f.close


#If executing as main, take command-line arguments to drive functions
if __name__ =='__main__':

    from optparse import OptionParser
    parser = OptionParser()
    
    parser.add_option('-i', '--input', help='Input pickle file name', dest = 'input',default = 'titles_and_mcssscorearray_charges.pickle')
    parser.add_option('-o', '--output', help='output pickle file name', dest = 'output',default = 'clusters.pickle')
    
    (options, args) = parser.parse_args()

    # Invoke function to run hier_cluster
 
    hier_cluster( options.input, options.output) 
