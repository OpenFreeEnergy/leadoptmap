#!/usr/bin/env python
import glob
import sys
import pickle
#from mmtools.moltools.ligandtools import *
#from openeye.oechem import *
from numpy import *
import copy
import networkx as nx
from schrodinger import structure


"""HANDLE VISUALIZATION OF CLUSTERS AND CALCULATION PLANNING"""


#=======================================================================
#INPUT
#=======================================================================

#Name of pickle file containing output of MCSS hierarchical clustering (that is to say, clusters we want to use as our starting point)
clusterfile = 'cluster_list.pickle' 
#ALTERNATIVELY: Just use each charge as a cluster by setting that to None
#clusterfile = None

#Name of pickle file containing titles and score array used for MCSS clustering
scorematrix = 'titles_and_mcssscorearray_charges.pickle'

#Score matrix from 'unstrict' checking to use rescuing orphaned clusters
unstrict_scorematrix = 'titles_and_mcssscorearray_charges_unstrict.pickle'

#Maximum distance (number of calculations/edges) allowed between nodes in a cluster. Setting this smaller will result in extra calculations being added within clusters to ensure that all nodes are with MAXDIST of one another
MAXDIST = 6

#Minimum similarity (common atoms) for connections between nodes
MINSIM = 5

#Print debugging information?
debug = True

#Intra-cluster weight (use this to keep individual clusters together)
intraweight = 1

#Minimum and maximum 'starting' cluster sizes -- the clusters coming out of hierarchical clustering aren't necessarily OK in size due to a balancing act we had to play in the clustering. That is, we had to settle with too many clusters some of which are way too small, to avoid ending up with clusters which are way too big. Here, as an initial step, we will merge some clusters. Set minimum and maximum cluster sizes for merging
minclustersize = 5 #Try to make every cluster have at least 5 if possible; soft cutoff (if similarity is low, clusters of size smaller than this will be tolerated)
maxclustersize = 1000 #No clusters larger than 45 molecules allowed


#=======================================================================
#LOAD DATA FILES AND MANAGE SOME DATA
#=======================================================================

#Load pickle file containing list of clusters and names in each cluster
if clusterfile: 
    file = open(clusterfile, 'r')
    clusterlist = pickle.load(file)
    file.close()
    print clusterlist #DEBUG
else: #Otherwise, we're just clustering by charge
    file = open('chargedict.pickle', 'r')
    charges = pickle.load(file)
    file.close()
    #Build a dictionary of molecules having each charge
    mols_by_chg = {}
    for mol in charges.keys():
        if mols_by_chg.has_key( charges[mol] ):
            mols_by_chg[ charges[mol] ].append(mol)
        else:
            mols_by_chg[ charges[mol] ] = [mol]
    clusterlist = [  mols_by_chg[chg] for chg in mols_by_chg ]

#Load titles and scores from MCSS structure comparison (score matrix modified to ensure that molecules of different charge end up in different clusters)
file = open(scorematrix, 'r')
(titles, mcss_scores, number_of_atoms) = pickle.load(file)
file.close()
print titles #DEBUG

#Load unstrict scores
file = open(unstrict_scorematrix, 'r')
(titles, mcss_scores_unstrict, dum) = pickle.load(file)
file.close()

#Build a dictionary to allow lookup of molecule size by title and molecule index by title
natoms = {}
molidx = {}
for (idx, mol) in enumerate(titles):
    natoms[mol] = number_of_atoms[idx]
    molidx[mol] = idx

#Number of molecules
Nmols = len(titles)

#Re-fill symmetric elements in score array
for i in range(len(titles)):
    for j in range(i):
        maxnonzero = max(mcss_scores[i,j], mcss_scores[j,i])
        maxnonzero_unstrict = max(mcss_scores_unstrict[i,j], mcss_scores_unstrict[j,i])
        mcss_scores[i,j] = maxnonzero
        mcss_scores[j,i] = maxnonzero
        mcss_scores_unstrict[i,j] = maxnonzero_unstrict
        mcss_scores_unstrict[j,i] = maxnonzero_unstrict




#=======================================================================
#OBTAIN SOME MOLECULAR PROPERTIES BY TITLE
#=======================================================================

#Read molecule file, obtain net charge
#Read molecule file, obtain net charge
#inputmolfiles = glob.glob('mol2_file/*.mol2')

#Read molecules into a dictionary by title; also store charges

file = open('chargedict_openeye.pickle', 'r')
charges = pickle.load(file)
file.close()

####### have different formal charge with openeye one. Here use old openeye one for charge 
#for molfile in inputmolfiles:
#    for mol1 in structure.StructureReader (molfile):
#    mol1 = structure.StructureReader (molfile)
#        title = mol1.title
#        charges[title] = mol1.formal_charge

#print charges
#file = open('chargedict.pickle', 'w')
#pickle.dump(charges, file)


#=======================================================================
#FUNCTION DEFINITIONS
#=======================================================================

def get_similarity( title1, title2, scores = mcss_scores):
    """Get similarity (number of atoms) between two titles."""
    idx1 = molidx[ title1]
    idx2 = molidx[ title2]
    return scores[ idx1, idx2 ]

def get_similarity_unstrict( title1, title2):
    """Get similarity (number of atoms) between two titles."""
    idx1 = molidx[ title1]
    idx2 = molidx[ title2]
    return mcss_scores_unstrict[ idx1, idx2 ]

def get_most_similar_clustermember( title, clus, debug = False, scores = mcss_scores):
    """Take a molecule title and 'clus', a list of (some of the) members of its cluster, and find the most similar cluster member, returning it. If the cluster only contains the original molecule, this is returned as the most similar member, otherwise another molecule is returned."""
    maxsim = 0
    besttitle = None
    idx1 = molidx[ title ]

    #If the cluster only has one entry which is self, no point in checking -- just return this
    if len(clus) < 2 and title in clus:
        return title

    for title2 in clus:
        if title != title2: #Only bother checking if we are not looking at the same molecule
            sim = scores[ idx1, molidx[ title2] ]
            if debug:
                print "    Similarity of %s-%s is %.2f..." % (title, title2, sim)
            if sim > maxsim:
                maxsim = sim
                besttitle = title2

    return besttitle
def get_most_similar_members_from_lists( clus1, clus2, debug = False, scores = mcss_scores):
    """Take lists of two sets of molecules and return the titles of the most similar pair.

    NOTE: Could be optimized by using array operations to do this, I think, but I'm going for now with the most straightforward implementation."""

    maxsim = 0
    besttitle1 = None
    besttitle2 = None

    for title1 in clus1:
        idx1 = molidx[ title1 ]
        for title2 in clus2:
            idx2 = molidx[ title2 ]
            sim = scores[ idx1, idx2]
            if sim > maxsim:
                maxsim = sim
                besttitle1 = title1
                besttitle2 = title2

    return besttitle1, besttitle2

def count_non_zero (score_dic):
    count = 0
    for i in score_dic.keys():
        if score_dic[i] > 0:
            count += 1
    return count

def incomplete_clusters_networkx(cluster,Gr):
    communicability = nx.communicability(Gr)
    inaccessible_cluster = []
    cluslen = len(cluster)
    
    accesslen = count_non_zero(communicability[cluster[0]])
    if accesslen < cluslen and cluslen > 1:
        inaccessible_cluster = cluster
    
    return inaccessible_cluster

def most_distant_members_networkx( cluster,  simcutoff = 0, scores= mcss_scores):
#def most_distant_members_networkx( cluster, Graph, simcutoff = 0):
    maxdist = 0
    for member in cluster:
        #Compute shortest path between this member node and all other nodes
        shortest_distance = nx.single_source_shortest_path_length(Gr, member)
        #Compute similarity between this node and the others and remove anything not meting our similarity cutoff from the shortest distance dictionary        
        sims = {}

        for member2 in cluster:
            if not member2 == member:
                sim = get_similarity( member, member2, scores = scores)
                if sim < simcutoff:
                    shortest_distance.pop( member2 )

        #Compute maximum distance between this node and any of the remaining elements meeting our similarity cutoff
        this_maxdist = max( shortest_distance.values() )
        #If that is larger than our existing max distance, find which nodes it's between
        if this_maxdist > maxdist:
            maxdist = this_maxdist #Store new maximum distance
            #Find node pair
            #print "check short dis", shortest_distance.keys()
            #raw_input()
            member2_list = []
            for member2 in shortest_distance.keys():
                if shortest_distance[member2] == maxdist:
                    member2_list.append(member2)
            member2_list.sort()
            title1 = member
            title2 = member2_list[0]
            if title1 == 'frag.vs.084':
    
                print "Title", title1, title2, maxdist, shortest_distance
    #Return the maximum distance and the cluster members
    return (title1, title2, maxdist)

def find_molecule_not_in_circle_networkx (cluster,Gr):
    cycle_list = nx.cycle_basis(Gr)
    node_in_cycle = []
    for cycle in cycle_list:
        for cycle_node in cycle:
            if cycle_node not in node_in_cycle:
                node_in_cycle.append(cycle_node)
#    for node_need_delete in node_in_cycle:
#        if node_need_delete in tem_cluster:
#            tem_cluster.remove(node_need_delete)
    single_node_list = []
    for single_node in cluster:
        if not single_node in node_in_cycle: 
            single_node_list.append(single_node)
    return single_node_list

#=======================================================================
#PLAN CALCULATIONS
#=======================================================================

#DESIGN CRITERIA:
# 0) Merge clusters which are too small, up to specified maximum cluster size (assuming similarity is high enough). 
# 1) Connect each node to the node in its cluster that it is most similar to as a starting point (which also ensures that each node is connected to at least one other node)
# 2) Every node in each cluster must be interconnected
# 3) Every cluster larger than 2 members should contain at least one cycle; also enforce a maximum number of calculations across a given cluster
# 4) Every cluster of the same net charge must be interconnected (preferably to another cluster it shares some similarity with)


lens = [ len(cluster) for cluster in clusterlist ]
print "Done consolidating clusters. Current cluster statistics: Max cluster size %s, min cluster size %s, and %s clusters:" % (max(lens), min(lens), len(clusterlist) )
print "Pausing..."
#raw_input()

#SET UP FOR STEP 1
#Storage for a dictionary of planned connections
connections = {}
for title in titles:
    connections[title]=[]


#Begin planning connections by connecting each node to the node in its cluster that it is most similar to
for clus in clusterlist:
    for title in clus:
        besttitle = get_most_similar_clustermember( title, clus)
        #If there is no best match, use unstrict checking
        if not besttitle:
            besttitle = get_most_similar_clustermember( title, clus, scores = mcss_scores_unstrict )
        if title!=besttitle: #Avoid self connections 
            connections[title].append(besttitle)
            connections[besttitle].append(title)

#Ensure uniqueness in lists
for title in titles:
    connections[title] = unique(connections[title])

#The remaining part of this will be done in pygraph after adding these connections

#=======================================================================
# VISUALIZE CLUSTERS AND PLAN INTRA-CLUSTER CALCULATIONS
#=======================================================================

#Obtain number of clusters
Nclus = len(clusterlist)
print "PROCESSING %s MOLECULES in %s CLUSTERS..." % (Nmols, Nclus)

#Obtain sizes of largest and smallest molecules (will be used in adjusting sizes of nodes)
maxatoms = max(number_of_atoms)
minatoms = min(number_of_atoms)

#We will color each cluster differently; color specifications are given here: http://graphviz.org/doc/info/colors.html
#Build a list of colors we will use for different clusters
color_list = [ 'antiquewhite', 'aquamarine', 'bisque4', 'blue', 'chartreuse', 'chocolate', 'crimson', 'cyan', 'darkgreen', 'darkorchid', 'deeppink', 'gold', 'maroon', 'palegreen', 'salmon', 'turquoise4', 'yellow', 'cornflowerblue', 'coral1', 'darkorange', 'darkslateblue', 'darkslategrey', 'magenta', 'orchid', 'purple', 'black', 'antiquewhite2', 'aquamarine3', 'azure3', 'bisque3', 'blue4', 'blueviolet', 'brown4', 'chartreuse4', 'cadetblue3', 'coral3', 'cornsilk3', 'cyan2', 'darkgoldenrod2', 'darkorchid3', 'darkseagreen2', 'darkslategray', 'deeppink4', 'deepskyblue3', 'dimgray', 'firebrick4', 'dodgerblue3', 'floralwhite', 'gainsboro', 'gray', 'gray100', 'greenyellow', 'hotpink3', 'indianred3', 'ivory2', 'khaki2', 'lavenderblush1', 'lemonchiffon', 'lightcyan2', 'lightgoldenrod4', 'lightpink2', 'lightskyblue3', 'magenta3', 'maroon3', 'mediumorchid2', 'mediumpurple2', 'mediumspringgreen', 'mistyrose1', 'navajowhite'] #Many more can be added if this is not enough; note colors have some similar repeats since usually each net charge does not have as many molecules

if len(color_list) < Nclus: print "WARNING: ADD MORE COLORS FOR CLUSTERS"


#Create empty graph
Gr = nx.Graph()

#Loop over clusters and build clusters
print "BUILDING INITIAL CLUSTERS..."
for (idx, clus) in enumerate(clusterlist):
    if debug: print "   Cluster %s of %s, size %s..." % (idx, len(clusterlist), len(clus))
    for elem in clus:
        Gr.add_node( elem )
#######need modify later       Gr.add_node_attribute( elem, ['color', color_list[idx] ] )
#######need modify later       Gr.add_node_attribute( elem, ['style', 'filled'] )
        #Calculate a size based on the number of atoms. Smallest molecule will have a size 1/(2*(maxatoms-minatoms+1)); largest molecule will have size 0.5

#Add edges based on nearest cluster memberships (#1 from above)
#Colorize edges added for this reason distinctively
print "ADDING INITIAL CLUSTER CONNECTIONS..."
for title in titles:
    for title2 in connections[title]:
        #Add connection
        if not Gr.has_edge(title, title2):
            Gr.add_edge( title, title2)

########check graph
file01 = open('after_step_one.pickle', 'w')
pickle.dump(Gr, file01)
file01.close()
 
#######need to change and check()            Gr.add_edge_attribute(title, title2, ['color', 'black']) #Green edges for those edges added due to greatest similarity

#=================================================================
# Step 2 of planning strategy: Every node of every cluster must be accessible from every other node within the cluster
# To do this:
# - Calculate accessibility (any node j you can get to from node i is flagged)
# - While not all members of a cluster are accessible:
#   - For each cluster with inaccessible members:
#     - Pick the least accessible member
#     - Connect it to another member of the cluster (the most similar) that is not already accessible to it

############debugcheck
#print "Neighbor of 018_1", Gr.neighbors('frag.vs.018_1'), len(clusterlist)
#for item in clusterlist:
#    if 'frag.vs.018_1' in item:
#        print item

#build subgragh 

Gr_sep = {}
for (idx,cluster) in enumerate(clusterlist):
    Gr_sep[idx] = Gr.subgraph(cluster)


for sub in Gr_sep.keys():
    #print sub
    clus = clusterlist[sub]

    #Determine imcomplete clusters (clusters with inaccessible members)
    incomplete = nx.is_connected( Gr_sep[sub])

    print "Check point one ", incomplete
    while incomplete is False:
        
        #Pick the least accessible member
        memberaccess = array([ len(nx.node_connected_component(Gr_sep[sub], title)) for title in clus ])
        minaccess = min(memberaccess)
        minaccess_idx = where( memberaccess == minaccess)
        mintitle = clus [ minaccess_idx[0][0] ] #Title of least connected member

        #Join it to the other member it is the most similar to that it is NOT already connected to
        #First build a list of what cluster members it is not connected to
        unconnected_clus = [ title for title in clus if not title in nx.node_connected_component(Gr_sep[sub], mintitle) ]
        newconnection = get_most_similar_clustermember( mintitle, unconnected_clus)
        #If there is no similar member according to this scoring, use the unstrict scoring
        if not newconnection:
            newconnection = get_most_similar_clustermember( mintitle, unconnected_clus, scores = mcss_scores_unstrict )

        #It is not always the case that the least connected member actually can connect to any others -- i.e. if there is a molecule in the cluster which only has similarity to one other molecule (but that other molecule is similar to several) the first molecule may be impossible to connect to any others. In such ases we need to back off and repeat the search using a better-connected member
        while not newconnection:
            print "Picking a better connected member to use in improving accessibility"
            tmpclus = copy.copy( clus )
            tmpclus.remove( mintitle )
            memberaccess = array([ len(nx.node_connected_component(Gr_sep[sub], title)) for title in tmpclus ])
            minaccess = min(memberaccess)
            mintitle = tmpclus[ minaccess_idx[0][0] ] #Title of least connected member
            unconnected_clus = [ title for title in tmpclus if not title in nx.node_connected_component(Gr_sep[sub], mintitle) ]
            newconnection = get_most_similar_clustermember( mintitle, unconnected_clus)
            if not newconnection:
                newconnection = get_most_similar_clustermember( mintitle, unconnected_clus, scores = mcss_scores_unstrict)

        #Add this connection
        print "FFFFFFFFind and add new edge %s %s" %(mintitle, newconnection)
        Gr_sep[sub].add_edge(mintitle, newconnection)
        Gr.add_edge( mintitle, newconnection) #Red edges for those edges added to flesh out cluster
#########need to change        Gr.add_edge_attribute( (mintitle, newconnection), ['color', 'red']) 

    #Recompute incomplete clusters
        incomplete = nx.is_connected( Gr_sep[sub])

########check graph
file02 = open('after_step_two.pickle', 'w')
pickle.dump(Gr, file02)
file02.close()

############debugcheck
#print "Neighbor of 018_1", Gr.neighbors('frag.vs.018_1'), len(clusterlist)
#for item in clusterlist:
#    if 'frag.vs.018_1' in item:
#        print item 
#=================================================================
# Step 3 of planning strategy: 
# (a) Every cluster should have no more than (MAXDIST) calculations to span across it -- that is to say, the maximum distance between any two members must be no more than MAXDIST
# (b) Every cluster member should be involved in at least one cycle (originally, I thought that I'd have every cluster contain at least one cycle, but this has the potential to leave lots of members stranded outside of cycles, though it's unlikely as the cluster csize increases. Insisting every member be in at least one cycle does not add that many additional connections (at least for larger clusters) but adds considerable robustness.) 


#Step 3a: Ensure each cluster has no more than MAXDIST calculations to span it
#To do so: Find the two elements which are the most distant, and if this distance exceeds MAXDIST, connect them IF they have at least MINSIM atoms in common, otherwise look for another match. Iterate until MAXDIST is satisfied.

new_edges_max = []
#Loop over clusters 
print "ENSURING EACH CLUSTER TAKES NO MORE THAN %s CALCULATIONS TO SPAN..." % MAXDIST
for (idx, clus) in enumerate(clusterlist):
    #Don't bother with single or two-member clusters since these are already fully connected
    if len(clus) > 2:
        #Find most distant cluster members 
        (title1, title2, maxdist) = most_distant_members_networkx( clus, simcutoff = MINSIM)
        if not title1 or not title2: #If there are no distant members, use alternative scoring
            print "FFFFFFFFFFFind use another score"
            (title1, title2, maxdist) = most_distant_members_networkx( clus, simcutoff = MINSIM, scores = mcss_scores_unstrict)
            
        
        #Iterate and join nodes
        iter = 0 #Iteration counter just to check we haven't gotten stuck
        while maxdist > MAXDIST and iter < Nmols:
            Gr.add_edge(title1, title2)
            new_edges_max.append((title1, title2, maxdist)) 
########need to change            Gr.add_edge_attribute( (title1, title2), ['color', 'blue']) #blue edges for those edges added to flesh out cluster

            #Recompute most distant members
            (title1, title2, maxdist) = most_distant_members_networkx( clus, simcutoff = MINSIM)
            if not title1 or not title2: #If there are no distant members, use alternative scoring
                (title1, title2, maxdist) = most_distant_members_networkx( clus, simcutoff = MINSIM,scores =  mcss_scores_unstrict)
            
            #NOTE: This could enter an infinite loop in some unusual cases, such as a cluster consisting of two large blocks, each of which has high internal similarity but none of the elements of block A has more than MINSIM atoms in common with any of the elements in block B. This is probably unlikely -- but check and print a warning if anything bad happens here
            iter+=1
            if iter > Nmols/2.:
                print "WARNING: At iteration %s within step 3a; possible cluster blocking problem..." % iter
        if iter >= Nmols:
            print 'WARNING: Terminated loop at 3a due to exceeding maximum iteration number; likely cluster blocking problem...' 

        if debug and iter>0: print "   For cluster %s, added %s connections due to MAXDIST" % (idx, iter)


print "New connection due to max dist %s" %len(new_edges_max), new_edges_max
f_max = open('max_dist.pickle','w')
pickle.dump(new_edges_max, f_max)
f_max.close()

f1 = open('graph_new_fuction.pickle','w')
pickle.dump(Gr,f1)
f1.close()

nx.draw_graphviz(Gr)
nx.write_dot(Gr,"before_cycle.dot")


############debugcheck
#print "Neighbor of 018_1", Gr.neighbors('frag.vs.018_1'), len(clusterlist)
#for item in clusterlist:
#    if 'frag.vs.018_1' in item:
#        print item

# Step 3b: Every cluster member should be in at least one cycle; if a cluster member is not,  create one by combining the two most distant elements with adequate similarity then check again. (NOTE: This could be streamlined by not doing this iteratively -- at the initial pass, just add one connection for every cluster member not already in a cycle).
print "ENSURING EVERY CLUSTER MEMBER IS IN AT LEAST ONE CYCLE..."

Gr_sep_cycle = {} 
for (idx, cluster) in enumerate(clusterlist):
    # build sub graph for cycle detection 
    Gr_sep_cycle[idx] = Gr.subgraph(cluster)

for sub in Gr_sep_cycle.keys():
    clus = clusterlist[sub]

#########debugcheck
#    if 'frag.vs.018_1' in clus:
        #print "ClusSSSSSSSSSSSSSSSSS ", clus, sub, len(Gr_sep_cycle.keys())
    #Don't bother with single or two-member clusters since cycles are meaningless there
    if len(clus) > 3:
######sliu 3/23/12
#    if len(clus) > 2:
##########
        #Obtain a list of those nodes occurring in only one cluster
        if debug: print "   Finding singleton nodes for cluster %s..." % idx
        singleton_nodes = find_molecule_not_in_circle_networkx (clus, Gr_sep_cycle[sub])
        if debug: print "   For cluster %s, there are initially %s singleton nodes which will be added to cycles..." % (idx, len(singleton_nodes))

        #Are we going to iteratively remove singleton nodes, or just do it all at once? Can get a big speed improvement by doing it all at once, but it will result in some extra connections. 
        #Iterative solution usually results in extremely slow clustering finalization and it not particularly recommended. 
        iteratively_solve = True 

        #As long as there are singleton nodes:
        if iteratively_solve:
          #Obtain a starting point for singleton nodes to work on
          worknode = singleton_nodes[0]
#          print "SSSSS" , singleton_nodes
#          print "WWWWW", worknode
          while singleton_nodes:
            
            #Build a list of cluster members it is not already connected to
            unconnected = copy.copy(clus)
#            print "UUUUUUUUU ", unconnected
#            print "Neighbors", Gr_sep_cycle[sub].neighbors( worknode )
            for nbr in Gr_sep_cycle[sub].neighbors( worknode ):
#                print "NNNNNNNNNN", nbr
########need to change


                if nbr in unconnected:
######
                    unconnected.remove( nbr )
            #Remove self from this list also -- we don't want to make self connections
            unconnected.remove( worknode )
            
           
############### sliu 3/22/12 in case unconnected is empty like 3 members which a-b a-c but b does not connect c, if worknode is a, it cause problem.
            
            if len(unconnected) > 0:
 
############################
            #Find most similar cluster member that is sufficiently distant which it is not already connected to
                title = get_most_similar_clustermember( worknode, unconnected)

            #Connect them IF we're not talking about a self-connection
            #Self connections can happen if we are looking at, for example, the middle member of a three-member cluster where the middle member is connnected to both other members -- the only unconnected member will be self.

################ sliu 03/22/12
                if not title:                                                                                           
                        title = get_most_similar_clustermember( worknode, unconnected, scores = mcss_scores_unstrict)

                if title and not worknode==title:

################
            #if not worknode==title:
                    Gr.add_edge(worknode, title)
                    Gr_sep_cycle[sub].add_edge(worknode, title)
                    print "ADD.............", worknode, title
############need to change                    gr.add_edge_attribute( (worknode, title), ['color', 'green']) 


            #Recompute singleton nodes
            singleton_nodes = find_molecule_not_in_circle_networkx (clus, Gr_sep_cycle[sub])

            #Obtain new work node
            if singleton_nodes:
                newworknode = singleton_nodes[0]
                if newworknode == worknode: #If we are stuck (i.e. this work node can't have new connections) go to the next one
                    if worknode in singleton_nodes: #If we are stuck on something that's still in our list 
                        worknode = singleton_nodes[ singleton_nodes.index(worknode)+1] #Look at the next entry in the list
                    else: #If we're not stuck on something still in our list, just use the first one in our list
                        worknode = singleton_nodes[0]
                else: #Otherwise go on to the next work node (the next singleton we found)
                    worknode = newworknode

        if not iteratively_solve:
            
            for worknode in singleton_nodes:
                if debug: print "      Singleton node %s..." % worknode
                #Build a list of cluster members it is not already connected to
                unconnected = copy.copy(clus)
############# 
                for nbr in Gr_sep_cycle[sub].neighbors(worknode):
#######need to change
                    #if nbr in unconnected:
############
                    unconnected.remove(nbr)

                #Remove self from this list also -- we don't want to make self connections
                unconnected.remove( worknode )
            
                #Find most similar cluster member that is sufficiently distant which it is not already connected to
                title = get_most_similar_clustermember( worknode, unconnected)
                if not title:
                    title = get_most_similar_clustermember( worknode, unconnected, scores = mcss_scores_unstrict)
                if debug: print "      For node %s, attempting to connect to member %s..." % (worknode, title)

                #Connect them IF we're not talking about a self-connection
                #Self connections can happen if we are looking at, for example, the middle member of a three-member cluster where the middle member is connnected to both other members -- the only unconnected member will be self.
                if title and not worknode==title:
                    #print "      Adding edge..."
                    Gr.add_edge( worknode, title)
                    Gr_sep_cycle[sub].add_edge(worknode, title)
                    #print "      Done adding edge, labeling..."
#######need to change                    gr.add_edge_attribute( (worknode, title), ['color', 'green']) 
                    #print "      Labeled...."

############debugcheck
#print "Neighbor of 018_1", Gr.neighbors('frag.vs.018_1'), len(clusterlist)
#for item in clusterlist:
#    if 'frag.vs.018_1' in item:
#        print item

print "Before graph............"
##########sliu 12/22/11 
#file = open('graph.pickle','w')
#pickle.dump(gr,file)
#file.close

########check graph
file03 = open('after_step_three.pickle', 'w')
pickle.dump(Gr, file03)
file03.close()

#WRITE OUT CLUSTERS ONCE BEFORE LINKING
nx.draw_graphviz(Gr)
nx.write_dot(Gr,"before_link.dot")
#==================================================================
# PLAN CALCULATIONS SPANNING ACROSS CLUSTERS AND UPDATE CLUSTER VISUALIZATION
#==================================================================
# Plan here:
# 4) Add relative calculations between each pair of clusters sharing the same net charge if they have members which are sufficiently similar
# 5) Add relative calculations to benzamidine if sufficiently similar and have same net charge
# 6) Track any clusters which don't get connected to others (for reasons such as not having enough similarity); possibly do dual topology calculations for these later.

#STEP 4: Link clusters of same net charge
#===================================================================


#Determine net charge of each cluster; store it in a list
cluster_charge = []
for cluster in clusterlist:
    cluster_charge.append(charges[ cluster[0]])

#For each cluster, link one member to the most similar member of each cluster having the same net charge
#First, build a dictionary of clusters having each charge
clusteridx_by_charge = {}
for (idx, cluster) in enumerate(clusterlist):
    chg = cluster_charge[idx]
    if not clusteridx_by_charge.has_key(chg):
        clusteridx_by_charge[chg] = [ idx ]
    else:
        clusteridx_by_charge[chg].append(idx)

#Now link each cluster to every other cluster sharing the same net charge, tracking a list of connections between clusters which have already been made to avoid duplication
connected_clusters = {}
for idx in range(len(clusterlist)):
    connected_clusters[idx]=[]

#Track edges added at this stage for weighting purposes
newedges=[]
print "LINKING CLUSTERS SHARING SAME NET CHARGE"
for chg in clusteridx_by_charge.keys():
    print "   Working on net charge %s..." % chg

    #Work on cluster to link
    for idx in clusteridx_by_charge[chg]:
        #Build list of other clusters to link it to
        other_clusters = copy.copy(clusteridx_by_charge[chg])
        other_clusters.remove(idx)
    
        #Loop over target clusters
        for targetclus in other_clusters:
            #Find a pair of molecules to (potentially) link
            title1, title2 = get_most_similar_members_from_lists( clusterlist[idx], clusterlist[ targetclus ])
            if not title1 or not title2:
                title1, title2 = get_most_similar_members_from_lists( clusterlist[idx], clusterlist[ targetclus], scores = mcss_scores_unstrict )

            #If the clusters are not already linked, and we get back a target set of molecules (which won't happen if the clusters are too disssimilar), link them
            if title1 and title2 and targetclus not in connected_clusters[ idx ]:
                #Add connection
                Gr.add_edge( title1, title2 )
########need to change                gr.add_edge_attribute( (title1, title2), ['color', 'magenta']) #Magenta edges for those added to link clusters
                #Track that we added it
                connected_clusters[idx].append( targetclus )
                connected_clusters[ targetclus].append(idx)
                newedges.append( (title1, title2) )
                newedges.append( (title2, title1))
                
                #NOTE: If there are more than two clusters, this will result in multiple paths between every pair of molecules, since every molecule is within at least one cycle (meaning that there are two paths between every pair of molecules within every cluster) and every cluster is connected to every other cluster (so there are multiple paths between each pair of clusters) 
                #However, if there are only two clusters, there will be only one path between the two clusters, and hence some molecules which are only connected by one path. So, in the special case of only two clusters, add an additional connection
                if len( clusteridx_by_charge[chg]) == 2:
                    unlinked_this = copy.copy( clusterlist[idx])
                    unlinked_other = copy.copy( clusterlist[ targetclus] )
                    unlinked_this.remove( title1 )
                    unlinked_other.remove( title2 )
                    ntitle1, ntitle2 = get_most_similar_members_from_lists( unlinked_this, unlinked_other)
                    if not ntitle1 or not ntitle2:
                        ntitle1, ntitle2 = get_most_similar_members_from_lists( unlinked_this, unlinked_other, scores = mcss_scores_unstrict)
                    if ntitle1 and ntitle2:
                        #Add connection
                        Gr.add_edge( ntitle1, ntitle2 )
#######need to change                        gr.add_edge_attribute( (ntitle1, ntitle2), ['color', 'magenta'])
                        #Track that we added it
                        connected_clusters[idx].append( targetclus )
                        connected_clusters[ targetclus].append(idx)
                        newedges.append( (ntitle1, ntitle2) )
                        newedges.append( (ntitle2, ntitle1) )

#==================================================================
#WRITE OUT CLUSTERS
#==================================================================

#PRINT SOME STATISTICS: Total number of nodes, total number of edges
print "\n"
edges = Gr.edges() #Note that having one edge connecting two nodes actually shows up in this list as two edges, i.e. an edge from 3-4 and an edge from 4-3 (since undirected). So for total number of edges we divide by 2
edgenum = len(edges)
molnum = len( Gr.nodes() )
print "%s molecules in %s clusters with %s edges, having charges ranging from %.2f to %.2f" % (molnum, len(clusterlist), edgenum, min(charges.values()), max( charges.values()))

print "Graph color scheme: Black -- edges between most similar cluster members. Red: Ensure all cluster members are accessible. Blue: Ensure maximum distance across cluster is small enough. Green: Ensure each member is in at least one cycle. Magenta: Intra-cluster connections (made only between clusters of the same charge)" 

#For visualization purposes, weight some edges (only those not added between clusters)
#######need to change
#for edge in Gr.edges():
#    if not edge in newedges:
#        Gr.set_edge_weight( edge, intraweight)

nx.draw_graphviz(Gr)
nx.write_dot(Gr,"final_cluster.dot")

#dot = write(gr)
#DotFile = open('clusters_iterative.dot', 'w')
#DotFile.write(dot)

#Save graph
#import pickle
#file = open('graph.pickle', 'w')
#pickle.dump(gr, file)
#file.close()

########check graph
file04 = open('final_graph.pickle', 'w')
pickle.dump(Gr, file04)
file04.close()

#Save planned calculations in an easy-to-use way. I want to flag calculations within a structural cluster differently than calculations across a structural cluster. Also, store by charge to make that aspect easy.
#Store this information in a couple dictionaries -- first, a cluster_calculations dictionary storing planned calculations by number of cluster; entires are a list with pairs of molecules. Second, an intra_cluster list storing all calculations between clusters. Also want to store indices of clusters by charge


#cluster_calculations = {}
#cross_calculations = []
#for (idx, cluster) in enumerate(clusterlist):
#    cluster_calculations[ idx ] = []
#    for node in cluster:
#        for connectednode in gr[ node ]:
            #Always do pairs in alphabetical order for simplicity and to catch redundancy
#            if node < connectednode:
#                nodepair = [ node, connectednode ]
#            else:
#                nodepair = [connectednode, node]
            #Store pair to appropriate list (either intra-cluster or inter-cluster)
#            if connectednode in cluster: #If intra-cluster:
#                if not nodepair in cluster_calculations[ idx]: cluster_calculations[ idx ].append( nodepair)
#            else:
#                if not nodepair in cross_calculations: cross_calculations.append( nodepair )


#file = open('planned_calculations.pickle', 'w')
#pickle.dump( (cluster_calculations, cross_calculations, clusteridx_by_charge), file)
#file.close()


#Identify all singleton molecules (which can't be connected to any others) and put their structures in a special directory. 
#if not os.path.isdir( 'singletons'): os.mkdir('singletons')
#for cluster in clusterlist:
#    if len(cluster)==1:
####### Sliu 12/20/11 change directory to copy the single moleclue
#        shutil.copy( os.path.join('mol2_file', cluster[0]+'.mol2'), os.path.join( 'singletons', cluster[0]+'.mol2') ) 
