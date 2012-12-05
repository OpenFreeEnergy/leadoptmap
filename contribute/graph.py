"""Functions and classes for creating a graph data structure of molecules
"""
 
from kbase import KBASE

import struc
import rule
import mcs
import similarity

import copy
import os
import subprocess
import hashlib
import networkx
import logging



def create( basic_graph, mcs_ids, rule, add_attr = True ) :
    """
    Returns a graph. Node = molecule's ID or name, edge = similarity score
    
    @type  mcs_ids : C{list} of C{str}
    @param mcs_ids : A list of common substructures' IDs
    @type  rule    : C{Rule}
    @param rule    : The rule to determine the similarity score between two structures
    """
    g = copy.deepcopy( basic_graph )
    for id in mcs_ids :
        id0, id1 = mcs.get_parent_ids( id )
        simi     = rule.similarity( id0, id1, mcs_id = id )
        if (simi > 0) :
            if (add_attr) :
                try :
                    partial_ring = int( KBASE.ask( id, "partial_ring" ) )
                except LookupError :
                    partial_ring = 0
                g.add_edge( id0, id1, similarity = simi, partial_ring = partial_ring, mcs_id = id )
            else :
                g.add_edge( id0, id1, similarity = simi )
    return g



def sign( x ) :
    return 1 if (x > 0) else (-1 if (x < 0) else 0)



def cmp_edge( g, x, y ) :
    return sign( g[x[0]][x[1]]["similarity"] - g[y[0]][y[1]]["similarity"] )



def cutoff_graph( g, simi_cutoff ) :
    """
    Generates a new graph by cutting off the similarity scores.

    @type            g: C{networkx.Graph}
    @param           g: Original graph
    @type  simi_cutoff: C{float}
    @param simi_cutoff: Similarity cutoff
    
    @return: A new graph
    """
    g = copy.deepcopy( g )
    edges_to_be_deleted = []
    for e in g.edges() :
        if (g[e[0]][e[1]]["similarity"] < simi_cutoff) :
            edges_to_be_deleted.append( e )
    g.remove_edges_from( edges_to_be_deleted )
    return g

    

def calc_cohesion( g, sg0, sg1, max_csize ) :
    """
    Cohesion score is defined as the sum of similarity scores of all boundary edges between two subgraphs.
    But if the sum of number of nodes of the two subgraphs is greater than C{max_csize}, the score is zero.

    @type          g: C{networkx.Graph}
    @param         g: Original graph
    @type        sg0: C{networkx.Graph}
    @param       sg0: First subgraph
    @type        sg1: C{networkx.Graph}
    @param       sg1: Second subgraph
    @type  max_csize: C{float}
    @param max_csize: Maximum size of a cluster
    
    @return: Cohesion score between the two subgraphs.
    """
    score = 0.0
    n0    = len( sg0 )
    n1    = len( sg1 )
    if (n0 + n1 <= max_csize) :
        boundary_edges = networkx.edge_boundary( g, sg0, sg1 )
        for e in boundary_edges :
            score += g[e[0]][e[1]]["similarity"]
    return score / max( n0, n1 )



def recluster( g, clusters, max_csize ) :
    """
    Regroups a list of clusters.
    
    @type          g: C{networkx.Graph}
    @param         g: The graph
    @type   clusters: C{list} of C{str}
    @param  clusters: A list of clusters to be collapsed
    @type  max_csize: C{int}
    @param max_csize: Maximum cluster size

    @return: A list of clusters.
    """
    while (len( clusters ) > 1) :
        # Step 1: Calculates the cohesion scores for all pairs of clusters.
        clusters = sorted( clusters, cmp = lambda x, y : len( x ) - len( y ) )
        cohesion = []   # Element = (cluster-index pair, score)
        n        = len( clusters )
        for i in range( n - 1 ) :
            for j in range( i + 1, n ) :
                cohesion_score = calc_cohesion( g, clusters[i], clusters[j], max_csize )
                if (cohesion_score > 0) :
                    cohesion.append( (i, j, cohesion_score,) )

        if (not cohesion) :
            break
                
        # Step 2: Finds the smallest cluster with the highest cohesion score. We can do this by sorting `cohesion' twice:
        #         (1) descendingly by scores, and then
        #         (2) ascendingly  by size of the smaller cluster in the pair.
        # Python's sort algorithm is stable.
        cohesion = sorted( cohesion, cmp = lambda x, y : sign( y[2] - x[2] ) )
        cohesion = sorted( cohesion, cmp = lambda x, y :       x[0] - y[0]   )

        # Step 3: Collapses the first pair of clusters in `cohesion'.
        i = cohesion[0][0]
        j = cohesion[0][1]
        collapsed = clusters[i] + clusters[j]
        del clusters[j]
        del clusters[i]
        clusters.append( collapsed )
        
    return clusters


    
def break_cluster( subgraph, orig_cutoff, max_csize ) :
    """
    @type     subgraph: C{networkx.Graph}
    @param    subgraph: The subgraph of the cluster to break
    @type  orig_cutoff: C{float}
    @param orig_cutoff: Original cutoff of similarity scores
    @type    max_csize: C{int}
    @param   max_csize: Maximum cluster size
    """
    # Figures out the number of heavy atoms corresponding to the original cutoff value.
    # This code asssumes the similarity scoring function is `similarity.exp_delta'.
    num_heavy = 64
    while (similarity.exp_delta( num_heavy, 0 ) <= orig_cutoff) :
        num_heavy -= 1

    # Decomposes the original cluster into smaller subclusters.
    cutoff       = similarity.exp_delta( num_heavy, 0 )
    simirule     = rule.Cutoff( cutoff )
    new_subgraph = cutoff_graph( subgraph, cutoff )
    raw_clusters = networkx.connected_components( new_subgraph )
    new_clusters = []
    for c in raw_clusters :
        if (len( c ) > max_csize) :
            new_clusters += break_cluster( subgraph.subgraph( c ), cutoff, max_csize )
        else :
            new_clusters.append( c )

    # Reclusters the subclusters.
    return recluster( subgraph, new_clusters, max_csize )



def trim_cluster( g, cluster, num_edges ) :
    """
    Reduces number of edges that each node can have.
    
    @type          g: C{networkx.Graph}
    @param         g: graph
    @type    cluster: C{list} of {str}
    @param   cluster: A list of nodes of a cluster
    @type  num_edges: C{int}
    @param num_edges: Number of edges that each node is wanted to have
    """
    edges = []
    for e in cluster :
        r = sorted( g.edges( e, data = True ), lambda x, y : cmp_edge( g, x, y ) )
        edges.extend( r[-num_edges:] )

    sg = networkx.Graph( networkx.subgraph( g, cluster ) )
    for e in sg.edges() :
        sg[e[0]][e[1]]["reversed_similarity"] = -sg[e[0]][e[1]]["similarity"]

    mst_edges = networkx.minimum_spanning_edges( sg, "reversed_similarity" )
    
    edges      = [(e[0], e[1],) for e in edges    ]
    edges     += [(e[0], e[1],) for e in mst_edges]
    del_edges  = []
    for e in g.edges( cluster ) :
        if (e not in edges and (e[1], e[0],) not in edges) :
            del_edges.append( e )
    g.remove_edges_from( del_edges )
    return g



def optimize_subgraph( complete, desired ) :
    """
    Optimize a given subgraph to minimize its number of edges. Constraints, such as the "maximum distance" and "circular
    connections", must be satisfied.

    The two arguments: C{complete} and C{desired} are two graphs. Their edges have an attribute called "similarity", whose
    value is the similarity score that we have calculated based on various rules. The score is a floating number within the
    range of [0.0, 1.0]. The greater the score is, the more similar the two compounds are, and in principle the more likely
    the two nodes should be connected.

    C{desired} was created from C{complete} (see below within this docstring), and it has much less number of edges, but they
    still are too many. So the main purpose of this function is to further minimize the C{desired} graph. The optimization
    might not really need C{complete}, but we pass it on just in case we need extra scores that happened to be cut off.
    
    @type  complete: C{networkx.Graph}
    @param complete: A graph where almost any pair of two nodes are connected.
    @type   desired: C{networkx.Graph}
    @param  desired: This graph is created from the C{complete} graph in the following way: First we do cutoff on all
                     similarity scores (meaning that if a score is less than a threshold value it will be set to zero,
                     otherwise it will be kept as is), then we delete edges with zero scores.
    @return: A optimized subgraph
    """
    pass

    

def gen_graph( mcs_ids, basic_rule, simi_cutoff, max_csize, num_c2c ) :
    """
    Generates and returns a graph according to the requirements.
    
    @type      mcs_ids: C{list} of C{str}
    @param     mcs_ids: A list of ids of the maximum substructures in C{KBASE}
    @type         rule: C{rule.Rule}
    @param        rule: The rule to determine the similarity score between two structures
    @type  simi_cutoff: C{float}
    @param simi_cutoff: Cutoff of similarity scores. Values less than the cutoff are considered as 0.
    @type    max_csize: C{int}
    @param   max_csize: Maximum cluster size
    @type      num_c2c: C{int}
    @param     num_c2c: Number of cluster-to-cluster edges
    """
    basic_graph = networkx.Graph()
    all_ids     = set()
    fh          = open( "simiscore", "w" ) if (logging.getLogger().getEffectiveLevel() == logging.DEBUG) else None
    logging.info( "  Calculating similarity scores..." )
    for id in mcs_ids :
        id0, id1 = mcs.get_parent_ids( id )
        simi     = basic_rule.similarity( id0, id1, mcs_id = id )
        KBASE.deposit_extra( id, "similarity", simi )
        all_ids.add( id0 )
        all_ids.add( id1 )
        if (fh) :
            print >> fh, simi
    logging.info( "  Calculating similarity scores... Done" )
        
    basic_graph.add_nodes_from( all_ids )

    complete = create( basic_graph, mcs_ids, rule.Cutoff( 0 ) )
    desired  = cutoff_graph( complete, simi_cutoff )
    clusters = sorted( networkx.connected_components( desired ), cmp = lambda x, y : len( y ) - len( x ) )

    logging.info( "  Original number of clusters: %d" % len( clusters ) )
    num_big_clusters = 0
    for i, c in enumerate( clusters ) :
        logging.info( "    size of cluster #%02d: %d" % (i, len( c ) ),)
        num_big_clusters += (len( c ) > max_csize)

    if (num_big_clusters and False) :
        logging.info( "  %d cluster(s) are too big. Break them into smaller ones. Reclustering..." % num_big_clusters )
        new_clusters = []
        for c in clusters :
            if (max_csize < len( c )) :
                new_clusters += break_cluster( desired.subgraph( c ), simi_cutoff, max_csize )
            else :
                new_clusters.append( c )
        clusters = new_clusters
        logging.info( "  Reclustering... Done" )
    clusters = sorted( clusters, cmp = lambda x, y : len( y ) - len( x ) )

    n = len( clusters )
    logging.info( "  %d clusters in total" % n )
    for i, c in enumerate( clusters ) :
        logging.info( "    size of cluster #%02d: %d" % (i, len( c ),) )

    # Optimizes the subgraphs.
    for e in clusters :
        # FIXME: Replcaes `trim_cluster' with `optimize_subgraph' when the latter is ready.
        trim_cluster( desired, e, 2 )
   
    # Connects the clusters.
    unconnected_clusters = set( range( n ) )
    while (unconnected_clusters) :
        c2c_edges      = []
        cluster_index  = unconnected_clusters.pop()
        this_cluster   = clusters[cluster_index]
        other_clusters = copy.copy( clusters )
        other_clusters.remove( this_cluster )
        for e in other_clusters :
            c2c_edges.extend( networkx.edge_boundary( complete, this_cluster, e ) )
        if (len( c2c_edges ) == 0) :
            logging.warn( "WARNING: Cannot connect cluster #%d with others." % (cluster_index,)      )
            logging.warn( "         If there should be connections, consider to adjust the rules to" )
            logging.warn( "         reduce 0-similarity assignments or loosen the MCS conditions."   )
            continue
        c2c_edges.sort( lambda x, y : cmp_edge( complete, x, y ) )
        connected_clusters = set()
        for k in range( -1, -num_c2c - 1, -1 ) :
            edge   = c2c_edges[k]
            node0  = edge[0]
            node1  = edge[1]
            simi   = complete[node0][node1]["similarity"]
            mcs_id = complete[node0][node1]["mcs_id"    ]
            desired.add_edge( node0, node1, similarity = simi, boundary = True, mcs_id = mcs_id )
            logging.warn( "  boundary similarity = %f between '%s' and '%s'" % (simi, KBASE.ask( node0 ), KBASE.ask( node1 ),) )
            for e in unconnected_clusters :
                if (node0 in clusters[e] or node1 in clusters[e]) :
                    connected_clusters.add( e )
        unconnected_clusters -= connected_clusters

    return desired, clusters



def annotate_nodes_with_smiles( g ) :
    """

    """
    for molid in g.nodes() :
        try :
            smiles = KBASE.ask( molid, "SMILES" )
        except LookupError :
            smiles = KBASE.ask( molid ).smiles()
            KBASE.deposit_extra( molid, "SMILES", smiles )
        g.node[molid]["SMILES"] = smiles



def annotate_nodes_with_title( g ) :
    """

    """
    for molid in g.nodes() :
        g.node[molid]["title"] = KBASE.ask( molid ).title()
        g.node[molid]["label"] = molid[:7]



def annotate_edges_with_smiles( g ) :
    """

    """
    for e in g.edges( data = True ) :
        try :
            mcs_id = e[2]["mcs_id"]
            try :
                smiles = KBASE.ask( mcs_id, "SMILES" )
            except LookupError :
                smiles = mcs.get_struc( mcs_id ).smiles()
            KBASE.deposit_extra( mcs_id, "SMILES", smiles )
            g[e[0]][e[1]]["SMILES"] = smiles
        except KeyError :
            pass



def annotate_edges_with_matches( g ) :
    """

    """
    for e in g.edges( data = True ) :
        try :
            mcs_id      = e[2]["mcs_id"]
            mol0        = KBASE.ask( e[0] )
            mol1        = KBASE.ask( e[1] )
            mcs_matches = KBASE.ask( mcs_id, "mcs-matches" )
            trimmed_mcs = KBASE.ask( mcs_id, "trimmed-mcs" )
            g[e[0]][e[1]]["original-mcs"] = {e[0]:mol0.smarts( mcs_matches[e[0]] ), e[1]:mol1.smarts( mcs_matches[e[1]] ),}
            g[e[0]][e[1]][ "trimmed-mcs"] = trimmed_mcs
        except KeyError :
            pass



def annotate_edges_with_hexcode( g ) :
    """

    """
    for e in g.edges() :
        g[e[0]][e[1]]["label"] = "%s-%s" % (e[0][:7], e[1][:7],)



if ("__main__" == __name__) :
    import pickle

    fh0 = open( "sampl3_complete.pkl" )
    fh1 = open( "sampl3_desired.pkl"  )
    complete = pickle.load( fh0 )
    desired  = pickle.load( fh1 )
    optimize_subgraph( complete, desired )
    break_cluster( desired, 0.4, 64 )
