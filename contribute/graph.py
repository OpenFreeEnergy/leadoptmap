"""Functions and classes for creating a graph data structure of molecules
"""


from kbase import KBASE

import struc
import rule
import mcs

import os
import subprocess
import hashlib
import networkx



def create( mcs_ids, rule ) :
    """
    Returns a graph. Node = molecule's ID or name, edge = similarity score
    
    @type  mcs_ids : C{list} of C{str}
    @param mcs_ids : A list of common substructures' IDs
    @type  rule    : C{Rule}
    @param rule    : The rule to determine the similarity score between two structures
    @type  use_name: C{bool}
    @param use_name: If true, node of the graph is molecule's name; otherwise, it is the molecule's ID.
    """
    # Collects all parents' IDs
    ids = set()
    for id in mcs_ids :
        id0, id1 = mcs.get_parent_ids( id )
        ids.add( id0 )
        ids.add( id1 )

    n = len( ids )
    i = 0
    g = networkx.Graph()
    g.add_nodes_from( ids )
        
    for id in mcs_ids :
        id0, id1 = mcs.get_parent_ids( id )
        simi = rule.similarity( id0, id1, mcs_id = id )
        if (simi > 0) :
            g.add_edge( id0, id1, similarity = simi * 0.1 )

    return g



def sign( x ) :
    return 1 if (x > 0) else (-1 if (x < 0) else 0)



def cmp_edge( g, x, y ) :
    return sign( g[x[0]][x[1]]["similarity"] - g[y[0]][y[1]]["similarity"] )



def trim_cluster( g, cluster, num_edges ) :
    """
    Reduces number of edges that each node can have.
    
    @type          g: C{networkx.Graph}
    @param         g: A cluster
    @type  num_edges: C{int}
    @param num_edges: Number of edges that each node is wanted to have
    """
    edges = []
    for e in cluster :
        r = sorted( g.edges( e, data = True ), lambda x, y : cmp_edge( g, x, y ) )
        edges.extend( r[-num_edges:] )
        
    edges     = [(e[0], e[1],) for e in edges]
    del_edges = []
    for e in g.edges( cluster ) :
        if (e not in edges and (e[1], e[0],) not in edges) :
            del_edges.append( e )
    g.remove_edges_from( del_edges )
    return g

    

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
    complete = create( mcs_ids, basic_rule )
    desired  = create( mcs_ids, rule.Cutoff( simi_cutoff, basic_rule ) )
    clusters = sorted( networkx.connected_components( desired ), cmp = lambda x, y : len( x ) - len( y ) )
    largest  = clusters[-1]

    # Does a binary search for an increased cutoff if the largest cluster is too big.
    if (max_csize < len( largest )) :
        high     = 1
        mid      = simi_cutoff
        low      = simi_cutoff
        trial    = 1
        pretrial = 0
        while (abs( trial - pretrial ) > 1E5) :
            n = len( largest )
            print n
            if (max_csize < n) :
                trial = 0.5 * (high + mid)
            elif (int( max_csize * 0.95 ) > n) :
                trial = 0.5 * (mid + low)
            else :
                break
            desired  = create( mcs_ids, rule.Cutoff( trial, basic_rule ) )
            clusters = sorted( networkx.connected_components( desired ), cmp = lambda x, y : len( x ) - len( y ) )
            largest  = clusters[-1]

    # Trims clusters.
    for e in clusters :
        trim_cluster( desired, e, 2 )
        
    # Connects the clusters.
    if (1 < len( clusters )) :
        i = 0
        n = len( clusters )
        while (i < n - 1) :
            j = i + 1
            while (j < n) :
                edges = networkx.edge_boundary( complete, clusters[i], clusters[j] )
                if (len( edges ) == 0) :
                    print "warning: Cannot connect clusters #%d and #%d." % (i, j,)
                    print "         If there should be connections, consider to adjust the rules to"
                    print "         reduce 0-similarity assignments or loosen the MCS conditions."
                    j += 1
                    continue
                edges.sort( lambda x, y : cmp_edge( complete, x, y ) )
                for k in range( -1, -num_c2c - 1, -1 ) :
                    edge  = edges[k]
                    node0 = edge[0]
                    node1 = edge[1]
                    desired.add_edge( node0, node1, similarity = complete[node0][node1]["similarity"], boundary = True )
                j += 1
            i += 1

    return desired, clusters
