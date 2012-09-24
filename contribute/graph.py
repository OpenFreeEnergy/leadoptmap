"""Functions and classes for creating a graph data structure of molecules
"""


from kbase import KBASE

import struc
import rule
import mcs

import copy
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
            g.add_edge( id0, id1, similarity = simi )

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

    # FIXME: Creates the `desired' from the `complete' instead of from scratch.
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

    n = len( clusters )
    print "%d clusters in total" % n
    for i, c in enumerate( clusters ) :
        print "  size of cluster #%02d: %d" % (i, len( c ),)
    
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
            print "WARNING: Cannot connect cluster #%d with others." % (cluster_index,)
            print "         If there should be connections, consider to adjust the rules to"
            print "         reduce 0-similarity assignments or loosen the MCS conditions."
            continue
        c2c_edges.sort( lambda x, y : cmp_edge( complete, x, y ) )
        connected_clusters = set()
        for k in range( -1, -num_c2c - 1, -1 ) :
            edge  = c2c_edges[k]
            node0 = edge[0]
            node1 = edge[1]
            simi  = complete[node0][node1]["similarity"]
            desired.add_edge( node0, node1, similarity = simi, boundary = True )
            print "boundary similarity = %f between '%s' and '%s'" % (simi, KBASE.ask( node0 ), KBASE.ask( node1 ),)
            for e in unconnected_clusters :
                if (node0 in clusters[e] or node1 in clusters[e]) :
                    connected_clusters.add( e )
        unconnected_clusters -= connected_clusters

    return desired, clusters



def annotate_nodes_with_smarts( g ) :
    """

    """
    for molid in g.nodes() :
        try :
            smarts = KBASE.ask( molid, "SMARTS" )
        except LookupError :
            smarts = KBASE.ask( molid ).smarts()
            KBASE.deposit_extra( molid, "SMARTS", smarts )
        g.node[molid]["SMART"] = smarts
