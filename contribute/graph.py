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
