"""Main program
"""


from kbase import KBASE

import struc
import mcs
import rule
import graph

import glob
import os
import sys
import hashlib
import networkx
import matplotlib.pyplot as plt



# Debug switch
DEBUG = True



def main( mol_fnames ) :
    """
    @type  mol_fnames: C{list} of C{str}'s
    @param mol_fnames: A list of structure file names
    """
    id_list = struc.read_n_files( mol_fnames )
    mols    = []
    for id in id_list :
        mols.append( KBASE.ask( id ) )
        
    mcs_engine = mcs.SchrodMcs( 3 )
    mcs_ids    = mcs_engine.search_all( mols )

    if (DEBUG) :
        for id in mcs_ids :
            mol_id = mcs.get_parent_ids( id )
            print "%.3f %s" % (rule.MCS.similarity( mol_id[0], mol_id[1], mcs_id = id ), KBASE.ask( id )[0].title(),)

    g = graph.create( mcs_ids, rule.MCS )
    l = networkx.spring_layout( g, iterations = 32, weight = "similarity", scale = 100 )
    networkx.draw_networkx( g, pos = l, with_labels = False )
    plt.savefig( "path.png" )

    

if ("__main__" == __name__) :
    argv = sys.argv
    try :
        dir = argv[1]
    except IndexError :
        print "$SCHRODINGER/run %s <%s>" % (argv[0], "structure-file dir",)
        sys.exit( 0 )

    n          = 0
    mol_fnames = glob.glob( dir + "/*.mol2" )
    print "Structure files:"
    for fname in mol_fnames :
        if (n < 16) :
            print os.path.basename( fname )
        elif (n == 16) :
            print "(more)..."
        n += 1
    print "--------------------------------------------"
    print "%d files in total" % len( mol_fnames )
    if (len( mol_fnames ) > 1) :
        main( mol_fnames )
        
