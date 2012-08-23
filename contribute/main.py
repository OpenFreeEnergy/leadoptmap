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
DEBUG = False



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
    basic_rule = rule.Mcs( rule.EqualCharge(), rule.MinimumNumberOfAtom() )

    if (DEBUG) :
        for id in mcs_ids :
            mol_id = mcs.get_parent_ids( id )
            print "%.3f %s" % (basic_rule.similarity( mol_id[0], mol_id[1], mcs_id = id ), KBASE.ask( id )[0].title(),)

    # Gets graph (`g') and clusters (`c').
    g, c = graph.gen_graph( mcs_ids, basic_rule, simi_cutoff = 0.2, max_csize = 128, num_c2c = 1 )

    if (DEBUG) :
        print "%d clusters (counted as the connected components in graph):" % len( c )
        c.sort( lambda x, y : len( x ) - len( y ) )
        for i, e in enumerate( c ) :
            print "cluser #%d, %d structures:" % (i, len( e ),)
            titles = [KBASE.ask( id ).title() for id in e]
            titles.sort()
            for t in titles :
                print "  %s" % t 
 
    l = networkx.spring_layout( g, iterations = 256, weight = "similarity", scale = 10 )
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
            break
        n += 1
    print "--------------------------------------------"
    print "%d files in total" % len( mol_fnames )
    if (len( mol_fnames ) > 1) :
        main( mol_fnames[:128] )
        
