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



def main( mol_fnames, output ) :
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
    g, c = graph.gen_graph( mcs_ids, basic_rule, simi_cutoff = 0.1, max_csize = 128, num_c2c = 1 )

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
    plt.savefig( output + ".png" )

    try :
        import graphviz
        
        ag = networkx.to_agraph( g )
        ag.node_attr["fixedsize"] = True
        ag.edge_attr["penwidth" ] = 2.0
        
        for e in ag.nodes_iter() :
            e.attr["label"] = str( e )[:7]

        simi  = [float( e.attr["similarity"] ) for e in ag.edges()]
        scale = 1.0 / max( simi )
        for e in ag.edges_iter() :
            saturation = float( e.attr["similarity"] ) * scale
            saturation = 0.0 if (saturation < 0) else (1.0 if (saturation > 1) else saturation)
            e.attr["color" ] = "0.8396,%f,0.8" % saturation
            e.attr["weight"] = saturation
            if (saturation < 0.01) :
                e.attr["style"] = "dashed"
        ag.write( output + ".dot" )
    except ImportError :
        pass
    
    

if ("__main__" == __name__) :
    from optparse import OptionParser

    parser = OptionParser( usage = "Usage: %prog [options] <structure-file dir>", version = "%prog v0.2" )
    parser.add_option( "-o", "--output", metavar = "BASENAME", default = "simimap",
                       help = "output files' base name (two files will be written: <basename>.png and <basename>.dot)" )

    (opt, args) = parser.parse_args()

    try :
        dir = args[0]
    except IndexError :
        parser.print_help()
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
        main( mol_fnames, output = opt.output )
        
