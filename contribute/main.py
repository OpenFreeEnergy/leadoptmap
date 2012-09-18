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

try :
    import graphviz
except ImportError :
    print "WARNING: Graphviz is not installed. Cannot write a .dot output file."



# Debug switch
DEBUG = False



def main( molid_list, opt ) :
    """
    @type  molid_list: C{list} of C{str}'s
    @param molid_list: A list of molecule IDs in the C{KBASE}
    """
    mols = []
    for id in molid_list :
        mols.append( KBASE.ask( id ) )

    if   (struc.infrastructure == "schrodinger") : mcs_engine = mcs.SchrodMcs( 3 )
    elif (struc.infrastructure == "oechem"     ) : mcs_engine = mcs.OeMcs()
    
    mcs_ids    = mcs_engine.search_all( mols )
    basic_rule = rule.Mcs( rule.EqualCharge(), rule.MinimumNumberOfAtom() )

    if (DEBUG) :
        for id in mcs_ids :
            mol_id = mcs.get_parent_ids( id )
            print "%.3f %s" % (basic_rule.similarity( mol_id[0], mol_id[1], mcs_id = id ), KBASE.ask( id )[0].title(),)

    # Gets graph (`g') and clusters (`c').
    g, c = graph.gen_graph( mcs_ids, basic_rule, simi_cutoff = 0.1, max_csize = 64, num_c2c = 1 )

    if (DEBUG) :
        print "%d clusters (counted as the connected components in the graph):" % len( c )
        c.sort( lambda x, y : len( x ) - len( y ) )
        for i, e in enumerate( c ) :
            print "cluser #%d, %d structures:" % (i, len( e ),)
            titles = [KBASE.ask( id ).title() for id in e]
            titles.sort()
            for t in titles :
                print "  %s" % t 
 
    l = networkx.spring_layout( g, iterations = 256, weight = "similarity", scale = 10 )
    networkx.draw_networkx( g, pos = l, with_labels = False )
    plt.savefig( opt.output + ".png" )

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
        ag.write( opt.output + ".dot" )
    except ImportError :
        print "WARNING: Graphviz is not installed. Cannot write a .dot output file."


    edges = g.edges( data = True )
    print "%d edges in total" % len( edges )
    
    if (opt.siminp) :
        if (opt.siminp_type == "gro") :
            raise NotImplementedError( "Support for writing Gromacs input files is not yet implemented." )
        if (opt.siminp_type == "mae") :
            import schrodinger.application.desmond.fep_mapping as dfm

            tmp_mae_fname = "__temp_file_ok_to_delete_after_running__.mae"
            for id0, id1, attr in edges :
                print id0[:7], id1[:7], attr['similarity']
                mol0 = KBASE.ask( id0 )
                mol1 = KBASE.ask( id1 )
                print mol0, mol1
                mol0._struc.property["s_fep_fragname"] = "none"
                mol1._struc.property["s_fep_fragname"] = "%s:%s" % (mol0.title(), mol1.title(),)
                mol0.write( tmp_mae_fname, mode = "w" )
                mol1.write( tmp_mae_fname, mode = "a" )
                try :
                    data = dfm.get_atom_mapping_data( tmp_mae_fname, atomtype = 3 )
                    dfm.write_fepsubst_to_file( data, "%s_%s_%s.mae" % (opt.siminp, id0[:7], id1[:7],) )
                except (RuntimeError, NameError,) :
                    print "WARNING: Failed to write the input files for '%s' and '%s'." % (mol0, mol1,)



if ("__main__" == __name__) :
    from optparse import OptionParser

    parser = OptionParser( usage = "Usage: %prog [options] <structure-file-dir | structure-file>...", version = "%prog v0.2" )
    parser.add_option( "-o", "--output", metavar = "BASENAME", default = "simimap",
                       help = "output files' base name (two files will be written: <basename>.png and <basename>.dot)" )
    parser.add_option( "-s", "--siminp", metavar = "BASENAME",
                       help = "simulation input files' base name" )
    parser.add_option( "-t", "--siminp_type", metavar = "TYPE", default = "mae",
                       help = "simulation input file type [mae | gro]" )

    (opt, args) = parser.parse_args()

    if (len( args ) == 0) :
        parser.print_help()
        sys.exit( 0 )

    molid_list = []
    for a in args :
        print "Reading structures from '%s'..." % a
        if (os.path.isfile( a )) :
            try :
                molid_list.extend( struc.read_n_files( [a,] ) )
            except :
                print "  Cannot read '%s', skip it." % a
        else :
            print "  It is a directory, reading *.mol2 and *.mae files in there..."
            n          = 0
            mol_fnames = glob.glob( a + "/*.mol2" ) + glob.glob( a + "/*.mae" )
            for fname in mol_fnames :
                if (n < 8) :
                    print "    %s" % os.path.basename( fname )
                elif (n == 8) :
                    print "    (more)..."
                    break
                n += 1
            print "    %d files found." % len( mol_fnames )
            if (len( mol_fnames ) > 1) :
                print "    Reading them..."
                molid_list.extend( struc.read_n_files( mol_fnames ) )
                print "      Done."
    print "--------------------------------------------"
    print "Finish reading structure input files. %d structures in total" % len( molid_list )
    if (len( molid_list ) > 1) :
        main( molid_list, opt )
        
