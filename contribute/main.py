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
import pickle
import logging

try :
    import graphviz
except ImportError :
    print "\nWARNING: Graphviz is not installed. Cannot write a .dot output file.\n"



# Logger setup
root = logging.getLogger()
if (root.handlers) :
    for handler in root.handlers :
        root.removeHandler( handler )
        
logging.basicConfig( format  = '%(asctime)s: %(message)s',
                     datefmt = '%m/%d/%y %I:%M:%S',
                     level   = logging.INFO )



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

    logging.info( "MCS searching..." )
    mcs_ids    = mcs_engine.search_all( mols )
    basic_rule = rule.Mcs( rule.EqualCharge(), rule.MinimumNumberOfAtom() )

    logging.info( "MCS searching... Done" )
    for id in mcs_ids :
        mol_id = mcs.get_parent_ids( id )
        logging.debug( "%.3f %s" % (basic_rule.similarity( mol_id[0], mol_id[1], mcs_id = id ), KBASE.ask( id )[0].title(),) )

    # Gets graph (`g') and clusters (`c').
    logging.info( "Creating graph..." )
    g, c = graph.gen_graph( mcs_ids, basic_rule, simi_cutoff = 0.06, max_csize = 64, num_c2c = 1 )
    graph.annotate_nodes_with_smiles ( g )
    graph.annotate_nodes_with_title  ( g )
    graph.annotate_edges_with_smiles ( g )
    graph.annotate_edges_with_hexcode( g )
    logging.info( "Creating graph... Done" )
    
    logging.debug( "%d clusters (counted as the connected components in the graph):" % len( c ) )
    c.sort( lambda x, y : len( x ) - len( y ) )
    for i, e in enumerate( c ) :
        logging.debug( "cluser #%d, %d structures:" % (i, len( e ),) )
        titles = [KBASE.ask( id ).title() for id in e]
        titles.sort()
        for t in titles :
            logging.debug( "  %s" % t )
 
   #l = networkx.spring_layout( g, iterations = 256, weight = "similarity", scale = 10 )
   #networkx.draw_networkx( g, pos = l, with_labels = False )
   #plt.savefig( opt.output + ".png" )

    pkl_fname = opt.output + ".pkl"
    pkl_fh    = open( pkl_fname, "w" )
    pickle.dump( g, pkl_fh )
    pkl_fh.close()
    
    try :
        import graphviz
        
        ag = networkx.to_agraph( g )
        ag.node_attr["fixedsize"] = True
        ag.edge_attr["penwidth" ] = 2.0
        
        simi  = [float( e.attr["similarity"] ) for e in ag.edges()]
        scale = 1.0 / max( simi )
        for e in ag.edges_iter() :
            try :
                partial_ring = int( e.attr["partial_ring"] )
            except (ValueError, TypeError) :
                partial_ring = 0
            saturation       = float( e.attr["similarity"] ) * scale
            saturation       = 0.0 if (saturation < 0) else (1.0 if (saturation > 1) else saturation)
            e.attr["color" ] = "0.8396,%f,0.8" % saturation
            e.attr["weight"] = saturation
            del e.attr["label"]
            if (saturation < 0.01 or partial_ring) :
                e.attr["style"] = "dashed"
        ag.write( opt.output + ".dot" )
    except ImportError :
        logging.warn( "WARNING: Graphviz is not installed. Cannot write a .dot output file." )


    edges = g.edges( data = True )
    logging.info( "%d edges in total" % len( edges ) )
    
    if (opt.siminp) :
        if (opt.siminp_type == "gro") :
            raise NotImplementedError( "Support for writing Gromacs input files is not yet implemented." )
        if (opt.siminp_type == "mae") :
            import schrodinger.application.desmond.fep_mapping as dfm

            tmp_mae_fname = "__temp_file_ok_to_delete_after_running__.mae"
            for id0, id1, attr in edges :
                mol0      = KBASE.ask( id0 )
                mol1      = KBASE.ask( id1 )
                out_fname = "%s_%s_%s.mae" % (opt.siminp, id0[:7], id1[:7],)
                mol0._struc.property["s_fepmap_fragname"] = "none"
                mol1._struc.property["s_fepmap_fragname"] = "%s:%s" % (mol0.title(), mol1.title(),)
                # mol0.write( out_fname, mode = "w" )
                # mol1.write( out_fname, mode = "a" )
                mol0.write( tmp_mae_fname, mode = "w" )
                mol1.write( tmp_mae_fname, mode = "a" )
                try :
                    data = dfm.get_atom_mapping_data( tmp_mae_fname, atomtype = 3 )
                    dfm.write_fepsubst_to_file( data, out_fname )
                except (RuntimeError, NameError,) :
                    logging.warn( "WARNING: Failed to write the input files for '%s' and '%s'." % (mol0, mol1,) )

    tmp_fnames = glob.glob( "__temp_file_ok_to_delete_after_running__.*" )
    for fname in tmp_fnames :
        os.remove( fname )
        


def startup() :
    """
    The startup function, which will handle the command line interface and call the `main' function.
    """
    from optparse import OptionParser

    parser = OptionParser( usage = "Usage: %prog [options] <structure-file-dir | structure-file>...", version = "%prog v0.2" )
    parser.add_option( "-o", "--output", metavar = "BASENAME", default = "simimap",
                       help = "output files' base name. The following files will be written: <basename>.dot, and "
                       "<basename>.pkl." )
    parser.add_option( "-s", "--siminp", metavar = "BASENAME",
                       help = "simulation input files' base name. When this option is specified, a number of input files "
                       "for FEP simulations will be written out." )
    parser.add_option( "-t", "--siminp_type", metavar = "TYPE", default = "mae",
                       help = "simulation input file type [mae | gro]" )

    (opt, args) = parser.parse_args()

    if (len( args ) == 0) :
        parser.print_help()
        sys.exit( 0 )

    molid_list = []
    for a in args :
        logging.info( "Reading structures from '%s'..." % a )
        if (os.path.isfile( a )) :
            try :
                molid_list.extend( struc.read_n_files( [a,] ) )
            except :
                logging.info( "  Cannot read '%s', skip it." % a )
        else :
            logging.info( "  It is a directory, reading *.mol2 and *.mae files in there..." )
            n          = 0
            mol_fnames = glob.glob( a + "/*.mol2" ) + glob.glob( a + "/*.mae" )
            for fname in mol_fnames :
                if (n < 8) :
                    logging.info( "    %s" % os.path.basename( fname ) )
                elif (n == 8) :
                    logging.info( "    (more)..." )
                    break
                n += 1
            logging.info( "    %d files found." % len( mol_fnames ) )
            if (len( mol_fnames ) > 1) :
                molid_list.extend( struc.read_n_files( mol_fnames ) )
                logging.info( "    Reading done." )
    logging.info( "--------------------------------------------" )
    logging.info( "Finish reading structure input files. %d structures in total" % len( molid_list ) )
    if (len( molid_list ) > 1) :
        main( molid_list, opt )


        
if ("__main__" == __name__) :
    startup()
