"""Maximum common substructure searching (MCSS) engine
"""


from kbase import KBASE

import struc

import os
import subprocess
import hashlib



class Mcs( object ) :
    def __init__( self ) :
        """

        """
        pass



    def search( self, mol0, mol1 ) :
        """

        """
        raise NotImplementedError( "`search' method not implemented in subclass" )

    

    def search_all( self, mols ) :
        """

        """
        raise NotImplementedError( "`search' method not implemented in subclass" )
        


try :
    import openeye.oechem as oechem

    class OeMcs( Mcs ) :
        def __init__( self, atom_expr = oechem.OEExprOpts_EqONS,
                      bond_expr = oechem.OEExprOpts_BondOrder|oechem.OEExprOpts_EqSingleDouble|oechem.OEExprOpts_EqAromatic,
                      ringchecking = 'Strict', min_num_atoms = 4, is_approximate = True ) :
            self._atom_expr      = atom_expr
            self._bond_expr      = bond_expr
            self._ringchecking   = ringchecking
            self._min_num_atoms  = min_num_atoms
            self._is_approximate = is_approximate
            


        def search( self, mol0, mol1 ) :
            """
            Returns a `Struc' object, which is maximum common substructure.
            """
            import MCSS_tool
            ret = MCSS_tool.determineMCSS_requireRings( mol1, mol0, debug = False, min_atoms = self._min_num_atoms,
                                                        atomexpr = self._atom_expr, bondexpr = self._bond_expr,
                                                        ringchecking = self._ringchecking, approximate = self._is_approximate )
            # Maybe does some other things.
            
            return ret[0]
        
        

        def search_all( self, mols ) :
            """
            O(N^2). Needs optimization.
            Returns a dictionary: key = pair of molecule indices in the given list, value = `Struc' object.
            """
            ret     = {}
            num_mol = len( mols )
            for i in range( num_mol ) :
                for j in range( i ) :
                    ret[(i, j,)] = self.search( mols[i], mols[j] )
            return ret
        
except ImportError :
    pass



try :
    class McsMatch( object ) :
        def __init__( self, mol_name0, mol_name1, mcs_atom0, mcs_atom1 ) :
            self.mol_name0 = mol_name0.strip()
            self.mol_name1 = mol_name1.strip()
            self.mcs_atom0 = mcs_atom0
            self.mcs_atom1 = mcs_atom1

            
            
    class SchrodMcs( Mcs ) :
        def __init__( self, atomtyping = 3 ) :
            self._cmd    = os.path.join( os.environ['SCHRODINGER'], "utilities", "canvasMCS" )
            self._typing = atomtyping
            


        def search( self, mol0, mol1 ) :
            return self.search_all( [mol0, mol1,] )



        def search_all( self, mols ) :
            tmp_fname = "__temp_file_ok_to_delete_after_running__.mae"
            out_fname = "__temp_file_ok_to_delete_after_running__.csv"
            log_fname = "__temp_file_ok_to_delete_after_running__.log"
            log_fh    = open( log_fname, "w" )

            if (os.path.isfile( tmp_fname )) :
                os.remove( tmp_fname )

            for mol in mols :
                mol.write( tmp_fname )
                
            cmd          = [self._cmd,
                            "-imae",     tmp_fname,
                            "-opw",      out_fname,
                            "-atomtype", str( self._typing ),
                            ]
            mcs_proc     = subprocess.Popen( cmd, stderr = subprocess.STDOUT, stdout = log_fh )
            null, stderr = mcs_proc.communicate()
            val          = mcs_proc.returncode

            if (val == 17) :
                raise RuntimeError( "Used a MCS feature that requires Schrodinger's CANVAS_ELEMENTS license." )
            if (val != 0 ) :
                msg = "CanvasMCS exited prematurely. This could be because the input molecules were too dissimilar" \
                      " or too numerous, or because the chosen atom-typing scheme was too general."
                with open( out_fname ) as fh:
                    msg += "\n\n"
                    msg += fh.read()
                raise RuntimeError( msg )
            
            with open( out_fname, "r" ) as fh :
                import csv
                
                lines     = fh.readlines()[1:]
                mcs_match = []
                for tokens in csv.reader( lines ) :
                    mcs_match.append( McsMatch( tokens[1], tokens[3], tokens[9], tokens[12] ) )

            ret = []
            for e in mcs_match :
                mcs_title = "mcs: %s, %s" % (e.mol_name0, e.mol_name1,)
                if (e.mol_name0 > e.mol_name1) :
                    mcs_title = "mcs: %s, %s" % (e.mol_name1, e.mol_name0,)
                    
                id0       = hashlib.sha1( e.mol_name0 ).hexdigest()
                id1       = hashlib.sha1( e.mol_name1 ).hexdigest()
                mcs_id    = hashlib.sha1( id0 + id1   ).hexdigest()
                atom_list = [int( i ) for i in e.mcs_atom0.split( ',' )]
                mcs_mol0  = KBASE.ask( id0 ).extract( atom_list )
                
                mcs_mol0.set_title( mcs_title )
                KBASE.deposit      ( mcs_id, [mcs_mol0,]                )
                KBASE.deposit_extra( mcs_id, "mcs_parents", (id0, id1,) )
                ret.append( mcs_id )

            return ret
        
except ImportError :
    pass



if ("__main__" == __name__) :
    filenames = ["xfer3.11.mol2", "xfer3.12.mol2",]
    id_list   = struc.read_n_files( filenames )
    mol0      = KBASE.ask( id_list[0] )
    mol1      = KBASE.ask( id_list[1] )
    mcs       = SchrodMcs( 3 )
    mcs_id    = mcs.search( mol0, mol1 )[0]
    mol_id    = KBASE.ask( mcs_id, "mcs_parents" )[0]
    mcs_struc = KBASE.ask( mcs_id )
    mol_struc = KBASE.ask( mol_id )

    out_fname = "out.mae"
    if (os.path.isfile( out_fname )) :
        os.remove( out_fname )
    mol_struc.write( out_fname )
    mcs_struc.write( out_fname )
