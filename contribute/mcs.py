"""Maximum common substructure searching (MCSS) engine
"""


from kbase import KBASE

import struc

import os
import subprocess
import hashlib



class Mcs( object ) :
    """
    Base class of MCSS engine
    """
    def __init__( self ) :
        """

        """
        pass

    

    @staticmethod
    def deposit_to_kbase( id0, id1, atom_match0, atom_match1, mcs_mol ) :
        """
        Deposits a MCS substructure and relevant information into the kbase and returns its ID in the C{KBASE}.

        @type         id0: C{str}
        @param        id0: ID of the first (reference) molecule in the C{KBASE}
        @type         id1: C{str}
        @param        id1: ID of the second molecule in the C{KBASE}
        @type  atom_match: C{list} of C{int}
        @param atom_match: A list of atom indices of matches atoms in the reference molecule
        @type     mcs_mol: C{Struc}
        @param    mcs_mol: C{Struc} object of the MCS substructure        
        """
        mol0      = KBASE.ask( id0 )
        mol1      = KBASE.ask( id1 )
        name0     = mol0.title()
        name1     = mol1.title()
        mcs_title = "mcs@%s..%s" % (name0, name1,)
        mcs_mol.set_title( mcs_title )
        mcs_mol.set_id   ( None      )
        
        mcs_id = KBASE.deposit( mcs_mol.id(), [mcs_mol,] )
        if (mcs_mol.id() != mcs_id) :
            mcs_mol.set_id( mcs_id )
            KBASE.deposit( mcs_id, [mcs_mol,], should_overwrite = True )
        
        KBASE.deposit_extra( mcs_id, "mcs-parents", (id0, id1,)                 )
        KBASE.deposit_extra( mcs_id, "mcs-matches", (atom_match0, atom_match1,) )
        return mcs_id
        
    

    def search( self, mol0, mol1 ) :
        """
        Finds out the maximum common substructure between C{mol0} and C{mol1} and deposits it in C{KBASE}. Returns the ID of
        the substructure in C{KBASE}.

        @type  mol0: C{Struc}
        @param mol0: First molecule (the reference molecule)
        @type  mol1: C{Struc}
        @param mol1: Second molecule
        """
        raise NotImplementedError( "`search' method not implemented in subclass" )

    

    def search_all( self, mols ) :
        """
        Finds out the maximum common substructures between any pair of the given structures and deposits them in C{KBASE}.
        Returns a list of IDs of the substructures in C{KBASE}.

        @type  mols: C{list} of C{Struc}
        @param mols: A list of molecules
        """
        raise NotImplementedError( "`search' method not implemented in subclass" )
        


try :
    import openeye.oechem as oechem

    class OeMcs( Mcs ) :
        def __init__( self, atom_expr = oechem.OEExprOpts_EqONS,
                      bond_expr = oechem.OEExprOpts_BondOrder|oechem.OEExprOpts_EqSingleDouble|oechem.OEExprOpts_EqAromatic,
                      ringchecking = 'Strict', min_num_atoms = 4, is_approximate = True ) :
            """
            (To be implemented)
            """
            self._atom_expr      = atom_expr
            self._bond_expr      = bond_expr
            self._ringchecking   = ringchecking
            self._min_num_atoms  = min_num_atoms
            self._is_approximate = is_approximate
            


        def search( self, mol0, mol1 ) :
            """
            (To be implemented)
            """
            import MCSS_tool
            ret = MCSS_tool.determineMCSS_requireRings( mol1, mol0, debug = False, min_atoms = self._min_num_atoms,
                                                        atomexpr = self._atom_expr, bondexpr = self._bond_expr,
                                                        ringchecking = self._ringchecking, approximate = self._is_approximate )
            # Maybe does some other things.
            
            return ret[0]
        
        

        def search_all( self, mols ) :
            """
            N.B.: O(N^2). Needs optimization.
            (To be implemented)
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
        """
        A class to temporarily store each entry (i.e., row) in Schrodinger's MCS result
        """
        def __init__( self, mol0_id, mol1_id, mcs_atom0, mcs_atom1 ) :
            """
            @type  mol0_id  : C{str}
            @param mol0_id  : First structure/molecule's ID
            @type  mol1_id  : C{str}
            @param mol1_id  : Second structure/molecule's ID
            @type  mcs_atom0: C{str}
            @param mcs_atom0: A list of comma separated integer numbers that are atom indices of the first molecule. If the
                              first molecule is trimmed such that only the specified atoms are left, you get the maximum
                              common substructure.
            @type  mol_name1: C{str}
            @param mcs_atom1: A list of comma separated integer numbers that are atom indices of the second molecule. If the
                              second molecule is trimmed such that only the specified atoms are left, you get the maximum
                              common substructure.
            """
            self.mol0_id   = mol0_id
            self.mol1_id   = mol1_id
            self.mcs_atom0 = mcs_atom0
            self.mcs_atom1 = mcs_atom1

            
            
    class SchrodMcs( Mcs ) :
        def __init__( self, atomtyping = 3 ) :
            """
            @type  atomtyping: C{int}
            @param atomtyping: Schrodinger Canvas' atom typing scheme (see below)
                                1 - All atoms equivalent; all bonds equivalent.
                                2 - Atoms distinguished by HB acceptor/donor; all bonds equivalent.
                                3 - Atoms distinguished by hybridization state; all bonds equivalent.
                                4 - Atoms distinguished by functional type: {H}, {C}, {F,Cl}, {Br,I}, {N,0},
                                    {S}, {other}; bonds by hybridization.
                                5 - Mol2 atom types; all bonds equivalent.
                                6 - Atoms distinguished by whether terminal, halogen, HB acceptor/donor;
                                    bonds distinguished by bond order.
                                7 - Atomic number and bond order.
                                8 - Atoms distinguished by ring size, aromaticity, HB acceptor/donor,
                                    ionization potential, whether terminal, whether halogen; bonds
                                    distinguished by bond order.
                                9 - Carhart atom types (atom-pairs approach); all bonds equivalent.
                               10 - Daylight invariant atom types; bonds distinguished by bond order.
                               11 - Same as 7, but distinguishing aromatic from non-aromatic.
                               12 - Same as 11, but distinguishing aliphatic atoms by ring/acyclic.
                               13 - Same as 12, but distinguishing rings by size.
                                C - Custom. Must be followed by location of a type definitions file.
            """
            self._cmd    = os.path.join( os.environ['SCHRODINGER'], "utilities", "canvasMCS" )
            self._typing = atomtyping
            


        def search( self, mol0, mol1 ) :
            return self.search_all( [mol0, mol1,] )



        def search_all( self, mols ) :
            mae_fname = "__temp_file_ok_to_delete_after_running__.mae"
            out_fname = "__temp_file_ok_to_delete_after_running__.csv"
            log_fname = "__temp_file_ok_to_delete_after_running__.log"
            log_fh    = open( log_fname, "w" )

            if (os.path.isfile( mae_fname )) :
                os.remove( mae_fname )

            for mol in mols :
                title = mol.title()
                mol.set_title( mol.id() )
                mol.write( mae_fname )
                mol.set_title( title )
                
            cmd          = [self._cmd,
                            "-imae",     mae_fname,
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
                id0  = e.mol0_id
                id1  = e.mol1_id
                mol0 = KBASE.ask( id0 )
                mol1 = KBASE.ask( id1 )

                atom_match0 = [int( i ) for i in e.mcs_atom0.split( ',' )]
                atom_match1 = [int( i ) for i in e.mcs_atom1.split( ',' )]
                mcs_mol0    = mol0.extract( atom_match0 )
                ret.append( self.deposit_to_kbase( id0, id1, atom_match0, atom_match1, mcs_mol0 ) )
                
            return ret
        
except ImportError :
    pass        



def get_parent_ids( mcs_id ) :
    """
    Returns a pair of IDs of the common substructure's parents.

    @type  mcs_id: C{str}
    @param mcs_id: ID of the common substructure
    """
    return KBASE.ask( mcs_id, "mcs-parents" )



if ("__main__" == __name__) :
    filenames = ["xfer3.11.mol2", "xfer3.12.mol2",]
    id_list   = struc.read_n_files( filenames )
    mol0      = KBASE.ask( id_list[0] )
    mol1      = KBASE.ask( id_list[1] )
    mcs       = SchrodMcs( 3 )
    mcs_id    = mcs.search( mol0, mol1 )[0]
    mol_id    = KBASE.ask( mcs_id, "mcs-parents" )[0]
    mcs_struc = KBASE.ask( mcs_id )[0]
    mol_struc = KBASE.ask( mol_id )

    out_fname = "out.mae"
    if (os.path.isfile( out_fname )) :
        os.remove( out_fname )
    mol_struc.write( out_fname )
    mcs_struc.write( out_fname )
