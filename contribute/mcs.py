"""Maximum common substructure searching (MCSS) engine
"""


from kbase import KBASE

import struc

import os
import subprocess
import hashlib
import logging
import tempfile
import pickle



__junk, __tempfile_rawname = tempfile.mkstemp( dir = "." )
os.remove( __tempfile_rawname )

tempfile_basename = __tempfile_rawname + "__temp_file_ok_to_delete_after_running__"



class Mcs( object ) :
    """
    Base class of MCSS engine
    """
    def __init__( self ) :
        """

        """
        pass

    

    @staticmethod
    def deposit_to_kbase( id0, id1, atom_match0, atom_match1 ) :
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
        mcs_id    = hashlib.sha1( mcs_title ).hexdigest()
        mcs_id    = KBASE.deposit( mcs_id, mcs_title )

        # Sorts the two lists according to the ascending order of atom indices of the first list.
        atom_match0, atom_match1 = zip( *sorted( zip( atom_match0, atom_match1 ), cmp = lambda x, y : x[0] - y[0] ) )
        atom_match0, atom_match1 = list( atom_match0 ), list( atom_match1 )

        KBASE.deposit_extra( mcs_id, "mcs-parents", (id0,             id1,            ) )
        KBASE.deposit_extra( mcs_id, "mcs-matches", {id0:atom_match0, id1:atom_match1,} )
        
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
        def __init__( self, atom_expr = oechem.OEExprOpts_IntType, bond_expr = 0, ringchecking = 'Strict', is_approximate = True ) :
            """
            (To be implemented)
            """
            self._atom_expr      = atom_expr
            self._bond_expr      = bond_expr
            self._ringchecking   = ringchecking
            self._is_approximate = is_approximate
            


        def search( self, mol0, mol1 ) :
            """
            
            """
            mol0 = mol0._struc
            mol1 = mol1._struc

            # Deletes hydrogen atoms before doing MCS search.
            p0 = mol0.CreateCopy()
            p1 = mol1.CreateCopy()
            for mol in (p0, p1,) :
                for atom in mol.GetAtoms() :
                    if (atom.IsHydrogen()) :
                        atom.SetIntType(1)
                    else:
                        atom.SetIntType(2)
            oechem.OESuppressHydrogens(p0)              
            oechem.OESuppressHydrogens(p1)     
            if (self._is_approximate) :
                mcss = oechem.OEMCSSearch( p1, self._atom_expr, self._bond_expr, oechem.OEMCSType_Approximate )
            else:
                mcss = oechem.OEMCSSearch( p1, self._atom_expr, self._bond_expr )

            # Shall we just use 1 as the mininum number of common atoms here?
            # We can later deal with the requirement of different minimum number of common atoms.
            mcss.SetMinAtoms( 1 )
            mcss.SetMCSFunc( oechem.OEMCSMaxAtomsCompleteCycles(1.5) )

            # There could be multiple matches. We select the one with the maximum number of atoms.
            # If there are more than 1 matches with the same maximum number of atoms, we arbitrarily select the first one.
            mcs_mol = None
            max_num = 0
            for match in mcss.Match( p0, True ) :
                # Gets the number of atoms in the MCS.
                num_atom = 0
                mcs_tmp  = oechem.OEMol()
                oechem.OESubsetMol( mcs_tmp, match, True )
                oechem.OEFindRingAtomsAndBonds(mcs_tmp)
                for atom in mcs_tmp.GetAtoms() :
                    if (not atom.IsHydrogen()) :
                        num_atom += 1

                if (num_atom > max_num) :
                    max_num = num_atom
                    mcs_mol = mcs_tmp
                    atom_match0 = []
                    atom_match1 = []
                    for matchpair in match.GetAtoms() :
                        atom_match0.append( matchpair.target .GetIdx()+1 )
                        atom_match1.append( matchpair.pattern.GetIdx()+1 )
                    
            if (mcs_mol) :
                mol0    = struc.OeStruc( mol0 )
                mol1    = struc.OeStruc( mol1 )
                mcs_mol = struc.OeStruc( mcs_mol )
                return self.deposit_to_kbase( mol0.id(), mol1.id(), atom_match0, atom_match1 )
            # Returns `None' if no MCS found.
            
        

        def search_all( self, mols, opt ) :
            """
            N.B.: O(N^2). Needs optimization.
            """
            title_vs_id = {}
            ret     = []
            num_mol = len( mols )
            print "check how many mols in total", num_mol
            for i in range( num_mol ) :
                title = mols[i]._struc.GetTitle()
                if not title_vs_id.has_key(title):
                    title_vs_id[title] = mols[i].id()
                for j in range( i+1, num_mol ) :
                
                    result = self.search( mols[i], mols[j] )
                    if (result) :
                        ret.append( result )
            file_title_v_id = open("title_vs_id.pickle","w")
            pickle.dump(title_vs_id, file_title_v_id)
            file_title_v_id.close()
            return ret

except ImportError :
    pass



try :
    import schrodinger.structutils.analyze as analyze
    
    class McsMatch( object ) :
        """
        A class to temporarily store each entry (i.e., row) in Schrodinger's MCS result
        """
        def __init__( self, mol0_id, mol1_id, mcs_smarts0, mcs_smarts1, mcs_atom0, mcs_atom1 ) :
            """
            @type  mol0_id    : C{str}
            @param mol0_id    : First structure/molecule's ID
            @type  mol1_id    : C{str}
            @param mol1_id    : Second structure/molecule's ID
            @type  mcs_smarts0: C{str}
            @param mcs_smarts0: SMARTS on the first molecule
            @type  mcs_smarts1: C{str}
            @param mcs_smarts1: SMARTS on the second molecule
            @type  mcs_atom0  : C{str}
            @param mcs_atom0  : A list of comma separated integer numbers that are atom indices of the first molecule. If the
                                first molecule is trimmed such that only the specified atoms are left, you get the maximum
                                common substructure.
            @type  mcs_atom1  : C{str}
            @param mcs_atom1  : A list of comma separated integer numbers that are atom indices of the second molecule. If the
                                second molecule is trimmed such that only the specified atoms are left, you get the maximum
                                common substructure.
            """
            self.mol0_id     = mol0_id
            self.mol1_id     = mol1_id
            self.mcs_smarts0 = mcs_smarts0
            self.mcs_smarts1 = mcs_smarts1
            self.mcs_atom0   = mcs_atom0
            self.mcs_atom1   = mcs_atom1

            
            
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



        def search_all( self, mols, opt ) :
            if (not opt.mcs) :
                mae_fname = tempfile_basename + ".mae"
                out_fname = tempfile_basename + ".csv"
                log_fname = tempfile_basename + ".log"
                log_fh    = open( log_fname, "w" )

                if (os.path.isfile( mae_fname )) :
                    os.remove( mae_fname )
                title_vs_id = {}

                for mol in mols :
                    title = mol.title()
                    mol.set_title( mol.id() )
                    mol.write( mae_fname )
                    mol.set_title( title )
                    if not title_vs_id.has_key(title):
                        title_vs_id[title] = mol.id()
                file_title_v_id = open("title_vs_id.pickle","w")
                pickle.dump(title_vs_id, file_title_v_id)
                file_title_v_id.close()
                    
                cmd          = [self._cmd,
                                "-imae",     mae_fname,
                                "-opw",      out_fname,
                                "-atomtype", str( self._typing ),
                                "-nobreakring",
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
            else :
                logging.debug( "DEBUG: Reuse previous MCS searching results: '%s'." % opt.mcs )
                out_fname = opt.mcs
            
            with open( out_fname, "r" ) as fh :
                import csv
                
                lines     = fh.readlines()[1:]
                mcs_match = []
                for tokens in csv.reader( lines ) :
                    mcs_match.append( McsMatch( tokens[1], tokens[3], tokens[11],  tokens[14], tokens[9], tokens[12] ) )

            ret = []
            for m in mcs_match :
                id0  = m.mol0_id
                id1  = m.mol1_id
                mol0 = KBASE.ask( id0 )
                mol1 = KBASE.ask( id1 )

                atom_match0 = [int( i ) for i in m.mcs_atom0.split( ',' )]
                atom_match1 = [int( i ) for i in m.mcs_atom1.split( ',' )]

                ret.append( self.deposit_to_kbase( id0, id1, atom_match0, atom_match1 ) )
                
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



def get_struc( mcs_id ) :
    """

    """
    title                    = KBASE.ask( mcs_id                )
    id0, id1                 = KBASE.ask( mcs_id, "mcs-parents" )
    mcs_matches              = KBASE.ask( mcs_id, "mcs-matches" )
    atom_match0, atom_match1 = mcs_matches[id0], mcs_matches[id1]
    mol0                     = KBASE.ask( id0 )
    mol1                     = KBASE.ask( id1 )
    mcs                      = KBASE.ask( id0 ).extract( atom_match0 )

    for i, e in enumerate( atom_match1, start = 1 ) :
        mcs.atom_prop[i]["mapped_index"] = e

    mcs.set_title( title  )
    mcs.set_id   ( mcs_id )
    return mcs
    


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
