"""We implement a minimalist rule engine, and a few basic rule classes.
"""



from kbase import KBASE

import mcs
import similarity

import hashlib
import logging



class Rule( object ) :
    """
    Base class of all rule classes.
    """
    def __init__( self, *subrules ) :
        self._subrules = subrules



    def _similarity( self, id0, id1, **kwarg ) :
        """
        Given the IDs of two molecular structures in the C{KBASE}, return a similarity score of the two molecules.
        By default, we return 1, assuming there's no difference at the most basic level.
        Subclass should normally override this method.
        
        @type  id0: C{str}
        @param id0: ID of the first molecule in the C{KBASE}
        @type  id1: C{str}
        @param id1: ID of the second molecule in the C{KBASE}
        """
        return 1.0


    
    def similarity( self, id0, id1, **kwarg ) :
        """
        Given the IDs of two molecular structures in the C{KBASE}, return a similarity score of the two molecules with all
        subrules combined.
        Similarity score is a floating number in the range of [0, 1].

        Each subrule will return a similarity score, and all the scores will be multiplied together to the score returned by
        the method C{_similarity}. And the product will be returned as the final result of this rule.

        @type  id0: C{str}
        @param id0: ID of the first molecule in the C{KBASE}
        @type  id1: C{str}
        @param id1: ID of the second molecule in the C{KBASE}
        """
        #print "inside rule.py main rule check id", id0,id1
        result = self._similarity( id0, id1, **kwarg )
        if (result > 0) :
            for e in self._subrules :
                if (isinstance( e, list )) :
                    result *= max( [e.similarity( id0, id1, **kwarg )] )
                else :
                    result *= e.similarity( id0, id1, **kwarg )
        return result



class MinimumNumberOfAtom( Rule ) :
    """
    Rule on minimum number of heavy atoms in the maximum common substructure

    Similarity score is 0 if the minimum number of heavy atoms in the maximum common substructure is less than a specified
    threshold value, or 1 if otherwise.
    """
    def __init__( self, threshold = 4, *subrules ) :
        """
        @type   threshold: C{int}
        @param  threshold: Threshold of number of atoms. Above this threshold, the similarity score is 1; otherwise it's 0.
        @type  heavy_only: C{bool}
        @param heavy_only: If the value is C{True}, we count only the heavy atoms; if it's C{False}, we count all atoms.
        """
        Rule.__init__( self, *subrules )
        self._threshold  = threshold


        
    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        num_atom_mcs  = KBASE.ask( kwarg["mcs_id"], "num_heavy_atoms" )
        num_atom_mol0 = len( KBASE.ask( id0 ).heavy_atoms() )
        num_atom_mol1 = len( KBASE.ask( id1 ).heavy_atoms() )

        return float( (num_atom_mcs  >= self._threshold    ) or
                      (num_atom_mol0 <  self._threshold + 3) or
                      (num_atom_mol1 <  self._threshold + 3) )



class Cutoff( Rule ) :
    """
    Rule that we ``cut off'' a similarity score. Cutting off here means that we consider a score to be zero if it is less than
    a threshold value.
    """
    def __init__( self, cutoff, *subrules ) :
        """
        @type  cutoff: C{float}
        @param cutoff: Cutoff threshold. Above this threshold, the similarity score is 1; otherwise it's 0.
        """
        Rule.__init__( self, *subrules )
        self._cutoff  = cutoff

        
        
    def similarity( self, id0, id1, **kwarg ) :
        try :
            mcs_id = kwarg["mcs_id"]
            simi   = KBASE.ask( mcs_id, "similarity" )
        except KeyError :
            simi = Rule.similarity( self, id0, id1, **kwarg )
        if (simi < self._cutoff) :
            simi = 0.0
        return simi



class EqualCharge( Rule ) :
    """
    The two molecules must be of the same net charge; otherwise, the similarity score is zero.
    """
    def __init__( self, *subrules ) :
        """
        """
        Rule.__init__( self, *subrules )

        
        
    def _similarity( self, id0, id1, **kwarg ) :
        # Retrieves the parent molecules of the MCS.
        mol0 = KBASE.ask( id0 )
        mol1 = KBASE.ask( id1 )
        if (mol0.total_charge() != mol1.total_charge()) :
            return 0.0
        return 1.0



class Mcs( Rule ) :
    """
    MCS-based rule
    Similarity is scored using the C{similarity.by_heavy_atom_count} (see the C{similarity} module).
    """
    def __init__( self, *subrules ) :
        Rule.__init__( self, *subrules )
        
        

    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        mcs_id = kwarg["mcs_id"]
        mcs0   = mcs.get_struc( mcs_id )
        mol0   = KBASE.ask( id0 )
        mol1   = KBASE.ask( id1 )

        num_heavy_atoms = len( mcs0.heavy_atoms() )
        num_light_atoms = len( mcs0.atom ) - num_heavy_atoms
        
        KBASE.deposit_extra( mcs_id, "num_heavy_atoms", num_heavy_atoms )
        KBASE.deposit_extra( mcs_id, "num_light_atoms", num_light_atoms )
        
        return similarity.by_heavy_atom_count( mol0, mol1, mcs0 )



class TrimMcs( Rule ) :
    """
    Delete chiral atoms and partial ring atoms from MCS, return a score.
    """
    def __init__( self, strict = True, *subrules ) :
        Rule.__init__( self, *subrules )
        self._strict = strict
        


    def _delete_broken_ring( self, mol0, mol1, mcs0 ) :
        mcs_ring_atoms = mcs0.ring_atoms()
        mcs_nonr_atoms = set( range( 1, len( mcs0.atom ) + 1 ) ) - mcs_ring_atoms
        mo0_ring_atoms = mol0.ring_atoms()
        mo1_ring_atoms = mol1.ring_atoms()

        mo0_conflict   = set( [mcs0.atom_prop[i][  "orig_index"] for i in mcs_nonr_atoms] ) & mo0_ring_atoms
        mo1_conflict   = set( [mcs0.atom_prop[i]["mapped_index"] for i in mcs_nonr_atoms] ) & mo1_ring_atoms

        def extend_conflict_to_whole_ring( mol, conflict ) :
            ring_agrps = mol.ring_atoms( aromaticity = 0, group = True )
            old_len    = 0
            while (old_len != len( conflict )) :
                old_len     = len( conflict )
                atom_to_add = set()
                for a in conflict :
                    for i, r in enumerate( ring_agrps ) :
                        if (a in r) :
                            atom_to_add   |= set( r )
                            ring_agrps[i]  = []
                conflict |= atom_to_add
                
        def extend_conflict_to_nonaromatic_ring( mol, conflict ) :
            # We have to be careful here. When we delete non-aromatic rings, we should avoid deleting atoms involved in both
            # aromatic and no-aromatic rings.
            ring_agrps = mol.ring_atoms( aromaticity = -1, group = True )
            arom_rings = mol.ring_atoms( aromaticity =  1, group = True )
            for r in ring_agrps :
                for ar in arom_rings :
                    r -= ar
            old_len = 0
            while (old_len != len( conflict )) :
                old_len     = len( conflict )
                atom_to_add = set()
                for a in conflict :
                    for i, r in enumerate( ring_agrps ) :
                        if (a in r) :
                            atom_to_add   |= set( r )
                            ring_agrps[i]  = []
                conflict |= atom_to_add

        if (self._strict) :
            extend_conflict_to_whole_ring( mol0, mo0_conflict )
            extend_conflict_to_whole_ring( mol1, mo1_conflict )
        else :
            extend_conflict_to_nonaromatic_ring( mol0, mo0_conflict )
            extend_conflict_to_nonaromatic_ring( mol1, mo1_conflict )
        
        # Now indices in mo0_conflict and mo1_conflict are indices in mo0 and mo1, respectively.
        # We need to map them back to the indices in mcs0.
        mo0_to_mcs = {}
        mo1_to_mcs = {}
        for i in range( 1, len( mcs0.atom ) + 1 ) :
            mo0_to_mcs[mcs0.atom_prop[i][  "orig_index"]] = i
            mo1_to_mcs[mcs0.atom_prop[i]["mapped_index"]] = i
        mo0_conflict = set( [mo0_to_mcs[i] for i in mo0_conflict if (i in mo0_to_mcs)] )
        mo1_conflict = set( [mo1_to_mcs[i] for i in mo1_conflict if (i in mo1_to_mcs)] )

        conflict = list( mo0_conflict | mo1_conflict )
        mcs0.delete_atom( conflict )
        return conflict
        
        

    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        mcs_id = kwarg["mcs_id"]
        mcs0   = mcs.get_struc( mcs_id ).copy()
        mol0   = KBASE.ask( id0 )
        mol1   = KBASE.ask( id1 )

        orig_num_heavy_atoms = len( mcs0.heavy_atoms() )
        partial_ring         = self._delete_broken_ring( mol0, mol1, mcs0 )
            
        # Deletes chiral atoms.
        chiral_atoms = mcs0.chiral_atoms()
        ring_atoms   = mcs0.  ring_atoms()
        chiral_atoms.sort( reverse = True )
        for atom_index in chiral_atoms :
            if (atom_index in ring_atoms) :
                bonded_atoms = set( mcs0.bonded_atoms( atom_index ) ) - ring_atoms
                if (bonded_atoms) :
                    i = 0
                    n = 0
                    for atom in bonded_atoms :
                        cp = mcs0.copy()
                        cp.delete_atom( atom )
                        m = len( cp.atom )
                        if (m > n) :
                            i = atom
                            n = m
                    mcs0.delete_atom( i )
                else :
                    logging.warn( "WARNING: Cannot delete chiral atom #%d in structure: %s" % (atom_index, mcs0.title(),) )
            else :
                # If the chiral atom is not a ring atom, we simply delete it.
                mcs0.delete_atom( atom_index )

        # If the deletion results in multiple unconnected fragments, we keep only the biggest one.
        atoms_to_delete = []
        for e in mcs0.molecules()[1:] :
            atoms_to_delete.extend( e )
        mcs0.delete_atom( atoms_to_delete )

        # Gets the SMARTS for the trimmed structure.
        atom_list0 = []
        atom_list1 = []
        for e in mcs0.atom_prop[1:] :
            atom_list0.append( e[  "orig_index"] )
            atom_list1.append( e["mapped_index"] )

        smarts0 = mol0.smarts( atom_list0 )
        try :
            smarts1 = mol1.smarts( atom_list1 )
        except ValueError :
            #print atom_list0, atom_list1
            #print mol0.title(), mol1.title()
            smarts1 = ""
        
        KBASE.deposit_extra( mcs_id, "trimmed-mcs",  {id0:smarts0, id1:smarts1,} )
        KBASE.deposit_extra( mcs_id, "partial_ring", len( partial_ring ) )

        num_heavy_atoms = len( mcs0.heavy_atoms() )
        num_light_atoms = len( mcs0.atom ) - num_heavy_atoms
        
        KBASE.deposit_extra( mcs_id, "num_heavy_atoms", num_heavy_atoms )
        KBASE.deposit_extra( mcs_id, "num_light_atoms", num_light_atoms )

        return similarity.exp_delta( 2 * (orig_num_heavy_atoms - num_heavy_atoms), 0 )



class TrimMcs_oe( Rule ) :
    """
    Delete chiral atoms and partial ring atoms from MCS, return a score.
    """
    def __init__( self, strict_ring_checking = True, *subrules ) :
        Rule.__init__( self, *subrules )
        self._strict_ring_checking = strict_ring_checking



    def _delete_broken_ring( self, mol0, mol1, mcs0, strict_ring_checking = True ) :
        mcs_ring_atoms = mcs0.ring_atoms()
        mcs_aromatic_atoms = mcs0.aromatic_atoms()
        #print "Check ring size inside rule", mcs0.ring_size()
        mcs0_ring_dic = mcs0.ring_size()
        mcs_nonr_atoms = set( range( 1, len( mcs0.atom ) + 1 ) ) - mcs_ring_atoms
        mo0_ring_atoms = mol0.ring_atoms()
        mo0_aromatic_atoms = mol0.aromatic_atoms()
        #print "Check ring size inside rule mol0", mol0.ring_size()
        mol0_ring_dic = mol0.ring_size()
        mo1_ring_atoms = mol1.ring_atoms()
        mo1_aromatic_atoms = mol1.aromatic_atoms()
        #print "Check ring size inside rule mol1", mol1.ring_size()
        mol1_ring_dic = mol1.ring_size()
        mo0_conflict = []
        mo1_conflict = []
        #Strict ring check 
        #Delete atoms which change ring size either from mol0 to mcs or from mol1 to mcs
        if strict_ring_checking:

            for i in mcs0_ring_dic.keys():
                mol0_key = mcs0.atom_prop[i][  "orig_index"]
                mol1_key = mcs0.atom_prop[i][  "mapped_index"]

                if mcs0_ring_dic[i] <> mol0_ring_dic[mol0_key]:
                    #print "The ring size is different mcs0 : %s, mol0 :%s i : %s \n" %(mcs0_ring_dic[i], mol0_ring_dic[mol0_key], i)
                    mo0_conflict.append(i)
                elif mcs0_ring_dic[i] <> mol1_ring_dic[mol1_key]:
                    #print "The ring size is different mcs0 : %s, mol1 :%s i : %s \n" %(mcs0_ring_dic[i], mol1_ring_dic[mol1_key], i)
                    mo1_conflict.append(i)
            #print "Check confilict idx" , mo0_conflict, mo1_conflict

        #Unstrict ring check 
        #Delete all atoms which either (a) in a ring in mol0 but not in a ring in mcs, or (b) in a ring in mol1 but not in a ring in mcs, or (c) in a ring in mcs and not in aromatic ring in either mcs, mol0, mol1 if the ring size is changed from mol0 to mcs or mol1 to mcs.  
        else:
            for i in mcs0_ring_dic.keys():                    
                mol0_key = mcs0.atom_prop[i][  "orig_index"]  
                mol1_key = mcs0.atom_prop[i][  "mapped_index"]
                                                              
                if mcs0_ring_dic[i] == 0 and mol0_ring_dic[mol0_key] > 0:   
                    mo0_conflict.append(i)                    
                elif mcs0_ring_dic[i] == 0 and mol1_ring_dic[mol1_key] >0:  
                    mo1_conflict.append(i)                    
                elif mcs0_ring_dic[i] > 0 and (not i in mcs_aromatic_atoms or not mol0_key in mo0_aromatic_atoms or not mol1_key in mo1_aromatic_atoms):
                    if mcs0_ring_dic[i] <> mol0_ring_dic[mol0_key]:         
                        mo0_conflict.append(i)                
                    elif mcs0_ring_dic[i] <> mol1_ring_dic[mol1_key]:       
                        mo1_conflict.append(i)

        conflict = list( set(mo0_conflict) | set(mo1_conflict) )
        mcs0.delete_atom( conflict )
        return conflict
        
    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        mcs_id = kwarg["mcs_id"]
        mcs0   = mcs.get_struc( mcs_id ).copy()
        mol0   = KBASE.ask( id0 )
        mol1   = KBASE.ask( id1 )

        orig_num_heavy_atoms = len( mcs0.heavy_atoms() )
        # Deletes chiral atoms.
        chiral_atoms = mcs0.chiral_atoms()
        ring_atoms   = mcs0.  ring_atoms()
        chiral_atoms.sort( reverse = True )
        for atom_index in chiral_atoms :
            if (atom_index in ring_atoms) :
                bonded_atoms = set( mcs0.bonded_atoms( atom_index ) ) - ring_atoms
                if (bonded_atoms) :
                    i = 0
                    n = 0
                    for atom in bonded_atoms :
                        cp = mcs0.copy()
                        cp.delete_atom( atom )
                        m = len( cp.atom )
                        if (m > n) :
                            i = atom
                            n = m
                    mcs0.delete_atom( i )
                    mcs0 = mcs0.copy()
                else :
                    logging.warn( "WARNING: Cannot delete chiral atom #%d in structure: %s" % (atom_index, mcs0.title(),) )
            else :
                # If the chiral atom is not a ring atom, we simply delete it.
                mcs0.delete_atom( atom_index )
                mcs0 = mcs0.copy()

        # If the deletion results in multiple unconnected fragments, we keep only the biggest one.
        #print "After the deletion of chiral atoms"
        #raw_input()
        mcs0 = mcs0.copy()
        atoms_to_delete = []
        for e in mcs0.molecules()[1:] :
            atoms_to_delete.extend( e )
        #print "Rule Atom to delete", atoms_to_delete
        mcs0.delete_atom( atoms_to_delete )
        mcs0 = mcs0.copy()
        partial_ring         = self._delete_broken_ring( mol0, mol1, mcs0 )
        mcs0 = mcs0.copy()
        atoms_to_delete_2 = []
        for e in mcs0.molecules()[1:] :
            atoms_to_delete_2.extend( e )
        #print "Rule Atom to delete", atoms_to_delete
        mcs0.delete_atom( atoms_to_delete_2 )
        smiles0 = mcs0.smiles()
        smiles1 = smiles0

        KBASE.deposit_extra( mcs_id, "trimmed-mcs",  {id0:smiles0, id1:smiles1,} )
        KBASE.deposit_extra( mcs_id, "partial_ring", len( partial_ring ) )

        num_heavy_atoms = len( mcs0.heavy_atoms() )
        num_light_atoms = len( mcs0.atom ) - num_heavy_atoms

        KBASE.deposit_extra( mcs_id, "num_heavy_atoms", num_heavy_atoms )
        KBASE.deposit_extra( mcs_id, "num_light_atoms", num_light_atoms )

        return similarity.exp_delta( 2 * (orig_num_heavy_atoms - num_heavy_atoms), 0 )



# Example of a complex rule: A combination of a few simple rules (in case, they are Mcs, MinimumNumberOfAtom, and Cutoff).
# cutoff_simi = Cutoff( 0.2, Mcs( MinimumNumberOfAtom( 4 ) ) )



if ("__main__" == __name__) :
    import struc
    from mcs   import SchrodMcs
    from kbase import KBASE
    
    filenames = ["xfer3.11.mol2", "xfer3.12.mol2",]
    id_list   = struc.read_n_files( filenames )
    mol0      = KBASE.ask( id_list[0] )
    mol1      = KBASE.ask( id_list[1] )
    mcs       = SchrodMcs( 3 )
    mcs_id    = mcs.search( mol0, mol1 )[0]
    mol_id    = KBASE.ask( mcs_id, "mcs_parents" )

    print MCS.similarity( mol_id[0], mol_id[1] )
    
