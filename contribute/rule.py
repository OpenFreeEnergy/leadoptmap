"""We implement a minimalist rule engine, and a few basic rule classes. With these infrastructure, we define a C{MCS} rule.
"""



from kbase import KBASE

import similarity

import hashlib



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
    Rule on minimum number of atoms in the maximum common substructure

    Similarity score is 0 if the minimum number of atoms in the maximum common substructure is less than a specified threshold
    number, or 1 if otherwise.
    """
    def __init__( self, threshold = 4, heavy_only = False, *subrules ) :
        """
        @type   threshold: C{int}
        @param  threshold: Threshold of number of atoms. Above this threshold, the similarity score is 1; otherwise it's 0.
        @type  heavy_only: C{bool}
        @param heavy_only: If the value is C{True}, we count only the heavy atoms; if it's C{False}, we count all atoms.
        """
        Rule.__init__( self, *subrules )
        self._threshold  = threshold
        self._heavy_only = heavy_only


        
    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        mcs0 = KBASE.ask( kwarg["mcs_id"] )[0]

        if (self._heavy_only) :
            num_atom = len( mcs0.heavy_atoms() )
        else :
            num_atom = mcs0.atom
        return float( num_atom >= self._threshold )



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
        simi = Rule.similarity( self, id0, id1, **kwarg )
        if (simi < self._cutoff) :
            simi = 0.0
        return simi



class Mcs( Rule ) :
    """
    MCS-based rule
    Similarity is scored using the C{similarity.by_heavy_atom_count} (see the C{similarity} module).
    """
    def __init__( self, *subrules ) :
        Rule.__init__( self, *subrules )

        

    def _similarity( self, id0, id1, **kwarg ) :
        # Uses the first common substructure.
        mcs0 = KBASE.ask( kwarg["mcs_id"] )[0]

        # Retrieves the parent molecules of the MCS.
        mol0 = KBASE.ask( id0 )
        mol1 = KBASE.ask( id1 )

        return similarity.by_heavy_atom_count( mol0, mol1, mcs0 )



# MCS rule: A combination of the following rules: Mcs, MinimumNumberOfAtom, and Cutoff.
MCS = Cutoff( 0.2, Mcs( MinimumNumberOfAtom() ) )



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
    
