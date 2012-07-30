"""Defines a rule engine.
"""



from kbase import KBASE

import similarity

import hashlib



class Rule( object ) :
    """
    Base class of all rules.
    """
    def __init__( self, *subrules ) :
        self._subrules = subrules



    def _similarity( self, id0, id1 ) :
        """
        
        """
        raise NotImplementedError( "`_similarity' method not implemented in subclass" )


    
    def similarity( self, id0, id1 ) :
        """
        Returns a floating number in the range of [0, 1].
        """
        result = self._similarity( id0, id1 )
        if (result > 0) :
            for e in self._subrules :
                if (isinstance( e, list )) :
                    result *= max( [e.similarity( id0, id1 )] )
                else :
                    result *= e.similarity( id0, id1 )
        return result



class MinimumNumberOfAtom( Rule ) :
    """
    Rule on mininum number of atoms
    """
    def __init__( self, threshold = 4, heavy_only = False, *subrules ) :
        Rule.__init__( self, *subrules )
        self._threshold  = threshold
        self._heavy_only = heavy_only


        
    def _similarity( self, id0, id1 ) :
        try :
            mcs = KBASE.ask( hashlib.sha1( id0 + id1 ).hexdigest() )
        except LookupError :
            mcs = KBASE.ask( hashlib.sha1( id1 + id0 ).hexdigest() )

        # Uses the first common substructure.
        mcs = mcs[0]
        
        if (self._heavy_only) :
            num_atom = len( mcs.heavy_atoms() )
        else :
            num_atom = mcs.atom
        return float( num_atom >= self._threshold )



class Mcs( Rule ) :
    """
    
    """
    def __init__( self, *subrules ) :
        Rule.__init__( self, *subrules )

        

    def _similarity( self, id0, id1 ) :
        try :
            mcs = KBASE.ask( hashlib.sha1( id0 + id1 ).hexdigest() )
        except LookupError :
            mcs = KBASE.ask( hashlib.sha1( id1 + id0 ).hexdigest() )

            # Swaps the id0 and id1 values. We presumes the order of them matters. id0 should be the reference molecule.
            id2 = id0
            id0 = id1
            id1 = id2
            
        # Uses the first common substructure.
        mcs = mcs[0]

        # Retrieves the parent molecules of the MCS.
        mol0 = KBASE.ask( id0 )
        mol1 = KBASE.ask( id1 )

        return similarity.by_heavy_atom_count( mol0, mol1, mcs )



MCS_RULE = Mcs( MinimumNumberOfAtom() )



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

    print MCS_RULE.similarity( mol_id[0], mol_id[1] )
    
