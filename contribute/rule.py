"""Defines a rule engine.
"""



import kbase



class Rule( object ) :
    """
    Base class of all rules.
    """
    def __init__( self, kbase, *subrules ) :
        self._kbase    = kbase
        self._subrules = *subrules



    def kbase( self ) :
        return self._kbase



    def set_kbase( self, kbase ) :
        self._kbase = kbase

        

    def _similarity( self, *arg, **kwarg ) :
        """

        """
        raise NotImplementedError( "`_similarity' method not implemented in subclass" )


    
    def similarity( self, *arg, **kwarg ) :
        """
        Returns a floating number in the range of [0, 1].
        """
        result = self._similarity( *arg, **kwarg )
        if (result > 0) :
            for e in self._subrules :
                if (isinstance( e, list )) :
                    result *= max( [e.similarity( *arg, **kwarg )] )
                else :
                    result *= e.similarity( *arg, **kwarg )
        return result



class MinimumNumberOfAtom( Rule ) :
    """
    Rule on mininum number of atoms
    """
    def __init__( self, threshold = 4, heavy_only = False, *subrules ) :
        Rule.__init__( self, *subrules )
        self._threshold  = threshold
        self._heavy_only = heavy_only


        
    def _similarity( self, struc ) :
        if (self._heavy_only) :
            num_atom = len( struc.heavy_atoms() )
        else :
            num_atom = struc.atoms()
        return float( num_atom >= self._threshold )



class Mcss( Rule ) :
    """
    
    """
    def __init__( self, subrules ) :
        Rule.__init__( self, *subrules )

        

    def _similarity( self, mol0, mol1 ) :
        raise NotImplementedError( "`_similarity' method not implemented in subclass of `Mcss'" )



class SchrodMcss( Mcss ) :
    """

    """
    def __init__( self, subrules ) :
        Mcss.__init__( self, *subrules )



    def _similarity( self, mol0, mol1 ) :
        pass
