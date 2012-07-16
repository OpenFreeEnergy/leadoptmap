"""Defines a simple dictionary-based knowledge base.
"""



class Kbase( object ) :
    def __init__( self ) :
        self._knowledge = {}



    def ask( self, quest ) :
        try :
            return self._knowledge[quest]
        except KeyError :
            raise LookupError( "Ignorance on %s" % quest )



    def deposit( self, key, knowlet ) :
        self._knowledge[key] = knowlet
        


def create_mcss_kbase( mols ) :
    """
    Creates a MCSS matrix for a given list of molecules. Stores the matrix into a C{Kbase} with key being a pair of molecule
    indices and value being a C{Struc} object of the common substructure between the two molecules.
    This function returns a C{Kbase} object containing the MCSS matrix data.
    
    @type  mols: C{list}
    @param mols: A list of molecular structures as C{Struc} objects
    """
    # Implementation to be added
    pass
