"""Defines a simple dictionary-based knowledge base.
"""



import hashlib



class Kbase( object ) :
    def __init__( self ) :
        self._knowledge = {}
        self._extra     = {}



    def ask( self, quest, tag = None ) :
        if (tag) :
            try :
                extra = self._extra[quest]
            except KeyError :
                raise LookupError( "Ignorance on %s" % quest )
            try :
                return extra[tag]
            except KeyError :
                raise LookupError( "Ignorance on %s with tag %s" % (quest, tag,) )
        try :
            return self._knowledge[quest]
        except KeyError :
            raise LookupError( "Ignorance on %s" % quest )



    def deposit( self, key, knowlet, should_overwrite = False ) :
        if (not should_overwrite) :
            while (key in self._knowledge) :
                key = hashlib.sha1( key ).hexdigest()
        self._knowledge[key] = knowlet
        return key

    

    def deposit_extra( self, key, tag, knowlet ) :
        if (key not in self._knowledge) :
            raise KeyError( "'%s' not found in the knowledge base." )
        if (key not in self._extra) :
            self._extra[key] = {}
        self._extra[key][tag] = knowlet

        

KBASE = Kbase()



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
