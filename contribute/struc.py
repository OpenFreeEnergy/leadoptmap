"""Defines a `Struc' class as a generic represention of molecular structure
"""



import os



class _AtomContainer( object ) :
    """

    """
    def __init__( self, struc ) :
        """
        Initializes the container.
        """
        self._struc = struc

        
        
    def __getitem__( self, index ) :
        """
        Index is 1-based.
        """
        return self._struc._atom( index )

    

    def __iter__(self):
        """
        Iteration access
        """
        for i in range( 1, len( self ) + 1 ) :
            yield _struc._atom( i )

            

    def __len__(self):
        """
        Returns the number of atoms.
        """
        return self._struc.num_atom()



class Struc( object ) :
    """

    """
    def __init__( self ) :
        self.atom = _AtomContainer( self )


        
    def _atom( self, index ) :
        """
        Returns the `index'-th atom.
        """
        raise NotImplementedError( "`_atom()' method not implemented by subclass" )


        
    def heavy_atoms( self ) :
        """
        Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
        """
        raise NotImplementedError( "`heavy_atoms' method not implemented by subclass" )



    def total_charge( self ) :
        """
        Returns the total charge of the structure.
        """
        raise NotImplementedError( "`total_charge' method not implemented by subclass" )



    def write( filename, format, mode = "a" ) :
        """
        Writes this structure into a file in the designated format.

        @type  mode: C{char}, 'a' | 'w'
        @param mode: When a file of the same name exists, this determines whether to overwrite ('w') or append ('a') to the
                     file.
        """
        raise NotImplementedError( "`write' method not implemented by subclass" )

    

try :
    import schrodinger.application.structure as structure
    
    class SchrodStruc( Struc ) :
        """
        A `Struc' subclass based on Schrodinger's infrastructure
        """
        def __init__( self, struc ) :
            """
            `struc' should be a `schrodinger.structure.Structure' object.
            """
            self._struc = struc

            # Public attributes:
            self.atom = self._struc.atom
            self.bond = self._struc.bond

            
            
        def heavy_atoms( self ) :
            """
            Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
            """
            ret = []
            for e in self.atom :
                if (e.atomic_number > 1) :
                    ret.append( int( e ) )



        def total_charge( self ) :
            """
            Returns the formal charge of the structure
            """
            return self._struc.formal_charge


        
        def write( filename, format = None, mode = "a" ) :
            """
            Writes this structure into a file in the designated format.

            @type  format: C{str} or C{None}
            @param format: If its value is C{None}, the file format is determined from the filename suffix. If specified, it
                           must be one of the following case-sensitive strings: "pdb", "mol2", "sd", "maestro", "smiles", and
                           "smilescsv".
            """
            if   ('a' == mode) :
                if (os.path.isfile( filename )) : self._struc.append( filename, format )
                else                            : self._struc.write ( filename, format )
            elif ('w' == mode) :
                self._struc.write ( filename, format )
            else :
                raise ValueError( "Invalid value for `mode' argument: '%s', should be one of 'a' and 'w'." )

            
            
    def read_maestro_file( filename ) :
        """
        Reads a .mae file and returns a list of `ScrhodStruc' objects.
        """
        ret    = []
        reader = structure.MaestroReader( filename )
        try :
            while (True) :
                ct = reader.next()
                ret.append( SchrodStruc( ct ) )
        except :
            pass
        return ret
        
except ImportError :
    pass

