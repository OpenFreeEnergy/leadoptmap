"""Defines a `Struc' class as a generic represention of molecular structure
"""



from kbase import KBASE

import os
import hashlib



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



    def title( self ) :
        raise NotImplementedError( "`title()' method not implemented by subclass" )
    

        
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
    import schrodinger.structure as structure
    
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



        def title( self ) :
            """
            Returns the title of this structure. (Normally title's a user-friendly description)
            """
            return self._struc.title
        
            
            
        def heavy_atoms( self ) :
            """
            Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
            """
            ret = []
            for e in self.atom :
                if (e.atomic_number > 1) :
                    ret.append( int( e ) )
            return ret



        def total_charge( self ) :
            """
            Returns the formal charge of the structure
            """
            return self._struc.formal_charge


        
        def write( self, filename, format = None, mode = "a" ) :
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
                self._struc.write( filename, format )
            else :
                raise ValueError( "Invalid value for `mode' argument: '%s', should be one of 'a' and 'w'." )

            
            
    def read_file( filename, format = None ) :
        """
        Reads a .mae file and returns a list of `ScrhodStruc' objects.

        @type  filename: C{str}
        @param filename: Name of the structure file
        @type  format  : C{str} or C{None}
        @param format  : Specify the format of the file. It must be one of the following case-sensitive strings: "pdb", "sd",
                         "mol2", and "maestro". If its value is C{None}, the format will be determined from extension name.
        """
        ret = []
        for ct in structure.StructureReader( filename, format = format ) :
            ret.append( SchrodStruc( ct ) )
        return ret



    def read_n_files( filenames ) :
        """
        `filenames' is a list of file names. The format of each file will be determined from the file's extension name. Reads
        the files and deposits them into the `KBASE'. Returns a list of keys.
        """
        strucid = []
        for fn in filenames :
            strucs = read_file( fn )
            for e in strucs :
                id = hashlib.sha1( e.title() ).hexdigest()
                strucid.append( KBASE.deposit( id, e ) )
        return strucid
    
except ImportError, e :
    raise e
    pass




if ("__main__" == __name__) :
    filenames = ["xfer3.10.mol2", "xfer3.11.mol2",]
    id_list = read_n_files( filenames )
    mol0 = KBASE.ask( id_list[0] )
    print mol0.title(), len( mol0.heavy_atoms() )
    mol1 = KBASE.ask( id_list[1] )
    print mol1.title(), len( mol1.heavy_atoms() )
