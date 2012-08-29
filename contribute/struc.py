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
        # Public attributes:
        self.atom = None

        # Private attributes:
        self._id = None

        

    def __str__( self ) :
        return self.title()
    

        
    def _atom( self, index ) :
        """
        Returns the `index'-th atom.
        """
        raise NotImplementedError( "`_atom()' method not implemented by subclass" )



    def id( self ) :
        """
        Returns the ID of this structure.
        """
        if (self._id is None) :
            self._id = hashlib.sha1( self.title() ).hexdigest()
        return self._id



    def set_id( self, id ) :
        """
        Sets the ID for this structure.
        """
        self._id = id



    def copy( self ) :
        raise NotImplementedError( "`copy' method not implemented by subclass" )

    

    def extract( self, indices ) :
        """
        Return a new structure object which contains the atoms of the current structure that appear in the specified list.
        """
        raise NotImplementedError( "`extract' method not implemented by subclass" )
    
    

    def title( self ) :
        raise NotImplementedError( "`title' method not implemented by subclass" )
    
    

    def set_title( self, new_title ) :
        raise NotImplementedError( "`set_title' method not implemented by subclass" )
    

        
    def heavy_atoms( self ) :
        """
        Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
        """
        raise NotImplementedError( "`heavy_atoms' method not implemented by subclass" )



    def is_chiral_atom( self, atom_index ) :
        """
        Returns true if the atom indicated by C{atom_index} is chiral; otherwise, false.

        @type  atom_index: C{int}
        @param atom_index: Atom index
        """
        raise NotImplementedError( "`is_chiral_atom' method not implemented by subclass" )
    
    
    
    def chiral_atoms( self ) :
        """
        Returns the indices of the chiral atoms.

        @rtype : C{list} of C{int}
        @return: A list of atom indices
        """
        raise NotImplementedError( "`chiral_atoms' method not implemented by subclass" )



    def ring_atoms( self ) :
        """
        Returns a set of ring atoms.

        @rtype : C{set} of C{int}
        @return: A set of atom indices
        """
        raise NotImplementedError( "`chiral_atoms' method not implemented by subclass" )



    def bonded_atoms( self, atom_index ) :
        """
        Returns a list of atom indices of atoms bonded to the indicated atom.
        
        @type  atom_index: C{int}
        @param atom_index: A single index or a list of indices of the atoms to be deleted
        
        @rtype : C{list} of C{int}
        @return: A list of atom indices of atoms bonded to the indicated atom
        """
        raise NotImplementedError( "`bonded_atoms' method not implemented by subclass" )



    def total_charge( self ) :
        """
        Returns the total charge of the structure.
        """
        raise NotImplementedError( "`total_charge' method not implemented by subclass" )



    def delete_atom( self, atom_index ) :
        """
        Deletes a atom.

        @type  atom_index: C{int} or C{list} of C{int}
        @param atom_index: A single index or a list of indices of the atoms to be deleted
        """
        raise NotImplementedError( "`delete_atoms' method not implemented by subclass" )



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
            Struc.__init__( self )
            self._struc = struc

            # Public attributes:
            self.atom = self._struc.atom
            self.bond = self._struc.bond



        def copy( self ) :
            """
            Returns a copy of this structure.
            """
            return SchrodStruc( self._struc.copy() )



        def extract( self, indices ) :
            """
            Return a new structure object which contains the atoms of the current structure that appear in the specified list.
            """
            return SchrodStruc( self._struc.extract( indices, True ) )

        
        
        def title( self ) :
            """
            Returns the title of this structure. (Normally title's a user-friendly description)
            """
            return self._struc.title
        
            
            
        def set_title( self, new_title ) :
            """
            Sets a new title to this structure.
            """
            self._struc.title = new_title
        
            
            
        def heavy_atoms( self ) :
            """
            Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
            """
            ret = []
            for e in self.atom :
                if (e.atomic_number > 1) :
                    ret.append( int( e ) )
            return ret



        def is_chiral_atom( self, atom_index ) :
            """
            Returns true if the atom indicated by C{atom_index} is chiral; otherwise, false.

            @type  atom_index: C{int}
            @param atom_index: Atom index
            """
            return atom.chirality in ["R", "S", "ANR", "ANS"]
        
        
        
        def chiral_atoms( self ) :
            """
            Returns the indices of the chiral atoms.

            @rtype : C{list} of C{int}
            @return: A list of atom indices
            """
            ret = []
            for atom in self._struc.atom :
                if (atom.chirality in ["R", "S", "ANR", "ANS"]) :
                    ret.append( int( atom ) )
            return ret



        def ring_atoms( self ) :
            """
            Returns a set of ring atoms.

            @rtype : C{set} of C{int}
            @return: A set of atom indices
            """
            rings = self._struc.find_rings()
            ret   = []
            for e in rings :
                ret.extend( e )
            return set( ret )



        def bonded_atoms( self, atom_index ) :
            """
            Returns a list of atom indices of atoms bonded to the indicated atom.
            
            @type  atom_index: C{int}
            @param atom_index: A single index or a list of indices of the atoms to be deleted
            
            @rtype : C{list} of C{int}
            @return: A list of atom indices of atoms bonded to the indicated atom
            """
            ret = []
            for e in self._struc.atom[atom_index].bonded_atoms :
                ret.append( int( e ) )
            return ret
        
                
        
        def total_charge( self ) :
            """
            Returns the formal charge of the structure
            """
            return self._struc.formal_charge



        def delete_atom( self, atom_index ) :
            """
            Deletes a atom.

            @type  atom_index: C{int} or C{list} of C{int}
            @param atom_index: A single index or a list of indices of the atoms to be deleted
            """
            if (not isinstance( atom_index, list )) :
                atom_index = [atom_index,]
            self._struc.deleteAtoms( atom_index )
            
            
        
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
                id = KBASE.deposit( e.id(), e )
                e.set_id( id )
                strucid.append( id )
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
