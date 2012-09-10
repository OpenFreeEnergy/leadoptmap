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

        # Private attributes:
        self._id = None


        
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



    def title( self ) :
        raise NotImplementedError( "`title()' method not implemented by subclass" )


    
    def set_title( self, new_title ) :                        
        raise NotImplementedError( "`set_title' method not implemented by subclass" )


        
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

    


try:
    from openeye.oechem import *
    from mmtools.moltools.ligandtools import *
 
    class OeStruc( Struc ) :
        """
        A `Struc' subclass based on Openeye OEMol's infrastructure
        """

        def __init__( self, struc ) :
            
            Struc.__init__( self )
            self._struc = struc

            #self.atom = self._struc.GetAtoms()
            #transfer from iterator type to property 
            self.atom = []
            for atom in self._struc.GetAtoms():
                self.atom.append( atom )



        def copy( self ) :
            return OeStruc( self._struc.CreateCopy() )

        

        def extract(self) :
            """
            need to add
            """
            pass
        
            

        def title( self ) :
            """     
            Returns the title of this structure. (Normally title's a user-friendly description)
            """
            return self._struc.GetTitle()

        
        
        def set_title( self, new_title ) :
            """                                               
            Sets a new title to this structure.
            """
            self._struc.SetTitle( new_title )

            

        def heavy_atoms(self) :
            """
            Returns a list of indices of heavy atoms (viz non-hydrogen atoms).
            """
            ret = []
            for e in self.atom:
                if not e.IsHydrogen():
                    ret.append( e.GetIdx() )
            return ret

        

        def total_charge( self ) :
            """
            Returns net charge of the structure
            """
            return OENetCharge( self._struc )

        

        def is_chiral_atom( self, atom_index ) :
            """
            need to add
            """
            pass

        

        def chiral_atoms( self ) :
            """
            Returns the indices of chiral atoms.
            """
            ret = []
            OEPerceiveChiral( self._struc )
            for e in self._struc.GetAtoms() :
                if (e.IsChiral()) :
                    ret.append( e.GetIdx() )
            return ret

        

        def ring_atoms( self ) :
            """
            Returns a set of ring atoms.
            """
            ret = []
            for e in self._struc.GetAtoms() :
                if (e.IsInRing()) :
                    ret.append( e.GetIdx() )
            return set( ret )



        def bonded_atoms( self, atom_index ) :
            """
            Returns a list of atom indices of atoms bonded to the indicated atom. 
            """
            ret = []
            for e in self._struc.atom[atom_index].GetAtoms() :
                ret.append( e.GetIdx() )
            return ret



        def delete_atom( self, atom_index ) :
            """
            Deletes an atom.
            """
            if (not isinstance( atom_index, list )) :
                atom_index = [atom_index,]
            atoms = []
            for i, e in enumerate( self._struc.GetAtoms() ) :
                if i in (atom_index) :
                    atoms.append( e )
            for a in atoms :
                self._struc.DeleteAtom( a )



        def write( self, filename, format, mode = "a" ) :
            """
            Writes this structure into a mol2 file.
            """
            return OEWriteMol2File( oemolostream( filename ), self._struc )

        
            
    def read_file_oe (filename):
        """
        Reads a .mol2 file and returns title of molecule , base molecule object and `OEMol' objects.
        """
        mol = readMolecule(filename)
        title = mol.GetTitle()
        oemol = OEMol( mol )
        return (title ,mol, oemol)



    def read_n_files_oe( filenames ) :
        """
        `filenames' is a list of file names. The format of each file will be determined from the file's extension name. Reads the files and deposits them into the `KBASE'. Returns a list of keys.
        """                                                   
        strucid = []                                          
        for fn in filenames :                                 
            strucs = read_file_oe( fn )[-1]                          
            e = OeStruc(strucs)                                
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
