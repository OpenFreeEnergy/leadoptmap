"""Defines a `Struc' class as a generic represention of molecular structure
"""



from kbase import KBASE

import os
import hashlib



# Flag indicating which infrastructure we are using.
infrastructure = None    # "schrodinger" | "oechem"



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



    def smarts( self ) :
        """
        Returns a SMARTS string of this structure.
        """
        raise NotImplementedError( "`smarts' method not implemented by subclass" )



    def smiles( self ) :
        """
        Returns a SMILES string of this structure.
        """
        raise NotImplementedError( "`smiles' method not implemented by subclass" )

        

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
            return self._struc.atom[atom_index].chirality in ["R", "S", "ANR", "ANS"]
        
        
        
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
            


        def smarts( self ) :
            """
            Returns a SMARTS string for this structure.
            """
            import schrodinger.structutils.analyze as analyze

            return analyze.generate_smarts( self._struc )


            
        def smiles( self ) :
            """
            Returns a SMILES string for this structure.
            """
            import schrodinger.structutils.analyze as analyze

            return analyze.generate_smiles( self._struc )
            
            
        
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

    infrastructure = "schrodinger"
    
except ImportError, e :
    pass



try:
    
    import openeye.oechem as oechem 
    class OeStruc( Struc ) :
        """
        A `Struc' subclass based on Openeye OEMol's infrastructure
        """

        def __init__( self, struc ) :
            
            Struc.__init__( self )
            self._struc = struc

            self.atom = {}
            for atom in self._struc.GetAtoms():
                oe_idx = atom.GetIdx() + 1
                if self.atom.has_key(oe_idx):
                    print "Struc has duplicate atom index : %s need to check"%(oe_idx)
                else:
                    self.atom[oe_idx] = atom


        def copy( self ) :
            return OeStruc( self._struc.CreateCopy() )

        

        def extract(self) :
            """                                               
            Return a new structure object which contains the atoms of the current structure that appear in the specified list.
            """                                               
            new_mol = self._struc.copy()
            for atom in new_mol.GetAtoms():                   
                oe_idx = atom.GetIdx() + 1                    
                if oe_idx not in indices:                     
                    new_mol.DeleteAtom(atom)                  
            return OeStruc(new_mol)

        
            

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
            for e in self.atom.values():                      
                if not e.IsHydrogen():                        
                    oe_idx = e.GetIdx() + 1
                    ret.append( oe_idx )                      
            return ret
        

        def total_charge( self ) :
            """
            Returns net charge of the structure
            """
            return OENetCharge( self._struc )

        

        def is_chiral_atom( self, atom_index ) :
            """                                               
            Returns true if the atom indicated by C{atom_index} is chiral; otherwise, false.
            """                                               
            return self.atom[atom_index].IsChiral()
        

        def chiral_atoms( self ) :
            """
            Returns the indices of chiral atoms.
            """
            ret = []                                          
            ret_atoms = []                                    
            oechem.OEPerceiveChiral(self._struc)              
            for e in self._struc.GetAtoms():                  
                ret_atoms.append(e)                           
            for atom in ret_atoms:                            
                if atom.IsChiral():                           
                    oe_idx = atom.GetIdx() + 1                
                    ret.append( oe_idx )                      
            return ret
        

        def ring_atoms( self ) :
            """
            Returns a set of ring atoms.
            """
            ret = []                                          
            ret_atoms = []                                    
            for e in self._struc.GetAtoms():                  
                ret_atoms.append(e)                           
            for atom in ret_atoms:                            
                if atom.IsInRing():                           
                    oe_idx = atom.GetIdx() + 1
                    ret.append( oe_idx )                      
            return set( ret )


        def bonded_atoms( self, atom_index ) :
            """
            Returns a list of atom indices of atoms bonded to the indicated atom. 
            """
            ret = []                                          
            #get atom object first to avoid interator problem 
            ret_atoms = []                                    
            for e in self.atom[atom_index].GetAtoms():        
                ret_atoms.append(e)                           
            for atom in ret_atoms:                            
                oe_idx = atom.GetIdx() + 1                    
                ret.append(atom.GetIdx())                     
            return ret


        def delete_atom( self, atom_index ) :
            """
            Deletes an atom.
            """
            if (not isinstance( atom_index, list )) :
                atom_index = [atom_index,]
            ret_atoms = []
            for (idx, e) in enumerate(self._struc.GetAtoms()):
                oe_idx = idx + 1
                for i in atom_index:
                    if i == oe_idx:
                        ret_atoms.append(e)

            for j in ret_atoms:                               
                self._struc.DeleteAtom(j)                     
                self.atom = {}                                
                for atom in self._struc.GetAtoms():           
                    oe_idx = atom.GetIdx() + 1                
                    if self.atom.has_key(oe_idx):
                        print "Struc has duplicate atom index :i %s need to check"%oe_idx                                                     
                    else:                                     
                        self.atom[oe_idx] = atom


        def write (self , filename, format = "mol2", mode = "w"):
            #print "Struc check write format", oechem.oemolostream(filename)
            if mode == "a":                                   
                raise ValueError("OeStruc write doesn't support append method")
            elif mode == "w":                                 
                if format == "pdb":                           
                    ofs = oechem.oemolostream(filename + '.pdb')
                    ofs.SetFormat(oechem.OEFormat_PDB)        
                elif format == "mol2":                        
                    ofs = oechem.oemolostream(filename + '.mol2')
                    ofs.SetFormat(oechem.OEFormat_MOL2)       
                elif format == "xyz":                         
                    ofs = oechem.oemolostream(filename + '.xyz')
                    ofs.SetFormat(oechem.OEFormat_XYZ)        
            else:                                             
                raise ValueError( "Invalid value for `mode' argument: '%s', should be one of 'a' and 'w'.")
            #format = OEFormat_MOL2                           
            #ofs.SetFormat(OEFormat_MOL2)                     
            #return oechem.OEWriteMol2File(ofs, self._struc  )
            return oechem.OEWriteMolecule(ofs, self._struc  )
        
            
    def read_file_oe (filename):
        """
        Reads a .mol2 file and returns title of molecule , base molecule object and `OEMol' objects.
        """
        istream = oechem.oemolistream()
        istream.open(filename)

        molecule = oechem.OEMol()

        oechem.OEReadMolecule(istream, molecule)

        istream.close()                                       
                                                              
        return molecule


    def read_n_files_oe( filenames ) :
        """
        `filenames' is a list of file names. The format of each file will be determined from the file's extension name.
        Reads the files and deposits them into the `KBASE'. Returns a list of keys.
        """                                                   
        strucid = []                                          
        for fn in filenames :                                 
            strucs = read_file_oe( fn )                         
            e = OeStruc(strucs)                                
            id = KBASE.deposit( e.id(), e )               
            e.set_id( id )                                
            strucid.append( id )
        return strucid        

    infrastructure = "oechem"

except ImportError, e :
    pass



if (infrastructure is None) :
    print "ERROR: Need either Schrodinger's or OEChem's infrastructure to run, but none is found."
    import sys
    sys.exit( 1 )

    


if ("__main__" == __name__) :
    filenames = ["xfer3.10.mol2", "xfer3.11.mol2",]
    id_list = read_n_files( filenames )
    mol0 = KBASE.ask( id_list[0] )
    print mol0.title(), len( mol0.heavy_atoms() )
    mol1 = KBASE.ask( id_list[1] )
    print mol1.title(), len( mol1.heavy_atoms() )
