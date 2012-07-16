"""Maximum common substructure searching (MCSS) engine
"""



import subprocess



class Mcs( object ) :
    def __init__( self ) :
        """

        """
        pass



    def search( self, mol0, mol1 ) :
        """

        """
        raise NotImplementedError( "`search' method not implemented in subclass" )

    

    def search_all( self, mols ) :
        """

        """
        raise NotImplementedError( "`search' method not implemented in subclass" )
        


try :
    import openeye.oechem as oechem

    class OeMcs( Mcs ) :
        def __init__( self, atom_expr = oechem.OEExprOpts_EqONS,
                      bond_expr = oechem.OEExprOpts_BondOrder|oechem.OEExprOpts_EqSingleDouble|oechem.OEExprOpts_EqAromatic,
                      ringchecking = 'Strict', min_num_atoms, is_approximate = True ) :
            self._atom_expr      = atom_expr
            self._bond_expr      = bond_expr
            self._ringchecking   = ringchecking
            self._min_num_atoms  = min_num_atoms
            self._is_approximate = is_approximate
            


        def search( self, mol0, mol1 ) :
            """
            Returns a `Struc' object, which is maximum common substructure.
            """
            import MCSS_tool
            ret = MCSS_tool.determineMCSS_requireRings( mol1, mol0, debug = False, min_atoms = self._min_num_atoms,
                                                        atomexpr = self._atom_expr, bondexpr = self._bond_expr,
                                                        ringchecking = self._ringchecking, approximate = self._is_approximate )
            # Maybe does some other things.
            
            return ret[0]
        
        

        def search_all( self, mols ) :
            """
            O(N^2). Needs optimization.
            Returns a dictionary: key = pair of molecule indices in the given list, value = `Struc' object.
            """
            ret     = {}
            num_mol = len( mols )
            for i in range( num_mol ) :
                for j in range( i ) :
                    ret[(i, j,)] = self.search( mols[i], mols[j] )
            return ret
        
except ImportError :
    pass



try :
    class SchrodMcs( Mcs ) :
        def __init__( self ) :
            self._cmd = os.path.join( os.environ['SCHRODINGER'], "utilities", "canvasMCS" )
            
            

        def search( self, mol0, mol1 ) :
            return self.search( [mol0, mol1,] )



        def search_all( self, mols ) :
            tmp_fname = "__temp_file_ok_to_delete_after_running__.mae"
            out_fname = "__temp_file_ok_to_delete_after_running__.cvs"

            for mol in mols :
                mol.write( temp_fname )
                mol.write( temp_fname )
            cmd          = [self._cmd,
                            "-imae",     tmp_fname,
                            "-opw",      out_fname,
                            "-atomtype", "11",
                            ]
            mcs_proc     = subprocess.Popen( cmd, stderr = subprocess.STDOUT, stdout = out_fname )
            null, stderr = mcs_proc.communicate()
            val          = mcs_proc.returncode

            if (val == 17) :
                raise RuntimeError( "Used a MCS feature that requires Schrodinger's CANVAS_ELEMENTS license." )
            if (val != 0 ) :
                msg = "CanvasMCS exited prematurely. This could be because the input molecules were too dissimilar" \
                      " or too numerous, or because the chosen atom-typing scheme was too general."
                with open( out_fname ) as fh:
                    msg += "\n\n"
                    msg += fh.read()
                raise RuntimeError( msg )


except ImportError :
    pass

        


            
