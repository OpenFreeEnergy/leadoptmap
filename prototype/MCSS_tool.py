#!/usr/bin/env python


"""Tools relating to scoring similarities between molecules based on the size of their (modified) maximum common substructures. The main functionality here is provided by determineMCSS_requireRings; see the documentation on that function for additional information.

All code here
Copyright (c) 2012, the University of New Orleans, written by David L. Mobley and Shuai Liu
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

   * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
   * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
   * Neither the name of the University of New Orleans nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

#IDEA: Require good matches to have at least one ring in common; the more atoms they have in addition to the ring the better. NOTE that it is not necessary that all of the atoms in the ring be of the same type, but there must be the same number.

from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
from openeye.oechem import *
import numpy
import pickle
import glob

def delete_unconnected_atoms( partialmol, debug = False):
    """Helper function. Take a (potentially partial) OE molecule after some atoms have been deleted, and clean it up -- any 'hanging' groups which are not connected to the rest of the molecule are removed and only a single connected molecule is returned. The main group of atoms is determined as the largest connected group of atoms; all smaller groups are deleted."""

    #Print initial atoms list
    if debug:
        print "Initial molecule consists of atoms:"
        atomlist = '   '
        for atom in partialmol.GetAtoms():
            atomlist += atom.GetName() + ' '
        print atomlist


    #Delete everything we don't want
    if partialmol.NumAtoms() > 0:
        #Build a list of connected groups of atoms 
        ngroups = 0

        #Begin doing this by enumerating all connections
        connections = {}
        for atom in partialmol.GetAtoms():
            for nbr in atom.GetAtoms():
                if not connections.has_key( atom ): connections[atom] = [ nbr ]
                else: connections[atom].append( nbr)

                if not connections.has_key( nbr ): connections[nbr] = [atom ]
                else: connections[nbr].append(atom)
        
        #Now build connected groups -- tracking visited atoms
        visited = []
        connected_groups = []
        for atom in partialmol.GetAtoms():
            #print atom.GetName(), visited
            #If we have not already visited this atom, it starts a new connected cluster and we need to visit all its neighbors, etc.
            if not atom in visited:
                clusidx = len( connected_groups )
                connected_groups.append( [ atom ] )
           
                #Track this as visited
                visited.append( atom )

                #Track connections to visit
                if connections.has_key( atom ):
                    to_visit = connections[ atom ]
                else: to_visit = []

                #Visit connections until there are no connected ones left to visit
                while len( to_visit ) > 0:
                    #Visit a new node
                    nbr = to_visit.pop()
                    #If we are seeing something new, add it to the visited list and add its nbrs to the list to visit, and add it to our cluster
                    if not nbr in visited:
                        connected_groups[clusidx].append( nbr )
                        visited.append( nbr )
                        for nbr2 in nbr.GetAtoms():
                            if not nbr2 in to_visit:
                                to_visit.append( nbr2 )
                        

            #Otherwise, we are done with this atom since we already visited it
                
        #Print debugging info
        if debug:
             print "Connected groups:"
             for group in connected_groups:
                 outstr = ''
                 for atom in group:
                     outstr += atom.GetName() + ' '
                 print outstr
        #If there is only one connected group, we are done
        #Otherwise we have to delete atoms
        if len(connected_groups) > 1:
            #Figure out what group to keep -- largest group
            grouplens = [ len(group) for group in connected_groups ]
            maxgroup = max(grouplens)
            keepgroup = grouplens.index( maxgroup ) #Index of group with max length

            #Delete atoms in the other (smaller) groups
            for (idx, group) in enumerate(connected_groups):
                if not idx==keepgroup:
                    if debug: print "   Deleting group %s with %s atoms..." % (idx, len(group))
                    for atom in group:
                        partialmol.DeleteAtom( atom )


    #When done, return molecule
    return OEMol(partialmol)


def determineMCSS_requireRings(ligand1, ligand2, debug = False, min_atoms = 4, atomexpr = OEExprOpts_EqONS, bondexpr = OEExprOpts_BondOrder | OEExprOpts_EqSingleDouble | OEExprOpts_EqAromatic, ringchecking = 'Strict', approximate = True  ):
    """Take two OEMol molecules, find the maximum common substructure and return it, along with the number of common heavy atoms. Minimum atoms constituting a match is 4 by default.
    
    Written by D. Mobley 1/27/11.

   Copyright (c) 2012, the University of New Orleans, written by David L. Mobley and Shuai Liu
   All rights reserved.
   Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

        * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
        * Neither the name of the University of New Orleans nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    
    The specific search done here is optimized around bonding patterns; atom matching is extremely loose. This gives great weight to finding rings with matching numbers of atoms, even if the composition of rings is not identical.
    
    The maximal substructure returned contains only WHOLE rings -- first the Maximal Common Substructure is calculated and then any partial rings created by deletion of atoms are removed, leaving only whole rings plus non-ring components.
    
    Arguments:
     - ligand 1: OEMol of molecule 1
     - ligand 2: OEMol of molecule 2
     Optional:
     - debug: False by default; if true, prints extra debugging info
     - min_atoms: 4 by default; if no match is found meeting this, the number of atoms returned will be zero and the match will be None.
     - Approximate: Default True; specify False to override. Switches between OEMCS types -- approximate or not. Approximate is a lot faster (and in general probably better for large sets) but for some small molecules can occasionally lead to problems. 
     - atomexpr:  OEExprOpts_EqONS by default, which gives a relatively strict matching when used in combination with the default bondexpr. Option 2 (use with option 2 bond expression) for a much less strict matching:  OEExprOpts_RingMember #Check only ring membership -- atoms match if they share common ring membership 
     - bondexpr:  OEExprOpts_BondOrder | OEExprOpts_EqSingleDouble | OEExprOpts_EqAromatic by default, which gives a relatively strict matching. Option 2 for a much less strict matching: OEExprOpts_ExactBonds | OEExprOpts_EqAromatic | OEExprOpts_EqSingleDouble #CHeck that bonding pattern matches -- either exactly, or counting single and double bonds as the same, or counting aromatic/nonaromatic as the same.
     - ringchecking: 'Strict' or 'Unstrict': Strict checking enforces the requirement that any remaining ring in the match must have a matching ring system of equivalent size in both ligands. This means that no deletion of a single ring in a joined/multiple ring system is possible. Unstrict checking allows single rings to be left behind from multiple ring systems, as long as these rings are aromatic (the idea being that these will be rigid and unlikely to be affected much by the conformation of any connected components which were deleted). Default: Strict checking.

     NOTE on atom and bond expressions: A less strict matching does not necessarily result in more matching atoms, because of conditions on the search. 


     Returns:
     - substructure: maximum substructure found
     - numatoms: Number of heavy atoms in match
    """

    #Make a log file for errors/warnings
    #oeerrs = oeofstream('mcss_log.txt')
    #OEThrow.SetOutputStream(oeerrs)

    #DLM 12/8/11 -- ethane is not being recognized as a match to itself since it has fewer than four heavy atoms. Attempt to work around by recognizing when the molecule size is smaller than the minimum number of match atoms and adjusting that threshold downward if necessary
    hvy1 = OECount( ligand1, OEIsHeavy() )
    hvy2 = OECount( ligand2, OEIsHeavy() )
    if hvy1 < min_atoms or hvy2 < min_atoms:
        min_atoms = min(hvy1, hvy2)

    #INTERNAL METHODS
    def delete_unconnected_components( oemol, debug = False):
        """Take an oemol with unconnected components due to deletion of atom(s) or bond(s). Identify the largest remaining component, and delete the rest."""

        #Identify the parts of the molecule    
        count, parts = OEDetermineComponents( oemol) #Count is the number of parts. 'parts' is a list specifying the part index, by atom number
        #if debug: print "Starting part number %s..." % count
        partnums = range( min(parts), max(parts)+1 )
        #Find how many atoms in each part
        atoms_per_part = [ parts.count( partnum ) for partnum in partnums ]
        maxpartnum = max( atoms_per_part )
        keepgrp = atoms_per_part.index( maxpartnum)
        keeppart = partnums[ keepgrp ]
        
        if debug:
            print "Molecule initially had %s parts. Retaining part number %s with %s atoms. Total atom number %s..." % (count, keeppart, maxpartnum, m.NumAtoms() )
        
        #Now delete all the atoms in other parts
        for atom in oemol.GetAtoms():
            partnum = parts[ atom.GetIdx() ]
            if partnum <> keeppart:
                oemol.DeleteAtom( atom )

        return OEMol( oemol )

    #To avoid some pathologies for ligands with multiple interconnected rings, track WHICH ring each atom is in. We want to only delete entire rings. 4/6/11
    def get_ring_indices_and_sizes( oemol, debug = debug):
        """Take an OEMol. Return two dictionaries -- one giving ring indices keyed by atom name, and another giving ring sizes keyed by atom name, where the ring size is the size of the "Ring System" (OE definition) they are in."""

        ring_indices = {}
        OEFindRingAtomsAndBonds( oemol )
        nrrings, rings = OEDetermineRingSystems( oemol )
        if debug: print "Molecule has %s rings..." % nrrings

        atomnames = []
        for (idx, atom) in enumerate( oemol.GetAtoms() ):
            atomname = atom.GetName()
            if rings[ idx ] > 0:
                ring_indices[ atomname ] = rings[idx]
            atomnames.append(atomname)

        #Now count how many atoms in each ring index
        ring_sizes = {}
        for ringindex in ring_indices.values():
            for (idx, atom) in enumerate( oemol.GetAtoms() ):
                atomname = atom.GetName()
                if not ring_sizes.has_key( atomname):
                    ring_sizes[atomname] = 0
                if rings[idx] == ringindex:
                    ring_sizes[atomname] +=1
        
        #Fill in ring indices and ring sizes for atoms which are not in rings
        for atom in atomnames:
            if not ring_indices.has_key( atom ):
                ring_indices[ atom ] = None
                ring_sizes[ atom ] = 0

        return ring_indices, ring_sizes



    #START OUR WORK
    #Initialize common substructure with first ligand
    common_substructure = ligand1.CreateCopy()
    ref = ligand2.CreateCopy()

    OEPerceiveChiral( common_substructure ) #Chirality is not automatically perceived when reading from a file
    OEPerceiveChiral( ligand2 ) #Chirality is not automatically perceived when reading from a file

    #Debug: Show original atoms
    if debug:
        atom_index = 1
        print "Initial ligand for the match:"
        for atom in common_substructure.GetAtoms():
            print "%5d %6s %12s" % (atom_index, atom.GetName(), atom.GetType())
            atom_index += 1
        print ""

    #Compute ring sizes and indices
    ring_indices, ring_sizes = get_ring_indices_and_sizes( common_substructure )
    #Compute them for the target too
    lig2_ring_indices, lig2_ring_sizes = get_ring_indices_and_sizes( ref )


    #Create mcss search from ligand 2

    #Configure search
    if not approximate:
        mcss = OEMCSSearch( ligand2, atomexpr, bondexpr )
    else:
        mcss = OEMCSSearch( ligand2, atomexpr, bondexpr, OEMCSType_Approximate ) #DLM 4/6/11 trying this for speed increase

    #Set minimum atoms
    mcss.SetMinAtoms(min_atoms)
    #Set max matches
    mcss.SetMaxMatches( 5 ) #DLM 4/6/11 modified to 5

    #mcss.SetMCSFunc( OEMCSMaxBondsCompleteCycles() ) #Modify scoring function to prefer completing cycles
    mcss.SetMCSFunc( OEMCSMaxAtomsCompleteCycles() ) #Modify scoring function to prefer completing cycles -- DLM 4/6/11 swapped to this

    if debug: print "Attempting to match %s with %s..." % (common_substructure.GetTitle(), ligand2.GetTitle() )

    nmatched = 0
    #Perform match and store atom mappings; also track names
    mappings = {}
    name_map = {}
    map_chiralities = {} #Store chiralities in the map (not original) for reference later
    map_aromaticities = {} #Store chiralities in the map (not original) for reference later
    map_isHeavy = {} #Dictionary storing whether BOTH of the atoms corresponding to the substructure are heavy or not.
    for match in mcss.Match(common_substructure, True): #Unique search for common substructure

        #Loop over matched atoms and store atom mappings for later use in checking ring membership
        for matchpair in match.GetAtoms():
            mappings[ matchpair.target ] = matchpair.pattern
            name_map[ matchpair.target.GetName() ] = matchpair.pattern.GetName()
            if debug: print "Proposed map %s->%s" % ( matchpair.target.GetName(), matchpair.pattern.GetName() )
            map_chiralities[ matchpair.target.GetName() ] = matchpair.pattern.IsChiral() 
            map_aromaticities[ matchpair.target.GetName() ] = matchpair.pattern.IsAromatic()
            #DLM 12/8/11 adding tracking of heavy atoms...
            map_isHeavy[ matchpair.target.GetName()] = True
            if matchpair.target.IsHydrogen() or matchpair.pattern.IsHydrogen():
                map_isHeavy[ matchpair.target.GetName() ] = False

        #Create a match substructure
        m = OEMol()
        OESubsetMol(m, match, True) #Get common substructure
    


        m_ring_indices, m_ring_sizes = get_ring_indices_and_sizes( m , debug = debug)

        #MCSS search as currently implemented ignores chirality (works with non-stereo smiles, it appears). To deal with this, break the molecule at any chiral centers and keep only the largest piece as our substructure. This is drastic (and we should come up with a better fix eventually) but if we don't do this, we will end up with a 'common substructure' that has different chirality than one of the molecules being compared. (ex. 326_0 vs 326_1).
        OEPerceiveChiral( m ) #Chirality is not automatically perceived when reading from a file
        for atom in m.GetAtoms():
            if atom.IsChiral() or map_chiralities[ atom.GetName() ]: 
                if debug: print "Chiral atom %s" % atom.GetName()
                #Deleting all chiral atoms is actually too restrictive as it will result in deletion of whole rings because this will often delete an atom in a ring. 
                #Better: Two scenarios: (1) Chiral atom is in a ring, in which case delete any non-ring atoms it's connected to, as well as any atoms in different neighboring rings (don't want to retain chiral connections between rings). (2) Chiral atom is not in a ring, in which case delete it.
                if not atom.IsInRing(): #If it is not in a ring, delete it
                    m.DeleteAtom(atom)
                    if debug: print "Deleting chiral atom %s as it is not in a ring..." % atom.GetName()
                #If it is in a ring, delete any connected atoms which are not in the same ring
                else:
                    thisring = m_ring_indices[ atom.GetName() ]
                    for nbr in atom.GetAtoms():
                        nbrname = nbr.GetName()
                        if m_ring_indices.has_key( nbrname): nbrring = m_ring_indices[ nbrname]
                        else: nbrring = None
                        if thisring <> nbrring:
                            m.DeleteAtom( nbr )
                            if debug: print "Deleting atom %s because it adjoins a chiral atom in a ring and is not in the same ring..." % nbr.GetName()

        m = delete_unconnected_components( m, debug = debug)

        #Recompute ring indices and sizes
        m_ring_indices, m_ring_sizes = get_ring_indices_and_sizes( m )

        if ringchecking=='Strict':
            #STRICT RING MEMBERSHIP CHECKING
            #Delete all atoms which either (a) underwent a change in ring size going to this molecule, or (b) are in a different ring size in the other molecule (DLM 4/6/11)
            for atom in m.GetAtoms():
                atomname = atom.GetName()
                matchname = name_map[ atomname ]
                #Case a: Different ring size here
                if m_ring_sizes[atomname] <> ring_sizes[ atomname]:
                    m.DeleteAtom(atom)
                    if debug: print "Deleting atom %s due to change of ring size..." % atom.GetName()
                elif m_ring_sizes[atomname] <> lig2_ring_sizes[ matchname]:
                    m.DeleteAtom(atom)
                    if debug: print "Deleting atom %s due to change of ring size relative to other ligand..." % atom.GetName()

        elif ringchecking == 'Unstrict':
            #NONSTRICT RING MEMBERSHIP CHECKING
            #Delete all atoms which either (a) underwent a change in ring membership going to this molecule, or (b) have different ring membership in the other molecule. ALSO, to prevent leaving nonplanar rings which would have had conformations affected by adjoining rings that are being deleted, (c) delete any non-aromatic ring (or ring which was non-aromatic in the other molecule) which underwent a change in ring size going to this molecule or (d) are in a different ring size in the other molecule
            for atom in m.GetAtoms():
                atomname = atom.GetName()
                matchname = name_map[ atomname]
                #Case (a) -- change in ring membership
                if m_ring_indices[atomname]==None and ring_indices[ atomname]<>None:
                    if debug: print "Atom %s was in a ring but is not now; deleting." % atom.GetName()
                    m.DeleteAtom(atom)
                #Case (b) -- change in ring membership relative to match
                elif m_ring_indices[atomname]==None and lig2_ring_indices[ matchname]<>None:
                    if debug: print "Atom %s is in a ring in the other pair member but is not here; deleting." % atom.GetName()
                    m.DeleteAtom(atom)
                #Case (c-d) -- atom is in a ring but not an aromatic one, or one that is not aromatic in the other molecule
                elif atom.IsInRing() and (not atom.IsAromatic() or not map_aromaticities[ atomname]):
                    #Case (c) -- if it underwent a change in ring size getting here
                    if m_ring_sizes[atomname] <> ring_sizes[atomname]:
                        if debug: print "Atom %s is in a non-aromatic ring that underwent a change of size. Deleting." % atom.GetName()
                        m.DeleteAtom(atom)
                    elif m_ring_sizes[atomname] <> lig2_ring_sizes[matchname]:
                        if debug: print "Atom %s is in a non-aromatic ring that had a different ring size in the other molecule. Deleting." % atom.GetName()
                        m.DeleteAtom(atom)

        else:
            print "Invalid ringchecking argument."
            raw_input()

        #THIS CURRENTLY LEAVES HANGING ATOMS IN MANY CASES; THESE NEED TO BE DELETED -- do the deletion
        if m.NumAtoms() >0: m =  delete_unconnected_components( m, debug = debug ) #DLM 4/6/11 swapped to delete_unconnected_components

        #DLM 4/6/11 updating to count number of heavy atoms instead of number of atoms
        #DLM 12/8/11 updating this to use preassigned isHeavy flags to avoid counting cases where it is a heavy atom in the scaffold but a hydrogen in the pattern
        nmatched = 0
        for atom in m.GetAtoms():
            if map_isHeavy[ atom.GetName() ]:
               nmatched+=1 



        if debug:
            #Print all retained atoms
            print "Atoms retained after deletion of connected atoms: "
            atom_index = 1
            
            print "%5s %6s %12s %6s" % ("Index", "Name", "Type", "Name in target molecule")
            for atom in  m.GetAtoms():
                thisname = atom.GetName()
                if name_map.has_key( thisname ):
                    origname = name_map[ thisname ]    
                else: origname = 'none'
                print "%5d %6s %12s %6s" % (atom_index, atom.GetName(), atom.GetType(), origname)
                atom_index += 1
            print ""

        #Only need to consider first match if it has any matching atoms:
        if nmatched >= min_atoms:
            break
            print "Exiting due to finding first match."

    #Upon deletion of some hydrogen atoms, they get transferred to being implicit hydrogen atoms. At this point we should have NO implicit hydrogen atoms -- remove them all.
    if nmatched >= min_atoms:
        for atom in m.GetAtoms():
            if atom.GetImplicitHCount() > 0:
                atom.SetImplicitHCount(0)

    #Ensure the match we are going to return has enough atoms
    if nmatched >= min_atoms: #DLM 12/8/11 switched this and the two immediately preceding to >= from >    
        return m, nmatched
    else:
        return None, 0

