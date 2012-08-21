#!/usr/bin/env python

"""Attempt to to build a score matrix of molecule similarities based on the number of atoms in their maximal common substructures."""

#IDEA: Require good matches to have at least one ring in common; the more atoms they have in addition to the ring the better. NOTE that it is not necessary that all of the atoms in the ring be of the same type, but there must be the same number.

from mmtools.moltools.ligandtools import *
from mmtools.moltools.relativefeptools import *
from openeye.oechem import *
import numpy
import pickle
import glob
from MCSS_tool import *

#inputfile = '../molecules_charged.mol2'
#ifs = oemolistream( inputfile )

inputmolfiles = glob.glob('mol2_file/*.mol2') 

#Make a log file for errors/warnings
oeerrs = oeofstream('mcss_log.txt')
OEThrow.SetOutputStream(oeerrs)

#Read molecules into a dictionary by title; also store charges
molecules = {}
charges = {}
#for mol1 in ifs.GetOEMols():
for molfile in inputmolfiles:
    mol1 = readMolecule( molfile )
    title = mol1.GetTitle()
    newmol = OEMol( mol1 )
    molecules[ title ] = newmol

    #Make alphabetical list of titles
    titles = molecules.keys()
    titles.sort()

    #Obtain and store charge
    chg = OENetCharge( newmol )
    charges[ title ] = chg


#Create storage for a matrix of scores
Nmols = len(titles)
mcss_scores = numpy.zeros( (Nmols, Nmols), float) #NxN array initialized with zeros, where N is the number of molecules




def match_pair( name1, name2, inputdir = 'molecules_aligned', outputdir = 'mcss', debug = True, writePairs = True):
    """Use MCSS search to match a pair of molecules whose filenames are provided, determining the maximal common substructure and writing it out to specified output directory with a filename composed of the titles of the two molecules joined by an underscore."""

    mol1 = readMolecule( os.path.join( inputdir, name1) )
    mol2 = readMolecule( os.path.join( inputdir, name2) )

    if not os.path.isdir(outputdir): os.mkdir( outputdir )

    #Determine title of resulting mcss overlay
    joint_title = mol1.GetTitle() +'_'+ mol2.GetTitle()
    joint_title=joint_title.replace(' ','_')

    #Determine maximal common substructure
    mcss, numatoms = determineMCSS_requireRings(mol1, mol2, debug = False)

    if mcss and writePairs: #If we get a match:
        #Write out mcss
        ostream = oemolostream( os.path.join( outputdir, joint_title+'.mol2'))
        OEWriteMol2File( ostream, mcss)
        ostream.close()

    if debug: print "Number of matching atoms for MCSS of %s-%s is %s..." % ( mol1.GetTitle(), mol2.GetTitle(), numatoms)

    return mcss, numatoms


def get_ring_number( mol ):
    """Attempt to calculate and return the number of rings in a molecule."""
    num_components, component_membership = OEDetermineComponents( mol )
    num_rings = mol.NumBonds() - mol.NumAtoms() + num_components

    return num_rings

number_of_atoms = []


#DEBUGGING
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.155'], molecules['frag.vs.277_0'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.013_0'], molecules['frag.vs.068_0'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.326_0'], molecules['frag.vs.326_1'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.086'], molecules['frag.vs.133'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.287_0'], molecules['frag.vs.326_0'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.013_0'], molecules['frag.vs.013_1'], debug = True)
#Above are problem cases I printed out then fixed
#The immediately following is a problem case which is now fixed
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.039_1'], molecules['frag.vs.323_0'], debug = True)
#Other testing
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.155'], molecules['frag.vs.156'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['frag.vs.147'], molecules['frag.vs.383'], debug = True)
#(mcss, natoms) = determineMCSS_requireRings(molecules['CRA-18305'], molecules['CRA-19858'], debug = True)
#print "Number of atoms in common substructure %s..." % natoms#

#print "Pausing before entering main loop..."
#raw_input()
#END DEBUGGING


#Loop and do MCSS comparison of all by all 
for (idx, title1) in enumerate(titles):

    mol1 = molecules[ title1 ]
    print "\nFor mol1=%s (%s/%s):" % (title1, idx, Nmols)

    mol1_natoms = mol1.NumAtoms()
    number_of_atoms.append( mol1_natoms )

    #Get ready to do a progress bar
    dots = OEDots( 50, 1, "  molecules of %s" % len( titles[idx:]) )
    
    #Charge for this molecule
    chg1 = charges[ title1 ]

    for title2 in titles[idx:]:
        mol2 = molecules[ title2 ]
        mol2_natoms = mol2.NumAtoms()
        #Charge for thsi molecule
        chg2 = charges[title2]

        dots.Update()

        #Do mcss comparison -- only if charges for both molecules are the same, otherwise we don't need to compare (since we only are going to cluster molecules of the same net charge)
        if chg1 == chg2:
            #Initially apply with default options, giving a relatively strict matching
            (mcss, natoms) = determineMCSS_requireRings(mol1, mol2, debug = False, ringchecking = 'Unstrict')
            #Apply with looser options if no match found
            if natoms==0:
                (mcss, natoms) = determineMCSS_requireRings( mol1, mol2, debug = False, atomexpr = OEExprOpts_RingMember, bondexpr = OEExprOpts_ExactBonds | OEExprOpts_EqAromatic | OEExprOpts_EqSingleDouble, ringchecking = 'Unstrict' )
        else:
            (mcss, natoms) = None, 0

        #Write out structure
        #Determine title of resulting mcss overlay
        if mcss: #Only if there is actually a match!
            joint_title = title1 + "_" + title2
            joint_title=joint_title.replace(' ','_')
            ostream = oemolostream( os.path.join( 'mcss_unstrict', joint_title+'.mol2'))
            OEWriteMol2File( ostream, mcss)
            ostream.close()

        #Save number of common atoms
        idx2 = titles.index( title2 )
        mcss_scores[ idx, idx2 ] = natoms
        mcss_scores[ idx2, idx ] = natoms

#Store scores
#file = open( 'titles_and_mcssscorearray.pickle', 'w')
file = open( 'titles_and_mcssscorearray_charges_unstrict.pickle', 'w')
pickle.dump( (titles, mcss_scores, number_of_atoms), file)
file.close()

