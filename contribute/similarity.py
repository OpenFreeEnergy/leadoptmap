"""Defines similarity scoring functions
"""



import math



# Think of the similarity score as a probability value, so it has the following characteristics:
# - Positive floating number
# - Range: [0, 1]. The values 0 and 1 represent the two extremes: alienity and identity, respectively.



def by_atom_count( mol0, mol1, mcs ) :
    """
    This function determines the similarity score by counting the number of extra atoms between two molecules. In the scenario
    of maximum common substructure, we compare each of the two molecule with the common substructure and then add up the two
    results. The final sum is considered as the difference in terms of atom counts between the two molecules.
    We tentatively design the scoring function to use the following formula:
        exp( -BETA * E )
    where BETA is a constant, which we give an arbitrary value 0.4, and E is a weighted atom count tentatively given by the
    formula:
        E = num-of-heavy-atoms + 0.25 * num-of-hydrogen-atoms

    @type  mol0: C{Struc}
    @param mol0: First molecule
    @type  mol1: C{Struc}
    @param mol1: Second molecule
    @type  mcs : C{Struc}
    @param mcs : Common substructure between C{mol0} and C{mol1}
    """
    delta_total_num_atom = len( mol0.atom          ) + len( mol1.atom          ) - 2 * len( mcs.atom          )
    delta_num_heavy_atom = len( mol0.heavy_atoms() ) + len( mol1.heavy_atoms() ) - 2 * len( mcs.heavy_atoms() )
    delta_num_hydrogen   = delta_total_num_atom - delta_num_heavy_atom
    BETA                 = 0.4
    energy               = delta_num_heavy_atom + 0.25 * delta_num_hydrogen
    return math.exp( -BETA * energy )



def by_heavy_atom_count( mol0, mol1, mcs ) :
    """
    Similar to C{by_atom_count} except that we only consider heavy atoms here
    """
    delta_num_heavy_atom = len( mol0.heavy_atoms() ) + len( mol1.heavy_atoms() ) - 2 * len( mcs.heavy_atoms() )
    BETA                 = 0.4
    energy               = delta_num_heavy_atom
    return math.exp( -BETA * energy )
    

