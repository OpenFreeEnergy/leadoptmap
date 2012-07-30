"""Defines similarity between two molecules
"""



import math



# Similarity definition:
# - Floating number
# - Range: [0, 1]. The values 0 and 1 represent the two extremes: alienity and identity, respectively.



def by_atom_count( mol0, mol1, mcs ) :
    """
    
    """
    delta_total_num_atom = len( mol0.atom          ) + len( mol1.atom          ) - 2 * len( mcs.atom          )
    delta_num_heavy_atom = len( mol0.heavy_atoms() ) + len( mol1.heavy_atoms() ) - 2 * len( mcs.heavy_atoms() )
    delta_num_hydrogen   = delta_total_num_atom - delta_num_heavy_atom
    BETA                 = -0.1
    energy               = delta_num_heavy_atom + 0.25 * delta_num_hydrogen
    return math.exp( BETA * energy )



def by_heavy_atom_count( mol0, mol1, mcs ) :
    """
    
    """
    delta_num_heavy_atom = len( mol0.heavy_atoms() ) + len( mol1.heavy_atoms() ) - 2 * len( mcs.heavy_atoms() )
    BETA                 = -0.1
    energy               = delta_num_heavy_atom
    return math.exp( BETA * energy )
    

