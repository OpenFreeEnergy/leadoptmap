#usr/bin/python

import pickle

def get_pair_MCSS_score(score_file = 'titles_and_mcssscorearray_charges.pickle', mol1 = 'APC-1-762', mol2 = 'APC-1144'):
"""
use score matrix, return a pair of MCSS score between mol1 and mol2. 
"""
    file1 = open( score_file,'r')
    (titles,scores,natoms) = pickle.load(file1)
    pos1 = 0
    for item1 in titles:
        if item1 == mol1:
            break
        else:
            pos1 += 1
    pos2 = 0
    for item2 in titles:
        if item2 == mol2:
            break
        else:
            pos2 += 1
    print pos1, pos2

print get_pair_MCSS_score (score_file = 'titles_and_mcssscorearray_charges.pickle', mol1 = 'APC-1-762', mol2 = 'APC-1144')
    

