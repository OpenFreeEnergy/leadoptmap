#! usr/bin/python

from pygraph.classes.graph import graph
from numpy import *  # for unique
import copy
import pickle

def unique_list(list):
    new_list = []
    mult_list = []
    for member in list:
        num = list.count(member)
        if num == 1:
            new_list.append(member)
        else:
            mult_list.append(member)

    mult_list.sort()
    if len(mult_list) != 0:
        new_list.append(mult_list[-1])
    return new_list



def find_molecule_not_in_circle (clus, gr, debug = False):
    clus_in_function = copy.copy(clus)
    if len(clus_in_function) < 3:
        molecule_not_in_circle = []
    else:
        molecule_not_in_circle = []
        list_level_dic = {}
        dead_son_dic = {}
        max_num = 0
        list_level_dic[1] = []
        list_total = clus_in_function
        lable_mol = {}
        dead_mol = {}
#find the one with most neighbors and put it in to first level
        for member in clus_in_function:
            if len(gr.node_neighbors[member]) >= max_num:
                max_num = len(gr.node_neighbors[member])
                max_neighbor = member
        list_level_dic[1].append(max_neighbor)
        level_num = 1
        level_factor = max_num + 1
        lable_mol[0] = max_neighbor

        count_loop = 0
#loop to find not in circle one level by level
        while True:
            available_son = []
            count_loop += 1
            for father in list_level_dic[level_num]:
                j = 1
                for k in lable_mol.keys():
                    if lable_mol[k] == father:
                        lable_of_member = k
                for son in gr.node_neighbors[father]:
                    if son in list_total and son != lable_mol[(int (lable_of_member/level_factor))]:
                        if son not in available_son:
                            available_son.append(son)
                            lable_mol[ lable_of_member * level_factor + j] = son
                            j += 1
            if len(available_son) == 0:
#########sliu 1/26/12  the last delete
                fake_dead = []
                def_dead_son_lable = []
                def_gen_son_lable = []                                                                     
                list_total_before = copy.copy(list_total)
                for last_member in list_total_before:                                                           
                    for l in lable_mol.keys():                                                             
                        if lable_mol[l] == last_member:
                            lable_of_last_member = l
                    for last_son in gr.node_neighbors[last_member]:                                       
                        if last_son not in list_total and last_son != lable_mol[(int (lable_of_last_member/level_factor))]:
                            fake_dead.append(last_son)                                                          
                            if last_member not in fake_dead:                                                    
                                fake_dead.append(last_member)                                                   
                    if len(fake_dead) != 0:                                                                     
                        for def_dead in fake_dead:                                                         
                            for m in lable_mol.keys():                                                     
                                if lable_mol[m] == def_dead:                                               
                                    def_dead_son_lable.append(m)                                           
                                                                                                           
                        for n in def_dead_son_lable:                                                       
                            while True:                                                                    
                                if n in lable_mol.keys():                                                  
                                    def_gen_son_lable.append(n)                                            
                                n = int(n/level_factor)                                                    
                                if n == 0: 
                                    def_gen_son_lable.append(n)                                                                
                                    break                                                                  
                        def_gen_son_lable = unique_list(def_gen_son_lable)                                 
                        for p in def_gen_son_lable:                                                        
                            if lable_mol[p] in list_total:                                                 
                                list_total.remove(lable_mol[p])
                                print "Found last Delete"                                            
                    fake_dead = []                                                                         
                    def_dead_son_lable = []                                                                
                    def_gen_son_lable = []                                                                 
###############################
                break

##############sliu 3/21/12
            available_son_stable = copy.copy(available_son)
            for lucky_son in available_son_stable:
##########
 
                dead_son_lable = []
                dead_generation_lable = []
                for q in lable_mol.keys():
                    if lable_mol[q] == lucky_son:
                        dead_son_lable.append(q)
               
                if len(dead_son_lable) == 1:
                    dead_son_lable = []
                else:
                    for r in dead_son_lable:
                        if lable_mol[r] in available_son:
                            available_son.remove(lable_mol[r])
                        while True:
                            if r in lable_mol.keys():
                                dead_generation_lable.append(r)
                            r = int(r/level_factor)
                            if r == 0:
                                dead_generation_lable.append(r)
                                break  
                    dead_generation_lable = unique_list(dead_generation_lable)
                    for t in dead_generation_lable:
                        if lable_mol[t] in list_total:
                            list_total.remove(lable_mol[t])
                    for lables in dead_generation_lable:
                        dead_mol[lables] = lable_mol[lables]

                if len(available_son) == 0:
                    break
            good_son = available_son
            level_num += 1
            list_level_dic[level_num] = good_son
            available_son = []
        if len(list_total) != 0:
    
            molecule_not_in_circle = list_total
            list_total = []
        lable_mol = {}
    return molecule_not_in_circle


