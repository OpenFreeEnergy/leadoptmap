import networkx as nx
from pygraph.classes.graph import graph
import pickle


def check_graph (net_file, py_file):
    """Take one networkx graph pickle file and one pygraph graph pickle file and get back similarity of these two graphs 
"""
    f1 = open( net_file,'r')
    net_graph = pickle.load(f1)
    f1.close()

    f2 = open(py_file ,'r')
    py_graph = pickle.load(f2)
    f2.close()

#####For the smae graph, edges in pygraph are ('mol1','mol2'), ('mol2', 'mol1'). While networkx graph only have one edge ('mol1','mol2'), or ('mol2', 'mol1'). So here expand networkx graph to compare.
    non_direct_net_edges = []
    ori_net_edges = []
    for net_edges in net_graph.edges():
        non_direct_net_edges.append(net_edges)
        ori_net_edges.append(net_edges)
        non_direct_net_edges.append((net_edges[1],net_edges[0]))
    
    py_edges = py_graph.edges()

    same_edges = []
    py_only_edges = []
    net_only_edges = []

    for edge in py_edges:
        if edge in non_direct_net_edges:
#        if edge in ori_net_edges:
            same_edges.append(edge)
        else: 
            py_only_edges.append(edge)
    
    for next_edge in non_direct_net_edges:
#    for next_edge in ori_net_edges:
        if next_edge not in same_edges:
            net_only_edges.append(next_edge)
    return(same_edges, py_only_edges, net_only_edges, len(py_edges), len(non_direct_net_edges))
#    return(same_edges, py_only_edges, net_only_edges, len(py_edges), len(ori_net_edges))

#(common, py_only, net_only, len_py, len_net) = check_graph('after_step_one.pickle' , './lee_pickle/after_step_one.pickle' ) 
#(common, py_only, net_only, len_py, len_net) = check_graph('after_step_two.pickle' , './lee_pickle/after_step_two.pickle' ) 
#(common, py_only, net_only, len_py, len_net) = check_graph('after_step_three.pickle' , './lee_pickle/after_step_three.pickle' )  
#(common, py_only, net_only, len_py, len_net) = check_graph('../use_py_to_test_net_trypsin/after_step_three.pickle' , '../trypsin_skip/after_step_three.pickle' )  
#(common, py_only, net_only, len_py, len_net) = check_graph('../use_py_to_test_net_trypsin/graph_new_fuction.pickle' , '../trypsin_skip/graph_new_fuction.pickle' )  
#(common, py_only, net_only, len_py, len_net) = check_graph('graph_new_fuction.pickle' , './lee_pickle/graph_new_fuction.pickle' )  
(common, py_only, net_only, len_py, len_net) = check_graph('final_graph.pickle' , './lee_pickle/final_graph.pickle' ) 

print "Common edges # %s \n"% len(common) , common
print "Pygraph only edges # %s \n" % len(py_only) , py_only
print "Networkx only edges # %s \n" % len(net_only) , net_only
print "num of edges of pygraph %s \n" % len_py
print "num of edges of networkx %s \n" % len_net
    

    
