import mcs
import struc
def matrix ( mols, mcs_ids, rule ):
    import numpy
    id_list = []
    id_vs_simi = {}
    id_vs_title = {}
    title_vs_simi = {}
    title_list = []
    for mol in mols:
        title = mol.title()
        id = mol.id()
        if id not in id_list:
            id_list.append(id)
        if title not in title_list:
            title_list.append(title)
        id_vs_title [id] = title
    for id in mcs_ids:
        id0, id1 = mcs.get_parent_ids(id)
        #if not id0 in id_list:
        #    id_list.append(id0)
        #if not id1 in id_list:
        #    id_list.append(id1)
        simi       = rule.similarity( id0, id1, mcs_id = id )        
        title0 = id_vs_title[id0]
        title1 = id_vs_title[id1]
        #id_vs_simi [(id0,id1)] = simi
        title_vs_simi [(title0,title1)] = simi
    #size          = len( id_list )
    size          = len( title_list )
    scores = numpy.zeros( (size, size,) )
    for i in range( size ) :
            scores[i, i] = 1.0
            for j in range( i + 1, size ) :
                #id_i = id_list[i]
                title_i = title_list[i]
                title_j = title_list[j]
                #simi = id_vs_simi[id_i,id_j]
                simi = title_vs_simi[(title_i,title_j)]
                scores[i, j] = simi                
                scores[j, i] = simi

    return (title_list, id_list, scores)                
