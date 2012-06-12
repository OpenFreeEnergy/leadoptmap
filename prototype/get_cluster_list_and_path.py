#!usr/bin/python

######### take clusters.pickle and titles_and_mcssscorearray_charges.pickle to get a clusters_by_path.pickle. This pickle files contain the path of the molecule based on the clusters they belong to. This pickle file make it  easy to check the cluster by pymol and get clustrer_list for next calculation planning strp########


import glob
import pickle


filenames = glob.glob('mol2_file/*.mol2')

correspond_dir = {}

for file in filenames:  

    f = open(file,'r')
    lines = f.readlines()
#    print lines[1] 
    target_line = lines[1].split()[0] 
    #correspond_dir[file] = lines[1]
    titles = file.split()[-1]    
    correspond_dir[titles] = target_line 

#print correspond_dir 
#print filename

#print correspond_dir

title_file = open ( 'titles_and_mcssscorearray_charges.pickle' , 'r')
(titles, scores, natoms) = pickle.load(title_file)

cluster_file = open('clusters.pickle','r')
(clusters) = pickle.load(cluster_file)

clusters_by_title = {}

for i in range (len(titles)):
    for j in clusters.keys():
        if i in clusters[j]:
            if not clusters_by_title.has_key(j):
                clusters_by_title[j] = []
            clusters_by_title[j].append(titles[i])

print clusters_by_title


clusters_by_path = {}

for i in clusters_by_title.keys():
    for name in (clusters_by_title[i]):
        for dir in correspond_dir.keys():
            if correspond_dir[dir] == name:
                if not clusters_by_path.has_key(i):
                    clusters_by_path[i] = []
                clusters_by_path[i].append(dir)

cluster_list = []

for i in clusters_by_title.keys():
    cluster_content = clusters_by_title[i]
    cluster_list.append(cluster_content)

#print clusters_by_path
    
f = open('clusters_by_path.pickle','w')
pickle.dump(clusters_by_path,f)
f.close

f2 = open('clusters_by_title.pickle','w')
pickle.dump(clusters_by_title, f2)
f2.close

f3 = open ('cluster_list.pickle','w')
pickle.dump(cluster_list, f3)
f3.close()
