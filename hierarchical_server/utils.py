# -*- coding: utf-8 -*-
import os
import numpy as np

def convert2List(pnp):
    temp = list(pnp)
    plist = []
    for i in range(len(temp)):
        plist.append(list(temp[i])) 

    return plist

def parseOTUOTU(otu_file, top_ratio):
    otu_otu = np.loadtxt(otu_file)
    otu_otu = np.tril(otu_otu, -1)
    ratio_value = np.percentile(np.abs(otu_otu[otu_otu != 0]), 100 - float(top_ratio))
    otu1 = []
    otu2 = []
    weight = []
    p = otu_otu.shape[0]
    for i in range(p):
        for j in range(i):
            if np.abs(otu_otu[i,j]) >= ratio_value:
                otu1.append(i)
                otu2.append(j)
                weight.append(otu_otu[i,j])

    sort_asso = sorted(zip(otu1, otu2, weight), key = lambda x: abs(x[2]), reverse = True)
    return sort_asso

def parseEFOTU(ef_file, top_ratio):
    ef_otu = np.loadtxt(ef_file)
    ratio_value = np.percentile(np.abs(ef_otu), 100 - float(top_ratio))
    otu1 = []
    ef1 = []
    weight = []
    q, p = ef_otu.shape
    for i in range(q):
        for j in range(p):
            if np.abs(ef_otu[i,j]) >= ratio_value:
                otu1.append(j)
                ef1.append(i)
                weight.append(ef_otu[i,j])

    sort_asso = sorted(zip(ef1, otu1, weight), key = lambda x: abs(x[2]), reverse = True)
    return sort_asso

def readAssociation(stat_dir, asso_info, top_ratio = 10):
    ef_otu_file = stat_dir + "/ef_otu_association"
    otu_otu_file = stat_dir + "/otu_otu_association"
    
    asso_info["otu_otu"] = parseOTUOTU(otu_otu_file, top_ratio)
    asso_info["ef_otu"] = parseEFOTU(ef_otu_file, top_ratio)
    asso_info["otu_otu_num"] = len(asso_info["otu_otu"])
    asso_info["ef_otu_num"] = len(asso_info["ef_otu"])
    asso_info["ef_otu_ratio"] = top_ratio
    asso_info["otu_otu_ratio"] = top_ratio

def combineMeanStd(mean_list, std_list):
    res_list = []
    num = len(mean_list)
    for i in range(num):
        per_str = "{0:.2g} ({1:.2g})".format(mean_list[i], std_list[i])
        res_list.append(per_str)

    return res_list

def combineRatios(ratio0, ratio1):
    res_list = []
    num = len(ratio0)
    for i in range(num):
        per_str = "{0:.2f}% / {1:.2f}%".format(100 * ratio0[i], 100 * ratio1[i])
        res_list.append(per_str)

    return res_list

def readMeanInfo(dir_name, mean_info):
    attr_file = dir_name + "/attr_info"
    meta_file = dir_name + "/meta_info"
    otu_file = dir_name + "/otu_info"

    # read meta 
    meta_info = np.loadtxt(meta_file)
    meta_res = combineMeanStd(meta_info[0,:], meta_info[1,:])
    # read otu
    otu_info = np.loadtxt(otu_file)
    otu_res = combineMeanStd(otu_info[0,:], otu_info[1,:])
    # read attr
    attr_info = np.loadtxt(attr_file)
    attr_res = combineRatios(attr_info[0,:], attr_info[1,:])

    mean_info["meta"] = meta_res
    mean_info["otu"] = otu_res
    mean_info["attr"] = attr_res

def readList(file):
    fnames = []
    for perl in open(file).readlines():
        fnames.append(perl.strip())
    return fnames

def readMatrixShape(matrix_file):
    matrix_shape = []
    for perl in open(matrix_file).readlines():
        matrix_shape.append(int(perl.strip()))

    return matrix_shape

def filterOTUNames(otu_name_list):
    new_name = []
    for per_name in otu_name_list:
        per_splits = per_name.split(";")
        split_num = len(per_splits)
        need_num = 2
        per_new_name = ""
        i = split_num - 1
        while i > -1:
            if len(per_splits[i]) > 3:
                need_num -= 1
                per_new_name = per_splits[i] + "; "+ per_new_name
                if need_num == 0:
                    break

            i -= 1;

        new_name.append(per_new_name)

    return new_name

def readOTUAndEF(otu_name_file, meta_name_file, attr_name_file):
    print("read OTU name !")
    otu_name = readList(otu_name_file)
    print("filter OTU name !")
    otu_name_filter = filterOTUNames(otu_name)
    print("read Meta name !")
    ef_name = readList(meta_name_file)
    print("read Attr name !")
    attr_name = readList(attr_name_file)
    
    return otu_name_filter, ef_name, attr_name

def getClusterNumber(dir_name):
    index = dir_name.rfind("-")
    return int(dir_name[(index+1):])

def getMergedInfo(merged_dir):
    stat_dict = {}
    stat_dir = merged_dir + "/stat/"
    basic_info_file = stat_dir + "/basic_info"
    basic_info = []
    for per_l in open(basic_info_file, 'r').readlines():
        basic_info.append(int(per_l.strip()))

    stat_dict["N"] = basic_info[0]
    stat_dict["P"] = basic_info[1]

    stat_dict["dir"] = stat_dir

    return stat_dict

def parseTree(node_dict, node_name, node_dir, layer, cluster_number):
    if layer == 1:
      node_dict["name"] = node_name
    else:
      node_dict["name"] = "C" + str(cluster_number) + "-" + node_name

    merged_dir = node_dir + "/merged/"
    added_dir = node_dir + "/added/"
    original_dir = node_dir + "/original/"

    if os.path.exists(merged_dir):
        merged_info = getMergedInfo(merged_dir)
        node_dict["merged"] = merged_info 

    if os.path.exists(original_dir):
        node_dict["original"] = []
        original_node_dict = {}
        parseTree(original_node_dict, "original-" + str(layer), original_dir, layer + 1, cluster_number)
        node_dict["original"].append(original_node_dict)

    if os.path.exists(added_dir):
        node_dict["added"] = []
        added_node_dict = {}
        parseTree(added_node_dict, "added-" + str(layer), added_dir, layer + 1, cluster_number) 
        node_dict["added"].append(added_node_dict)

def getSampleIndexInfo(index_dir):
  index_file = index_dir + "/Sample_Index_1"
  sample_index = np.loadtxt(index_file)
  return sample_index.astype(np.integer)

def parseIndexTree(index_dict, node_dir, node_name, layer, cluster_number):
    if layer != 1:
      node_name = "C" + str(cluster_number) + "-" + node_name

    merged_dir = node_dir + "/merged/"
    added_dir = node_dir + "/added/"
    original_dir = node_dir + "/original/"

    if os.path.exists(merged_dir):
        index_info = getSampleIndexInfo(merged_dir)
        index_dict[node_name] = index_info

    if os.path.exists(original_dir):
        parseIndexTree(index_dict, original_dir, "original-" + str(layer), layer + 1, cluster_number)

    if os.path.exists(added_dir):
        parseIndexTree(index_dict, added_dir, "added-" + str(layer), layer + 1, cluster_number) 


def parseSampleIndex(tree_dir):
    index_dict = {}
    layer = 1
    first_layer_name = []
    for per_dir in os.listdir(tree_dir):
        per_full_dir = tree_dir + "/" + per_dir
        if os.path.isdir(per_full_dir):
            cluster_number = getClusterNumber(per_dir)
            node_name = "Cluster-" + str(cluster_number)
            parseIndexTree(index_dict, per_full_dir, node_name, layer, cluster_number)
    
    return index_dict


def parseHier(tree_dir, meta_pres, otu_pres, attr_pres, cluster_names):
    tree_dict = {}
    tree_dict["root"] = []
    tree_dict["first_layer"] = []
    k = 0
    pvalue_cnames = []
    for per_dir in os.listdir(tree_dir):
        per_dir_full = tree_dir + "/" + per_dir
        if os.path.isdir(per_dir_full):
            cluster_number = getClusterNumber(per_dir)
            k += 1
            per_cluster_dict = {}
            parseTree(per_cluster_dict, "Cluster-" + str(cluster_number), per_dir_full, 1, cluster_number)
            tree_dict["root"].append(per_cluster_dict)
            pvalue_cnames.append("Cluster-" + str(cluster_number))

            # read mean info
            per_stat_dir = per_dir_full + "/merged/stat/"
            per_mean_info = {}
            per_merged_info = getMergedInfo(per_dir_full + "/merged/")
            per_mean_info["cluster_name"] = "Cluster-" + str(cluster_number) + ", P={0}, N={1}".format(per_merged_info["P"], per_merged_info["N"])
            readMeanInfo(per_stat_dir, per_mean_info)
            per_mean_info["stat_dir"] = per_stat_dir
            # read assocation info
            per_asso_info = {}
            readAssociation(per_stat_dir, per_asso_info)
            per_mean_info["association"] = per_asso_info

            tree_dict["first_layer"].append(per_mean_info)

    tree_dict["K"] = k
    search_res = {}
    searchPvalues(pvalue_cnames, meta_pres, otu_pres, attr_pres, cluster_names, search_res)
    tree_dict["Pvalue"] = search_res

    return tree_dict

def mapCname2Index(cnames, cluster_names):
    # print(cnames)
    # print(cluster_names)
    index = []
    for i in range(len(cnames)):
        index.append(cluster_names.index(cnames[i]))
    index.sort()
    index = np.array(index).astype(np.integer)
    cluster_names = np.array(cluster_names)
    # print(index)
    return index, cluster_names[index]

def getPvalue(pvalues, i, j, N):
    #print(pvalues)
    #print(i)
    #print(j)
    index = (2*N - i - 1) * i / 2 + (j - (i + 1))
    #print(index)
    return pvalues[index]

def combinePvalueStr(name1, name2, pres):
  return name1 + " v.s. " + name2 + ": {0:.3g}".format(pres)

def searchPerPvalue(pvalues, cindex, cnames, C):
    cnum = len(cindex)
    res_list = []
    for i in range(cnum):
        for j in range(i+1, cnum):
            per_p = getPvalue(pvalues, cindex[i], cindex[j], C)
            per_res = combinePvalueStr(cnames[i], cnames[j], per_p)
            res_list.append(per_res)

    return "; ".join(res_list)

def searchPvalues(cnames, meta_pres, otu_pres, attr_pres, cluster_names, search_res):
    cindex, cnames = mapCname2Index(cnames, cluster_names)
    C = len(cluster_names)
    mres = []
    ores = []
    ares = []
    Q = len(meta_pres)
    P = len(otu_pres)
    attr_num = len(attr_pres)
    for i in range(Q):
        per_p = meta_pres[i]
        per_sm = searchPerPvalue(per_p, cindex, cnames, C)
        mres.append(per_sm)

    for i in range(P):
        per_p = otu_pres[i]
        per_so = searchPerPvalue(per_p, cindex, cnames, C)
        ores.append(per_so) 

    for i in range(attr_num):
        per_p = attr_pres[i]
        per_sa = searchPerPvalue(per_p, cindex, cnames, C)
        ares.append(per_sa) 
         
    search_res["meta"] = mres
    search_res["otu"] = ores 
    search_res["attr"] = ares

    return

def parsekLDMResult(result_dir, otu_name_file, meta_name_file, attr_name_file, shape_file, meta_pres, otu_pres,attr_pres, cluster_names):
    # result_dir = target_dir + "/result_dir/"
    print("start read name file !")
    otu_name, ef_name, attr_name = readOTUAndEF(otu_name_file, meta_name_file, attr_name_file)
    print("read name file success ~")
    n, p, q = readMatrixShape(shape_file)
    print("read matrix shape success ~")
    print("dataset info: P:{0} Q:{1} N:{2}".format(p,q,n))
    # print(otu_name)
    # print(ef_name)
    hier_dir = result_dir + "/hierarchical_results/"
    # print(os.listdir(hier_dir)) 
    print("parse recursively ~")
    tree_dict = parseHier(hier_dir, meta_pres, otu_pres, attr_pres, cluster_names)
    print("parse ok ~")
    tree_dict["P"] = p
    tree_dict["Q"] = q
    tree_dict["N"] = n
    tree_dict["otu_name"] = otu_name
    tree_dict["meta_name"] = ef_name
    tree_dict["attr_name"] = attr_name
    tree_dict["attr_num"] = len(attr_name)

    # finish association info
    first_layer = tree_dict["first_layer"]
    for per_layer in first_layer:
        per_asso = per_layer["association"]
        per_otu_otu = per_asso["otu_otu"]
        per_ef_otu = per_asso["ef_otu"]
        per_otu_otu_new = supplyNamesOTUOTU(per_otu_otu, otu_name)
        per_ef_otu_new = supplyNamesEFOTU(per_ef_otu, ef_name, otu_name)
        per_asso["otu_otu"] = per_otu_otu_new
        per_asso["ef_otu"] = per_ef_otu_new
    
    return tree_dict 

def supplyNamesOTUOTU(otu_otu, otu_name):
    otu_otu_new = []
    # print(otu_name)
    for pert in otu_otu:
        # print(pert)
        otu1 = otu_name[pert[0]]
        otu2 = otu_name[pert[1]]
        weight = pert[2]
        otu_otu_new.append([otu1, otu2, weight])

    return otu_otu_new

def supplyNamesEFOTU(ef_otu, ef_name, otu_name):
    ef_otu_new = []
    for pert in ef_otu:
        ef1 = ef_name[pert[0]]
        otu1 = otu_name[pert[1]]
        weight = pert[2]
        ef_otu_new.append([ef1, otu1, weight])
    
    return ef_otu_new
