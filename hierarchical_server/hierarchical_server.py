# -*- coding: utf-8 -*-
from flask import Flask, render_template, request, jsonify
import utils
import argparse
import os
import numpy as np  
import json
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon 
from scipy.stats import mannwhitneyu

parser = argparse.ArgumentParser(description = "visualize the result of kLDM")
parser.add_argument("-on", "--otu_name", type = str, required = True, help = "OTU Name File")
parser.add_argument("-mn", "--meta_name", type = str, required = True, help = "Meta Name File")
parser.add_argument("-an", "--attr_name", type = str, required = True, help = "Attribute Name File")
parser.add_argument("-ms", "--matrix_shape", type = str, required = True, help = "Matrix Shape File")
parser.add_argument("-res", "--kldm_res", type = str, required = True, help = "kLDM result directory")
parser.add_argument("-ot", "--otu_table", type = str, required = True, help = "OTU Table File")
parser.add_argument("-mt", "--meta_table", type = str, required = True, help = "Meta Table File")
parser.add_argument("-at", "--attr_table", type = str, required = True, help = "Attribute Table File")

args = parser.parse_args()
otu_name_file = args.otu_name
meta_name_file = args.meta_name
attr_name_file = args.attr_name
shape_file = args.matrix_shape
result_dir = args.kldm_res
otu_table_file = args.otu_table
meta_table_file = args.meta_table
attr_table_file = args.attr_table
otu_name = None
meta_name = None
attr_name = None
tree_info = None

# read otu and meta table
otu_count = np.loadtxt(otu_table_file)
otu_log = np.log(otu_count + 1)
otu_table = (otu_log.T / otu_log.sum(axis = 1)).T
meta_table = np.loadtxt(meta_table_file)
attr_table = np.loadtxt(attr_table_file)

print("Load OTU, Meta and Attributes Tables Success ~")
hier_dir = result_dir + "/hierarchical_results/"
sample_index_dict = utils.parseSampleIndex(hier_dir)
meta_pvalue_result = []
otu_pvalue_result = []
attr_pvalue_result = []
N, Q = meta_table.shape
_, P = otu_table.shape
_, attr_num = attr_table.shape

meta_p_file = "./meta_pvalue.txt"
otu_p_file = "./otu_pvalue.txt"
attr_p_file = "./attr_pvalue.txt"
cluster_names = sample_index_dict.keys()
cluster_num = len(cluster_names)

if os.path.exists(meta_p_file):
    meta_p = np.loadtxt(meta_p_file)
    meta_pvalue_result = utils.convert2List(meta_p) 
    print("Load EF Pvalue OK~")
else:
    for i in range(Q):
        per_col = meta_table[:,i]
        per_pvalue = []
        for k1 in range(cluster_num):
            for k2 in range(k1 + 1, cluster_num):
                # print(sample_index_dict[cluster_names[k1]])
                group1 = per_col[sample_index_dict[cluster_names[k1]]]
                # print(sample_index_dict[cluster_names[k2]])
                group2 = per_col[sample_index_dict[cluster_names[k2]]]
                per_res = ttest_ind(group1, group2).pvalue
                per_pvalue.append(per_res)
            
        meta_pvalue_result.append(per_pvalue) 

    np.savetxt("./meta_pvalue.txt", np.array(meta_pvalue_result))
    print("Compute EF Pvalue OK~")

if os.path.exists(otu_p_file):
    otu_p = np.loadtxt(otu_p_file)
    otu_pvalue_result = utils.convert2List(otu_p)
    print("Load OTU Pvalue OK~")
else:
    for i in range(P):
        per_col = otu_table[:,i]
        per_pvalue = []
        for k1 in range(cluster_num):
            for k2 in range(k1 + 1, cluster_num):
                group1 = per_col[sample_index_dict[cluster_names[k1]]]
                group2 = per_col[sample_index_dict[cluster_names[k2]]]
                # per_res = wilcoxon(group1, group2).pvalue
                per_res = mannwhitneyu(group1, group2).pvalue
                per_pvalue.append(per_res)
            
        otu_pvalue_result.append(per_pvalue)

    np.savetxt("./otu_pvalue.txt", np.array(otu_pvalue_result))
    print("Compute OTU Pvalue OK~")

if os.path.exists(attr_p_file):
    attr_p = np.loadtxt(attr_p_file)
    attr_pvalue_result = utils.convert2List(attr_p)
    print("Load Attr Pvalue OK~")
else:
    for i in range(attr_num):
        per_col = attr_table[:,i]
        per_pvalue = []
        for k1 in range(cluster_num):
            for k2 in range(k1 + 1, cluster_num):
                group1 = per_col[sample_index_dict[cluster_names[k1]]]
                group2 = per_col[sample_index_dict[cluster_names[k2]]]
                per_res = ttest_ind(group1, group2).pvalue
                per_pvalue.append(per_res)
            
        attr_pvalue_result.append(per_pvalue)
                
    np.savetxt("./attr_pvalue.txt", np.array(attr_pvalue_result))
    print("Compute Attributes Pvalue OK~")

app = Flask(__name__)

@app.route("/getAssoInfo/", methods = ["GET"])
def getAssoInfo():
    dir_name = request.args.get("dir")
    asso_type = request.args.get("type")
    top_ratio = float(request.args.get("ratio"))
    final_asso = None
    global otu_name 
    global meta_name
    if asso_type == "ef":
        ef_file = dir_name + "/ef_otu_association"
        per_asso = utils.parseEFOTU(ef_file, top_ratio)
        final_asso = utils.supplyNamesEFOTU(per_asso, meta_name, otu_name)
    elif asso_type == "otu":
        otu_file = dir_name + "/otu_otu_association"
        per_asso = utils.parseOTUOTU(otu_file, top_ratio)
        final_asso = utils.supplyNamesOTUOTU(per_asso, otu_name)
    elif asso_type == "all":
        ef_file = dir_name + "/ef_otu_association"
        per_ef_asso = utils.parseEFOTU(ef_file, top_ratio)
        final_ef_asso = utils.supplyNamesEFOTU(per_asso, meta_name, otu_name)
        otu_file = dir_name + "/otu_otu_association"
        per_otu_asso = utils.parseOTUOTU(otu_file, top_ratio)
        final_otu_asso = utils.supplyNamesOTUOTU(per_asso, otu_name)
        final_asso = {}
        final_asso["ef_otu"] = final_ef_asso
        final_asso["otu_otu"] = final_otu_asso

    result = {"asso": final_asso}
    return jsonify(result)

@app.route("/getMeanInfo/", methods = ['POST'])
def getMeanInfo():
    origin_data = request.form["value"]
    # print(origin_data)
    json_data = json.loads(origin_data.decode("utf-8"))
    # print(json_data)

    cluster_name = json_data.get("name").strip()
    dir_name = json_data.get("dir")
    # print(cluster_name)
    # print(dir_name)
    current_cnames = json_data.get("cluster_names")
    current_cnames = [per_name.strip() for per_name in current_cnames]
    print("current_names")
    print(current_cnames)
    # read mean info
    mean_info = {}
    utils.readMeanInfo(dir_name, mean_info)
    mean_info["cluster_name"] = cluster_name
    ef_file = dir_name + "/ef_otu_association"
    top_ratio = 10 
    # read association
    global otu_name
    global meta_name
    per_ef_asso = utils.parseEFOTU(ef_file, top_ratio)
    final_ef_asso = utils.supplyNamesEFOTU(per_ef_asso, meta_name, otu_name)
    otu_file = dir_name + "/otu_otu_association"
    per_otu_asso = utils.parseOTUOTU(otu_file, top_ratio)
    final_otu_asso = utils.supplyNamesOTUOTU(per_otu_asso, otu_name)
    # read pvalue 
    search_res = {}
    global meta_pvalue_result
    global otu_pvalue_result 
    global attr_pvalue_result 
    global cluster_names

    utils.searchPvalues(current_cnames, meta_pvalue_result, otu_pvalue_result, attr_pvalue_result, cluster_names, search_res)

    mean_info["ef_otu"] = final_ef_asso
    mean_info["otu_otu"] = final_otu_asso
    mean_info["ratio"] = top_ratio
    mean_info["dir"] = dir_name
    mean_info["meta_p"] = search_res["meta"]
    mean_info["otu_p"] = search_res["otu"]
    mean_info["attr_p"] = search_res["attr"]

    return jsonify(mean_info)

@app.route("/", methods = ['POST', 'GET'])
def welcome():
    print("start parsing ~")
    global tree_info
    if tree_info == None:
        tree_info = utils.parsekLDMResult(result_dir, otu_name_file, meta_name_file, attr_name_file, shape_file, meta_pvalue_result, otu_pvalue_result, attr_pvalue_result, cluster_names)
    print("parsing finished !")
    global otu_name
    global meta_name
    otu_name = tree_info["otu_name"]
    meta_name = tree_info["meta_name"]
    # print(tree_info)
    return render_template("welcome.html", tree_info = tree_info)

if __name__ == "__main__":
    app.debug = True
    app.run(host = "0.0.0.0", port = 8087)
