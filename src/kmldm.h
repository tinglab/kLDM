#ifndef KMLDM_H
#define KMLDM_H

#include "utils.h"
#include "gmm2.h"
#include "mldm.h"

#define SMALL_NOISE 1e-4

// record the otu id and its average size
struct OTUSize {
    int otu_index;
    int otu_zero;
    double otu_size;
};

// record the matched child nodes' index
struct ChildInfo {
    // the index of split result
    int node_index;
    bool is_left;
};

// record previous result before merge 
struct OldCluster {
    VectorXd meta_mean;
    MatrixXd meta_cov;
    std::vector<int> sample_index;
    std::vector<int> otu_index;
    mLDMResult mldm_res;
};

// record the intermediate result of the split process
struct SplitResult {
    int parent_index;
    int child1_index;
    int child2_index;
    // if merge with other cluster, record 
    int layer_num;
    // show if the child is discarded
    bool child1_discard;
    bool child2_discard;
    // store the 2-GMM results
    GMM2Result gmm2_res;
    // store the mldm results
    mLDMResult child1_mldm_res;
    mLDMResult child2_mldm_res;
    // record merge infomation
    std::vector<int> child1_old_list;
    std::vector<int> child2_old_list;
    // selected OTUs of child1
    std::vector<int> child1_otus_list;
    // selected OTUs of child2
    std::vector<int> child2_otus_list;
    // record the other cluster index when merge
    std::vector<ChildInfo> child1_merge_list;
    std::vector<ChildInfo> child2_merge_list;
};

class KmLDM {
private:
    // N * P OTU table
    MatrixXi X;
    // N * Q meta table
    MatrixXd M;
    // samples number
    int N;
    // otus number
    int P;
    // environmental factors number
    int Q;
    // minimum number of samples of one cluster
    int N_min;
    // maximum zero count in sampels of per OTU
    float zero_ratio;
    // noise cluster size
    // when the size of cluster smaller than the noise_size,
    // this cluster is discarded
    // int noise_size;
    // store intermediate results of split process via 2-Gaussian Mixture Model
    std::vector<SplitResult> split_results_list;
    // point to the next position of new split_reuslts_list
    int split_tail;
    // store intermediate results of merge process 
    std::vector<OldCluster> old_results_list;
    // point to the next position of old_results_list
    int old_tail;
    // the threshold of kmeans
    double delta_kmeans;
    // max iterations of kmeans
    int max_iterations_kmeans;
    // repeated search times of kmeans
    // select center randomly
    int repeat_times_kmeans;
    // the threshold of 2-gmm
    double delta_2gmm;
    // max iterations of 2-gmm
    int max_iterations_2gmm;
    // the threshold of mldm
    double delta_mldm;
    // max iterations of mldm
    int max_iterations_mldm;
    // remove sample
    // bool enable_remove_sample;
    // discard the cluster (when the size of cluster < noise_size)
    // bool enable_discard_cluster;

public:
    // constructor
    KmLDM(float zero_r = 0.9, double del_kmeans = 1e-4, int max_iters_kmeans = 500,
    int repeat_kmeans = 10, double del_2gmm = 1e-4, int max_iters_2gmm = 1000, double del_mldm = 1e-3,
    int max_iters_mldm = 2000 /*, int noise_s = 50 *//*, bool remove_sample = false */);

    // enable down sample
    // void enableRemoveSample(bool flag);

    // get otu number
    int getOTUsNumber();

    // set the threshold of min cluster
    // when the size of one cluster small than N_min,
    // split stop
    void setClusterThreshold(int min_size);

    // split and merge process
    void work();
    // split the samples according to the environmental factors, via 2-GMM
    int splitEFPoints(MatrixXd& efs, std::vector<int>& real_index_list, int parent_index, int layer);
    // filter the zero column from matrix
    void filterZeroVarColumn(MatrixXi& x, std::vector<int>& feature_list, MatrixXi& x_filter);

    // read matrix from otu_file
    bool loadOTUs(char* otu_file);
    // read matrix from meta_file
    bool loadMetas(char* meta_file);
    // read otu table and meta table, load into otus and efs respectively
    void loadData(char* otu_file, char* meta_file, char* shape);
    // export the estimated association networks into directory
    void exportResult(char* result_dir);
    // choose the childs candidate for merge
    void findMatchedChilds(int cur_index, VectorXd& direction, VectorXd& center, std::vector<ChildInfo>& candidates);
    // merge childs candidate
    void mergeMatchedChilds(VectorXd& left_center, VectorXd& right_center, std::vector<ChildInfo>& left_nodes, std::vector<ChildInfo>& right_nodes);

    // get merged otus when merge clusters
    void getMergedOTUsList(std::vector<int>& left_otus_list, std::vector<int>& right_otus_list,
    std::vector<int>& left_samples_list, std::vector<int>& right_samples_list, std::vector<int>& merged_otus_list,
    bool& left_reestimate, bool& right_reestimate, bool& otu_not_match);

    // construct otu and meta data according to otu list and sample list
    void constructData(std::vector<int>& samples_list, std::vector<int>& otus_list, MatrixXi& otu_data, MatrixXd& ef_data);

    // export hierarchical result 
    void exportHierarchicalNodes(int cluster_number, int node_index, bool is_left, std::string result_dir_name);
};

#endif
