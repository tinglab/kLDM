#include "kmldm.h"

#include <fstream>

#define DEBUG
#define VERBOSE

KmLDM::KmLDM(float zero_r, double del_kmeans, int max_iters_kmeans,
    int repeat_kmeans, double del_2gmm, int max_iters_2gmm, double del_mldm,
    int max_iters_mldm /*, int noise_s*//*, bool remove_sample*/) {
    // set the next position of the splitResults
    split_tail = 0;
    old_tail = 0;
    // init parameters
    zero_ratio = zero_r;
    delta_kmeans = del_kmeans;
    max_iterations_kmeans = max_iters_kmeans;
    repeat_times_kmeans = repeat_kmeans;
    delta_2gmm = del_2gmm;
    max_iterations_2gmm = max_iters_2gmm;
    delta_mldm = del_mldm;
    max_iterations_mldm = max_iters_mldm;
    // noise_size = noise_s;
   // enable_remove_sample = remove_sample;
}

// enable remove samples
/*
void KmLDM::enableRemoveSample(bool flag) {
    enable_remove_sample = flag;
}*/

// set the split threshold of clusters
void KmLDM::setClusterThreshold(int min_size) {
    N_min = min_size;
}

// split and merge process
// Step 1: split whole datapoints according to environmental factors, via 2-GMM
// Step 2: estimate parameters within every clusters
// Step 3: merger similar clusters
void KmLDM::work() {
    //Step 1: split datapoints recursively
#ifdef VERBOSE
    std::cout << "====================== Split Meta Data Recursively ======================" << std::endl;
    std::cout << "When the number of samples of the cluster < " << N_min << ", the split stop !" << std::endl;
#endif
    // The true index of environmental factors within matrix M
    std::vector<int> initial_index_list;
    for (int i = 0; i < M.rows(); ++ i) {
        initial_index_list.push_back(i);
    }

    // begin split process
    int root_index = splitEFPoints(M, initial_index_list, -1, 0);
#ifdef VERBOSE
    std::cout << "Split Meta Data Finished ~" << std::endl;
#endif

}

// split the samples according to the environmental factors, via 2-GMM
int KmLDM::splitEFPoints(MatrixXd& efs, std::vector<int>& real_index_list, int parent_index, int layer) {
    int split_index = split_tail;
    split_tail ++;
    GMM2Result cur_result;
    // 2 gaussian mixture model clustering
    GMM2(efs, cur_result, delta_2gmm, delta_kmeans, Q, max_iterations_2gmm, max_iterations_kmeans, repeat_times_kmeans);
    std::cout << "after gmm2";
    int cluster1_size = cur_result.cluster1_index_list.size();
    int cluster2_size = cur_result.cluster2_index_list.size();
    // convert the index in GMM2Result into real index
    for (int i = 0; i < cluster1_size; ++i) {
        cur_result.cluster1_index_list[i] = real_index_list[cur_result.cluster1_index_list[i]];
    }
    for (int i = 0; i < cluster2_size; ++i) {
        cur_result.cluster2_index_list[i] = real_index_list[cur_result.cluster2_index_list[i]];
    }

#ifdef DEBUG
    if (parent_index != -1) {
        std::cout << "Cluster [" << split_index << "] (size: " << cluster1_size + cluster2_size << "), its parent is Cluster [" << parent_index << "]" << std::endl;
    } else {
        std::cout << "Cluster [" << split_index << "] (size: " << cluster1_size + cluster2_size << "), it is the Root" << std::endl;
    }
#endif

    // save the split result into SplitResult
    SplitResult cur_split;
    cur_split.parent_index = parent_index;
    cur_split.gmm2_res = cur_result;
    cur_split.child1_index = -1;
    cur_split.child1_discard = false;
    cur_split.child2_index = -1;
    cur_split.child2_discard = false;
    cur_split.layer_num = layer;
    split_results_list.push_back(cur_split);

    // split cur cluster into two sub-clusters
    // generate partial efs for left child
    MatrixXd efs_1 = MatrixXd::Zero(cluster1_size, Q);
    for (int i = 0; i < cluster1_size; ++ i) {
        efs_1.row(i) = M.row(cur_result.cluster1_index_list[i]);
    }

    if (cluster1_size > N_min) {
        cur_split.child1_index = splitEFPoints(efs_1, cur_result.cluster1_index_list, split_index, layer + 1);
    } else { // left child is the leaf
#ifdef DEBUG
        std::cout << "Left child of cluster [" << split_index << "] (size: " << cluster1_size << ") is the leaf, split stop !" << std::endl;
#endif
        // judge if it is the noise cluster
        /* if (cluster1_size <= noise_size) {
            cur_split.child1_discard = true;
        } else { */

        cur_split.child1_discard = false;

        // estimate association networks of left child
        // generate partial otus for left child
        MatrixXi otus_1 = MatrixXi::Zero(cluster1_size, P);
        for (int i = 0; i < cluster1_size; ++ i) {
            otus_1.row(i) = X.row(cur_result.cluster1_index_list[i]);
        }

        // test if there is zero-variance OTU
        std::vector<int> otus_1_list;
        MatrixXi otus_1_filter;
        filterZeroVarColumn(otus_1, otus_1_list, otus_1_filter);

        // if all otus or all samples are filtered, discard the cluster
        if (otus_1_list.size() <= 0 /*|| otus_1_filter.rows() <= noise_size */) {
            cur_split.child1_discard = true;
        } else {
            cur_split.child1_discard = false;

#ifdef DEBUG
            std::cout << "Left child of cluster [" << split_index << "]'s efficient otus' num is " << otus_1_list.size() << std::endl;
#endif

            // estimate association network via mLDM
            mLDM child1_mldm(otus_1_filter, efs_1, otus_1_filter.rows(), otus_1_filter.cols(), Q, max_iterations_mldm, delta_mldm, false /*, noise_size*/);
            // child1_mldm.enableRemoveSample(enable_remove_sample);
            child1_mldm.work();
            // save result for left child
            // if can not estimate associations, discard the cluster
            if (! child1_mldm.best_result.exists) {
#ifdef DEBUG
                std::cout << "child1 not exist" << std::endl;
#endif
                cur_split.child1_discard = true;
            } else {
#ifdef DEBUG
                std::cout << "child1 exist ~" << std::endl;
#endif
                // update the remained sampels
                /*int child1_removed_size = child1_mldm.removed_samples_index.size();
                if (child1_removed_size > 0) {
                    for (int i = 0; i < child1_removed_size; ++i) {
                        cur_split.gmm2_res.cluster1_index_list.erase(cur_split.gmm2_res.cluster1_index_list.begin()
                            + child1_mldm.removed_samples_index[i]);
                    }
                    // update gmm2 result
                    // update mean and cov
                    MatrixXd child1_meta_data = MatrixXd::Zero(cluster1_size - child1_removed_size, Q);
                    for (int i = 0; i < cur_split.gmm2_res.cluster1_index_list.size(); ++i) {
                        child1_meta_data.row(i) = M.row(cur_result.cluster1_index_list[i]);
                    }
                    VectorXd child1_meta_mean = child1_meta_data.colwise().mean().transpose();
                    MatrixXd child1_meta_cov = MatrixXd::Zero(Q, Q);
                    computeCov(child1_meta_data, child1_meta_cov);
                    cur_split.gmm2_res.cluster1_mean = child1_meta_mean;
                    cur_split.gmm2_res.cluster1_cov = child1_meta_cov;
                }/
                std::cout << "child1 remove " << child1_removed_size << " samples" << std::endl; */

                cur_split.child1_discard = false;
                cur_split.child1_otus_list = otus_1_list;
                cur_split.child1_mldm_res = child1_mldm.best_result;
            }
        }

        // }

        // show if this child is discarded
        if (cur_split.child1_discard) {

#ifdef DEBUG
            std::cout << "Left child of cluster [" << split_index << "] (size: " << cluster1_size << ") is discarded !" << std::endl;
#endif
        }
    }

    // generate partial efs for right child
    MatrixXd efs_2 = MatrixXd::Zero(cluster2_size, Q);
    for (int i = 0; i < cluster2_size; ++ i) {
        efs_2.row(i) = M.row(cur_result.cluster2_index_list[i]);
    }

    if (cluster2_size > N_min) {
        cur_split.child2_index = splitEFPoints(efs_2, cur_result.cluster2_index_list, split_index, layer + 1);
    } else { // right child is the leaf
#ifdef DEBUG
        std::cout << "Right child of cluster [" << split_index << "] (size: " << cluster2_size << ") is the leaf, split stop !" << std::endl;
#endif

        /* if (cluster2_size < noise_size) {
            cur_split.child2_discard = true;

        } else { */

        cur_split.child2_discard = false;
        // estimate association networks of right child
        // generate partial otus for right child
        MatrixXi otus_2 = MatrixXi::Zero(cluster2_size, P);
        for (int i = 0; i < cluster2_size; ++ i) {
            otus_2.row(i) = X.row(cur_result.cluster2_index_list[i]);
        }
        // test if there is zero-variance OTU
        std::vector<int> otus_2_list;
        MatrixXi otus_2_filter;
        filterZeroVarColumn(otus_2, otus_2_list, otus_2_filter);

        // if all otus or all samples are filtered, discard the cluster
        if (otus_2_list.size() <= 0 /*|| otus_2_filter.rows() <= noise_size */) {
            cur_split.child2_discard = true;
        } else {
            cur_split.child2_discard = false;
#ifdef DEBUG
            std::cout << "Right child of cluster [" << split_index << "]'s efficient otus' num is " << otus_2_list.size() << std::endl;
#endif

            // estimate association network via mLDM
            mLDM child2_mldm(otus_2_filter, efs_2, otus_2_filter.rows(), otus_2_filter.cols(), Q, max_iterations_mldm, delta_mldm, false /*, noise_size*/);
            // child2_mldm.enableRemoveSample(enable_remove_sample);
            child2_mldm.work();
            // save result for right child
            if (! child2_mldm.best_result.exists) {
#ifdef DEBUG
                std::cout << "child2 not exist" << std::endl;
#endif
                cur_split.child2_discard = true;
            } else {
#ifdef DEBUG
                std::cout << "child2 exist ~" << std::endl;
#endif
                // update the remained sampels
                /*int child2_removed_size = child2_mldm.removed_samples_index.size();
                if (child2_removed_size > 0) {
                    for (int i = 0; i < child2_removed_size; ++i) {
                        cur_split.gmm2_res.cluster2_index_list.erase(cur_split.gmm2_res.cluster2_index_list.begin()
                            + child2_mldm.removed_samples_index[i]);
                    }
                    // update gmm2 result
                    // update mean and cov
                    MatrixXd child2_meta_data = MatrixXd::Zero(cluster2_size - child2_removed_size, Q);
                    for (int i = 0; i < cur_split.gmm2_res.cluster2_index_list.size(); ++i) {
                        child2_meta_data.row(i) = M.row(cur_result.cluster2_index_list[i]);
                    }
                    VectorXd child2_meta_mean = child2_meta_data.colwise().mean().transpose();
                    MatrixXd child2_meta_cov = MatrixXd::Zero(Q, Q);
                    computeCov(child2_meta_data, child2_meta_cov);
                    cur_split.gmm2_res.cluster2_mean = child2_meta_mean;
                    cur_split.gmm2_res.cluster2_cov = child2_meta_cov;
                }
                std::cout << "child2 remove " << child2_removed_size << " samples" << std::endl; */

                cur_split.child2_discard = false;
                cur_split.child2_otus_list = otus_2_list;
                cur_split.child2_mldm_res = child2_mldm.best_result;
            }
        }

        // }

        if (cur_split.child2_discard) {
#ifdef DEBUG
            std::cout << "Right child of cluster [" << split_index << "] (size: " << cluster2_size << ") is discarded !" << std::endl;
#endif
        }
    }

    // if it is the leaf cluster, then run mLDM to estimate association network
    if (cluster1_size < N_min && cluster2_size < N_min) {
#ifdef DEBUG
        std::cout << "Cluster [" << split_index << "] is the leaf, split stop !" << std::endl;
#endif
    }

    // update the cur split result
    split_results_list[split_index] = cur_split;

    // sample list 
 /*   std::cout << "cluster child1 " << "sample list" << std::endl;
    for (int i = 0; i < cluster1_size; ++i) {
        std::cout << cur_split.gmm2_res.cluster1_index_list[i] << " ";
    }
    std::cout << "cluster child2 " << "sample list" << std::endl;
    for (int i = 0; i < cluster2_size; ++i) {
        std::cout << cur_split.gmm2_res.cluster2_index_list[i] << " ";
    }
*/
    // when split stop, merge the children recursively
    // left nodes and right nodes
    std::vector<ChildInfo> left_nodes;
    VectorXd left_center = cur_split.gmm2_res.cluster1_mean;
    std::vector<ChildInfo> right_nodes;
    VectorXd right_center = cur_split.gmm2_res.cluster2_mean;

    // left direction
    VectorXd left_direction = right_center - left_center;
    // right direction
    VectorXd right_direction = - left_direction;

#ifdef DEBUG
    std::cout << "/////////////////// find matched process begin: " << std::endl;
#endif
    // find left candidate
    if (cur_split.child1_index == -1) {
        if (! cur_split.child1_discard) {
            ChildInfo left_leaf;
            left_leaf.node_index = split_index;
            left_leaf.is_left = true;
            left_nodes.push_back(left_leaf);
#ifdef DEBUG
            std::cout << "find left child of Cluster [" << left_leaf.node_index << "] " << std::endl;
#endif
        } else {
#ifdef DEBUG
            std::cout << "left child of Cluster [" << split_index << "] not match" << std::endl;
#endif
        }
    } else {
        findMatchedChilds(cur_split.child1_index, left_direction, left_center, left_nodes);
    }
    // find right candidate
    if (cur_split.child2_index == -1) {
        if (! cur_split.child2_discard) {
            ChildInfo right_leaf;
            right_leaf.node_index = split_index;
            right_leaf.is_left = false;
            right_nodes.push_back(right_leaf);
#ifdef DEBUG
            std::cout << "find right child of Cluster [" << right_leaf.node_index << "] "<< std::endl;
#endif
        } else {
#ifdef DEBUG
            std::cout << "right child of Cluster [" << split_index << "] not match"<< std::endl;
#endif
        }
    } else {
        findMatchedChilds(cur_split.child2_index, right_direction, right_center, right_nodes);
    }
#ifdef DEBUG
    std::cout << "left candidates's size: " << left_nodes.size() << std::endl;
    std::cout << "right candidates's size: " << right_nodes.size() << std::endl;
#endif
#ifdef DEBUG
    std::cout << "/////////////////// find matched process end~ " << std::endl;
#endif
    // try to merge small clusters to the big
    // pop two nodes from left and right list respectively
    mergeMatchedChilds(left_center, right_center, left_nodes, right_nodes);

    return split_index;
}

// try to merge two clusters' otus
void KmLDM::getMergedOTUsList(std::vector<int>& left_otus_list, std::vector<int>& right_otus_list,
    std::vector<int>& left_samples_list, std::vector<int>& right_samples_list, std::vector<int>& merged_otus_list,
    bool& left_reestimate, bool& right_reestimate, bool& otu_not_match) {
    std::vector<int> extra_otus_1;
    std::vector<int> extra_otus_2;
    std::vector<int> intersect_otus;
    otu_not_match = false;

    // compute the difference set of cluster1 and cluster2 respectively
    // for cluster 1
    getDiffList(right_otus_list, left_otus_list, extra_otus_1);
    std::cout << "extra otus 1 size is " << extra_otus_1.size() << std::endl;
    // for cluster 2
    getDiffList(left_otus_list, right_otus_list, extra_otus_2);
    std::cout << "extra otus 2 size is " << extra_otus_2.size() << std::endl;
    // get intersection set
    getIntersectList(left_otus_list, right_otus_list, intersect_otus);
    std::cout << "intersection otus size is " << intersect_otus.size() << std::endl;

    // filter zero variance otus
    // if filter more than 90% otus of extra_otus_1 or extra_otus_2
    // stop merge
    // filter cluster 1
    int n1 = left_samples_list.size();
    int n2 = right_samples_list.size();
    int extra1_remove = 0;
    int extra2_remove = 0;
    for (int i = 0; i < extra_otus_1.size(); ++i) {
        int otu_index = extra_otus_1[i];
        VectorXd xi = VectorXd::Zero(n1);
        for (int j = 0; j < n1; ++j) {
            xi[j] = X(left_samples_list[j], otu_index);
        }
        // compute the variance and zero num
        double var_i = computeVar(xi);
        int zero_i = zeroNum(xi);

        if (isZero(var_i) || (zero_i * 1.0 / n1 > zero_ratio)) {
            std::cout << "left node remove otu: " << i << " zero: " << zero_i << " var: " << var_i << std::endl;
            extra_otus_1.erase(extra_otus_1.begin() + i);
            extra1_remove += 1;
            i -= 1;
            otu_not_match = true;
            // return;
        }
    }
    // if remaining otus less than zero_ratio, don't merge cluster
    std::cout << "extra 1 remove " << extra1_remove << std::endl;
    /*if (extra1_remove * 1.0 / n1 > (1 - zero_ratio)) {
        std::cout << "right cluster's extra otus may not exist in left cluster ! can not merge them !" << std::endl;
        return;
    }*/
    if (extra_otus_1.size() > 0) {
        left_reestimate = true;
    }
    if (extra1_remove != 0) {
        right_reestimate = true;
    }
    // filter cluster 2
    for (int i = 0; i < extra_otus_2.size(); ++i) {
        int otu_index = extra_otus_2[i];
        VectorXd xi = VectorXd::Zero(n2);
        for (int j = 0; j < n2; ++j) {
            xi[j] = X(right_samples_list[j], otu_index);
        }
        // compute the variance and zero num
        double var_i = computeVar(xi);
        int zero_i = zeroNum(xi);

        if (isZero(var_i) || (zero_i * 1.0 / n2 > zero_ratio)) {
            // std::cout << i << " " << zero_i << " " << var_i << std::endl;
            std::cout << "right node remove otu: " << i << " zero: " << zero_i << " var: " << var_i << std::endl;
            extra_otus_2.erase(extra_otus_2.begin() + i);
            extra2_remove += 1;
            i -= 1;
            otu_not_match = true;
            // return;
        }
    }
    // if remaining otus less than zero_ratio, don't merge cluster
    std::cout << "extra 2 remove " << extra2_remove << std::endl;
    /*if (extra2_remove * 1.0 / n2 > (1 - zero_ratio)) {
        std::cout << "left cluster's extra otus may not exist in right cluster ! can not merge them !" << std::endl;
        return;
    }*/
    if (extra_otus_2.size() > 0) {
        right_reestimate = true;
    }
    if (extra2_remove != 0 ){
        left_reestimate = true;
    }

    // get merged otus
    std::vector<int> extra_otus;
    getUnionList(extra_otus_1, extra_otus_2, extra_otus);
    getUnionList(intersect_otus, extra_otus, merged_otus_list);

    return;
}

// construct otu and meta data according to otu list and sample list
void KmLDM::constructData(std::vector<int>& samples_list, std::vector<int>& otus_list, MatrixXi& otu_data, MatrixXd& ef_data) {
    int merged_n = samples_list.size();
    int merged_p = otus_list.size();
    for (int k = 0; k < merged_n; ++k) {
        // new meta data
        ef_data.row(k) = M.row(samples_list[k]);
        // new otu data
        for (int g = 0; g < merged_p; ++g) {
            otu_data(k, g) = X(samples_list[k], otus_list[g]);
        }
        // add noise when the sample is very small
        double per_sample_num = otu_data.row(k).sum();
        if ((per_sample_num < (1 - zero_ratio) * merged_p)) {
            // generate random number to decide adding 1 count to which otu
            VectorXd random_ratio = VectorXd::Random(merged_p);
            // set the prob of add 1 count to zero_ratio
            for (int j = 0; j < random_ratio.size(); ++j) {
                if (random_ratio[j] > zero_ratio || random_ratio[j] < -zero_ratio) {
                    otu_data(k, j) += 1;
                }
            }
        }
    }
}

// pick one node from left and one from right
// merge their data to run mLDM
// if bic(merge) < bic(1) + bic(2), merge two nodes
// else continue
void KmLDM::mergeMatchedChilds(VectorXd& left_center, VectorXd& right_center, std::vector<ChildInfo>& left_nodes, std::vector<ChildInfo>& right_nodes) {
    std::vector<int> left_otus_list;
    std::vector<int> left_samples_list;
    bool left_position;
    bool right_position;
    std::vector<int> right_otus_list;
    std::vector<int> right_samples_list;
    std::vector<int> merged_otus_list;
    std::vector<int> merged_samples_list;

    // record EBIC value
    double EBIC_left;
    double EBIC_right;
    double EBIC_merged;
    // record obj value
    double OBJ_left;
    double OBJ_right;
    double OBJ_merged;

#ifdef DEBUG
    std::cout << "+++++++++++++++++++ merge process begin: " << std::endl;
#endif

    while (! left_nodes.empty()) {
        ChildInfo left_node = left_nodes.back();
#ifdef DEBUG
        if (left_node.is_left) {
            std::cout << "Left Node: " << " left child of Cluster [" << left_node.node_index << "]" << std::endl;
        } else {
            std::cout << "Left Node: " << " right child of Cluster [" << left_node.node_index << "]" << std::endl;
        }
#endif
        // clear the lists
        left_otus_list.clear();
        left_samples_list.clear();
        // get left node otus and samples list
        if (left_node.is_left) {
            left_position = true;
            left_otus_list = split_results_list[left_node.node_index].child1_otus_list;
            left_samples_list = split_results_list[left_node.node_index].gmm2_res.cluster1_index_list;
            EBIC_left = split_results_list[left_node.node_index].child1_mldm_res.EBIC;
	        OBJ_left = split_results_list[left_node.node_index].child1_mldm_res.obj;
        } else {
            left_position = false;
            left_otus_list = split_results_list[left_node.node_index].child2_otus_list;
            left_samples_list = split_results_list[left_node.node_index].gmm2_res.cluster2_index_list;
            EBIC_left = split_results_list[left_node.node_index].child2_mldm_res.EBIC;
	        OBJ_left = split_results_list[left_node.node_index].child2_mldm_res.obj;
        }

        bool exist = false;
        for (int i = 0; i < right_nodes.size(); ++ i) {
            ChildInfo right_node = right_nodes[i];
            // clear right lists
            right_otus_list.clear();
            right_samples_list.clear();
            // get right node otus and samples list
            if (right_node.is_left) {
                right_position = true;
                right_otus_list = split_results_list[right_node.node_index].child1_otus_list;
                right_samples_list = split_results_list[right_node.node_index].gmm2_res.cluster1_index_list;
                EBIC_right = split_results_list[right_node.node_index].child1_mldm_res.EBIC;
	            OBJ_right = split_results_list[right_node.node_index].child1_mldm_res.obj;
            } else {
                right_position = false;
                right_otus_list = split_results_list[right_node.node_index].child2_otus_list;
                right_samples_list = split_results_list[right_node.node_index].gmm2_res.cluster2_index_list;
                EBIC_right = split_results_list[right_node.node_index].child2_mldm_res.EBIC;
	            OBJ_right = split_results_list[right_node.node_index].child2_mldm_res.obj;
            }

            // merge their data to compute the BIC
            // get merged otus list
            merged_otus_list.clear();
            // getUnionList(left_otus_list, right_otus_list, merged_otus_list);
            bool left_reestimate = false, right_reestimate = false;
            bool otu_not_match = false;
            getMergedOTUsList(left_otus_list, right_otus_list, left_samples_list,
                right_samples_list, merged_otus_list, left_reestimate, right_reestimate, otu_not_match);
            // otu must match
            if (otu_not_match) {
                std::cout << "left and right otu not match !" << std::endl;
                continue;
            }
            if (merged_otus_list.size() <= 0) {
                continue;
            }
            std::cout << "reestimate result: left: " << left_reestimate << " right: " << right_reestimate << std::endl;
            // get merged samples list
            merged_samples_list.clear();
            getUnionList(left_samples_list, right_samples_list, merged_samples_list);

            // construct otu and meta data
            int merged_n = merged_samples_list.size();
            int merged_p = merged_otus_list.size();
            MatrixXd merged_meta_data = MatrixXd::Zero(merged_n, Q);
            MatrixXi merged_otu_data = MatrixXi::Zero(merged_n, merged_p);
            constructData(merged_samples_list, merged_otus_list, merged_otu_data, merged_meta_data);

            // run mLDM for the dataset
            // estimate association network via mLDM
            mLDM merged_mldm(merged_otu_data, merged_meta_data, merged_n, merged_p, Q, max_iterations_mldm, delta_mldm, false /*, noise_size*/);
            // when merge cluster, and allow remove sample
            // merged_mldm.enableRemoveSample(true);
            merged_mldm.work();
            if (! merged_mldm.best_result.exists) {
                continue;
            }
            // check if has removed some samples
            /*if (merged_mldm.removed_samples_index.size() > 0) {
                for (int k = 0 ;k < merged_mldm.removed_samples_index.size(); ++k) {
                    merged_samples_list.erase(merged_samples_list.begin() +
                        merged_mldm.removed_samples_index[k]);
                }
                std::cout << "merged cluster has removed [" << merged_mldm.removed_samples_index.size() << "] samples" << std::endl;
                // check if left samples need to be updated
                std::vector<int> new_left_samples;
                // get intersection between left samples and merged samples
                getIntersectList(left_samples_list, merged_samples_list, new_left_samples);
                int new_left_size = new_left_samples.size();
                if (new_left_size != left_samples_list.size()) {
                    left_reestimate = true;
                    if (new_left_size * 1.0 / left_samples_list.size() < zero_ratio) {
                        std::cout << "left remove too many samples ! stop merge !" << std::endl;
                        continue;
                    } else {
                        left_samples_list = new_left_samples;
                        std::cout << "update left cluster samples list" << std::endl;
                        std::cout << "new left samples' size is: " << left_samples_list.size() << std::endl;
                    }
                }
                // check if right samples need to be updated
                std::vector<int> new_right_samples;
                // get intersection between left samples and merged samples
                getIntersectList(right_samples_list, merged_samples_list, new_right_samples);
                int new_right_size = new_right_samples.size();
                if (new_right_size != right_samples_list.size()) {
                    right_reestimate = true;
                    if (new_right_size * 1.0 / right_samples_list.size() < zero_ratio) {
                        std::cout << "right remove too many samples ! stop merge !" << std::endl;
                        continue;
                    } else {
                        right_samples_list = new_right_samples;
                        std::cout << "update right cluster samples list" << std::endl;
                        std::cout << "new right samples' size is: " << right_samples_list.size() << std::endl;
                    }
                }
            }*/

            EBIC_merged = merged_mldm.best_result.EBIC;
	        OBJ_merged = merged_mldm.best_result.obj;

            // check left or right cluster if they need to be re-estimated
            if (left_reestimate) {
                std::cout << "re estimate left cluster begin~" << std::endl;
                MatrixXd left_efs_data = MatrixXd::Zero(left_samples_list.size(), Q);
                MatrixXi left_otus_data = MatrixXi::Zero(left_samples_list.size(), merged_p);
                constructData(left_samples_list, merged_otus_list, left_otus_data, left_efs_data);

                mLDM left_mldm(left_otus_data, left_efs_data, left_samples_list.size(), merged_p, Q, max_iterations_mldm, delta_mldm, false /*, noise_size*/);
                // when merge cluster, not allow remove sample
                // left_mldm.enableRemoveSample(false);
                left_mldm.work();
                if (! left_mldm.best_result.exists) {
                    std::cout << "re estimate left cluster failed !!!" << std::endl;
                    continue;
                }
                std::cout << "re estimate left cluster finished = = " << std::endl;
                EBIC_left = left_mldm.best_result.EBIC;
                OBJ_left = left_mldm.best_result.obj;
            }
            if (right_reestimate) {
                std::cout << "re estimate right cluster begin~" << std::endl;
                MatrixXd right_efs_data = MatrixXd::Zero(right_samples_list.size(), Q);
                MatrixXi right_otus_data = MatrixXi::Zero(right_samples_list.size(), merged_p);
                constructData(right_samples_list, merged_otus_list, right_otus_data, right_efs_data);

                mLDM right_mldm(right_otus_data, right_efs_data, right_samples_list.size(), merged_p, Q, max_iterations_mldm, delta_mldm, false /*, noise_size*/);
                // when merge cluster, not allow remove sample
                // right_mldm.enableRemoveSample(false);
                right_mldm.work();
                if (! right_mldm.best_result.exists) {
                    std::cout << "re estimate right cluster failed !!!" << std::endl;
                    continue;
                }
                std::cout << "re estimate right cluster finished = = " << std::endl;
                EBIC_right = right_mldm.best_result.EBIC;
                OBJ_right = right_mldm.best_result.obj;
            }

	        std::cout << "EBIC merged: " << EBIC_merged << " EBIC left: " << EBIC_left << " EBIC right: " << EBIC_right << std::endl;
	        std::cout << "OBJ merged: " << OBJ_merged << " OBJ left: " << OBJ_left << " OBJ right: " << OBJ_right << std::endl;
            // if can merge, update split_results_list, left_nodes and right_nodes
            // break the loop
            if (EBIC_merged <= EBIC_left + EBIC_right) {
#ifdef DEBUG
        if (right_node.is_left) {
            std::cout << "Matched Right Node: " << " left child of Cluster [" << right_node.node_index << "]" << std::endl;
        } else {
            std::cout << "Matched Right Node: " << " right child of Cluster [" << right_node.node_index << "]" << std::endl;
        }
#endif
                exist = true;
                // if left and right nodes don't belong to the same parent
                // if (left_node.node_index != right_node.node_index) {
                    // compute mean and covariance of new cluster
                    VectorXd merged_meta_mean = merged_meta_data.colwise().mean().transpose();
                    MatrixXd merged_meta_cov = MatrixXd::Zero(Q, Q);
                    computeCov(merged_meta_data, merged_meta_cov);
                    // put the cluster into cloest parent
                    if ((merged_meta_mean - left_center).squaredNorm()
                        <= (merged_meta_mean - right_center).squaredNorm()) { // merge into left

                        // discard the right node
#ifdef DEBUG
                        if (right_position) {
                            std::cout << "Merge left child of node " << right_node.node_index << " into ";
                        } else {
                            std::cout << "Merge right child of node " << right_node.node_index << " into ";
                        }
                        if (left_position) {
                            std::cout << "left child of node " << left_node.node_index << std::endl;
                        } else {
                            std::cout << "right child of node " << left_node.node_index << std::endl;
                        }
#endif
        
                        if (right_position) {
                            split_results_list[right_node.node_index].child1_discard = true;
                        } else {
                            split_results_list[right_node.node_index].child2_discard = true;
                        }

                        right_nodes.erase(right_nodes.begin() + i);

                        // update the left node
                        if (left_position) {
                            // record merge info 
                            // record the other cluster index
                            split_results_list[left_node.node_index].child1_merge_list.push_back(right_node);
                            // record current cluster old result 
                            OldCluster left_old;
                            left_old.sample_index = left_samples_list;
                            left_old.otu_index = left_otus_list;
                            left_old.meta_mean = split_results_list[left_node.node_index].gmm2_res.cluster1_mean;
                            left_old.meta_cov = split_results_list[left_node.node_index].gmm2_res.cluster1_cov;
                            left_old.mldm_res = split_results_list[left_node.node_index].child1_mldm_res;
                            old_results_list.push_back(left_old);
                            split_results_list[left_node.node_index].child1_old_list.push_back(old_tail);
                            old_tail ++;

                            // update otus list
                            split_results_list[left_node.node_index].child1_otus_list = merged_otus_list;
                            // update kmldm result
                            split_results_list[left_node.node_index].child1_mldm_res = merged_mldm.best_result;
                            // update the gmm result
                            split_results_list[left_node.node_index].gmm2_res.cluster1_mean = merged_meta_mean;
                            split_results_list[left_node.node_index].gmm2_res.cluster1_cov = merged_meta_cov;
                            split_results_list[left_node.node_index].gmm2_res.cluster1_index_list = merged_samples_list;
                            

                        } else {
                            // record merge info 
                            // record the other cluster index
                            split_results_list[left_node.node_index].child2_merge_list.push_back(right_node);
                            // record current cluster old result 
                            OldCluster left_old;
                            left_old.sample_index = left_samples_list;
                            left_old.otu_index = left_otus_list;
                            left_old.meta_mean = split_results_list[left_node.node_index].gmm2_res.cluster2_mean;
                            left_old.meta_cov = split_results_list[left_node.node_index].gmm2_res.cluster2_cov;
                            left_old.mldm_res = split_results_list[left_node.node_index].child2_mldm_res;
                            old_results_list.push_back(left_old);
                            split_results_list[left_node.node_index].child2_old_list.push_back(old_tail);
                            old_tail ++;

                            // update otus list
                            split_results_list[left_node.node_index].child2_otus_list = merged_otus_list;
                            // update kmldm result
                            split_results_list[left_node.node_index].child2_mldm_res = merged_mldm.best_result;
                            // update the gmm result
                            split_results_list[left_node.node_index].gmm2_res.cluster2_mean = merged_meta_mean;
                            split_results_list[left_node.node_index].gmm2_res.cluster2_cov = merged_meta_cov;
                            split_results_list[left_node.node_index].gmm2_res.cluster2_index_list = merged_samples_list;
                        }
                    } else { // merge into right
                        // discard the left node
#ifdef DEBUG
                        if (left_position) {
                            std::cout << "Merge left child of node " << left_node.node_index << " into ";
                        } else {
                            std::cout << "Merge right child of node " << left_node.node_index << " into ";
                        }
                        if (right_position) {
                            std::cout << "left child of node " << right_node.node_index << std::endl;
                        } else {
                            std::cout << "right child of node " << right_node.node_index << std::endl;
                        }
#endif
                        if (left_position) {
                            split_results_list[left_node.node_index].child1_discard = true;
                        } else {
                            split_results_list[left_node.node_index].child2_discard = true;
                        }
                        left_nodes.pop_back();

                        // update the right node
                        if (right_position) {
                            // record merge info 
                            // record the other node index
                            split_results_list[right_node.node_index].child1_merge_list.push_back(left_node);
                            // record current cluster old result 
                            OldCluster right_old;
                            right_old.sample_index = right_samples_list;
                            right_old.otu_index = right_otus_list;
                            right_old.meta_mean = split_results_list[right_node.node_index].gmm2_res.cluster1_mean;
                            right_old.meta_cov = split_results_list[right_node.node_index].gmm2_res.cluster1_cov;
                            right_old.mldm_res = split_results_list[right_node.node_index].child1_mldm_res;
                            old_results_list.push_back(right_old);
                            split_results_list[right_node.node_index].child1_old_list.push_back(old_tail);
                            old_tail ++;

                            // update otus list
                            split_results_list[right_node.node_index].child1_otus_list = merged_otus_list;
                            // update kmldm result
                            split_results_list[right_node.node_index].child1_mldm_res = merged_mldm.best_result;
                            // update the gmm result
                            split_results_list[right_node.node_index].gmm2_res.cluster1_mean = merged_meta_mean;
                            split_results_list[right_node.node_index].gmm2_res.cluster1_cov = merged_meta_cov;
                            split_results_list[right_node.node_index].gmm2_res.cluster1_index_list = merged_samples_list;
                        } else {
                            // record merge info 
                            // record the other cluster index 
                            split_results_list[right_node.node_index].child2_merge_list.push_back(left_node);
                            // record current cluster old result 
                            OldCluster right_old;
                            right_old.sample_index = right_samples_list;
                            right_old.otu_index = right_otus_list;
                            right_old.meta_mean = split_results_list[right_node.node_index].gmm2_res.cluster2_mean;
                            right_old.meta_cov = split_results_list[right_node.node_index].gmm2_res.cluster2_cov;
                            right_old.mldm_res = split_results_list[right_node.node_index].child2_mldm_res;
                            old_results_list.push_back(right_old);
                            split_results_list[right_node.node_index].child2_old_list.push_back(old_tail);
                            old_tail ++;

                            // update otus list
                            split_results_list[right_node.node_index].child2_otus_list = merged_otus_list;
                            // update kmldm result
                            split_results_list[right_node.node_index].child2_mldm_res = merged_mldm.best_result;
                            // update the gmm result
                            split_results_list[right_node.node_index].gmm2_res.cluster2_mean = merged_meta_mean;
                            split_results_list[right_node.node_index].gmm2_res.cluster2_cov = merged_meta_cov;
                            split_results_list[right_node.node_index].gmm2_res.cluster2_index_list = merged_samples_list;
                        }
                    }
                // }

                break;
            }
            // else continue
        }

        // if there is not right nodes merged with the left node
        // remove the left node from list
        if (! exist) {
#ifdef DEBUG
            std::cout << "No matched nodes: Left Node is removed from list" << std::endl;
#endif
            left_nodes.pop_back();
        }
    }

#ifdef DEBUG
    std::cout << "+++++++++++++++++++ merge process end~ " << std::endl;
#endif

    return;
}

// choose the childs candidate for merge
void KmLDM::findMatchedChilds(int cur_index, VectorXd& direction, VectorXd& center, std::vector<ChildInfo>& candidates) {
/*#ifdef DEBUG
    std::cout << "left child of cluster ["<< cur_index << "]: child index: "<< split_results_list[cur_index].child1_index << " discard: " << split_results_list[cur_index].child1_discard << std::endl;
    std::cout << "right child of cluster ["<< cur_index << "]: child index: "<< split_results_list[cur_index].child2_index << " discard: " << split_results_list[cur_index].child2_discard << std::endl;
#endif */
    // if left child is not leaf
    if (split_results_list[cur_index].child1_index != -1) {
        findMatchedChilds(split_results_list[cur_index].child1_index, direction, center, candidates);

    // if the left child is leaf and it is not discarded
    } else if (split_results_list[cur_index].child1_index == -1 && ! split_results_list[cur_index].child1_discard) {
        VectorXd left_child_center = split_results_list[cur_index].gmm2_res.cluster1_mean;

#ifdef DEBUG
    /*    double test_zero = (left_child_center - center).norm();
        double test_dot = (left_child_center - center).dot(direction);
        std::cout << "left child center: " << left_child_center << std::endl;
        std::cout << "origin center: " << center << std::endl;
        std::cout << "near zero value: " << test_zero << std::endl;
        std::cout << "direction dot: " << test_dot << std::endl; */
#endif
        if (isZero((left_child_center - center).norm()) || (left_child_center - center).dot(direction) >= 0) {
            ChildInfo left_node;
            left_node.node_index = cur_index;
            left_node.is_left = true;
#ifdef DEBUG
            std::cout << "find left child of Cluster [" << left_node.node_index << "]" << std::endl;
#endif
            candidates.push_back(left_node);
        } else {
#ifdef DEBUG
            std::cout << "left child of Cluster [" << cur_index << "] far from the center, not match" << std::endl;
#endif
        }
    } else {
#ifdef DEBUG
        std::cout << "left child of Cluster [" << cur_index << "] discarded, not match" << std::endl;
#endif
    }

    // if right child is not leaf
    if (split_results_list[cur_index].child2_index != -1) {
        findMatchedChilds(split_results_list[cur_index].child2_index, direction, center, candidates);

    // if the right child is leaf and it is not discarded
    } else if (split_results_list[cur_index].child2_index == -1 && ! split_results_list[cur_index].child2_discard) {
        VectorXd right_child_center = split_results_list[cur_index].gmm2_res.cluster2_mean;

#ifdef DEBUG
/*        double test_zero = (right_child_center - center).norm();
        double test_dot = (right_child_center - center).dot(direction);
        std::cout << "right child center: " << right_child_center << std::endl;
        std::cout << "origin center: " << center << std::endl;
        std::cout << "near zero value: " << test_zero << std::endl;
        std::cout << "direction dot: " << test_dot << std::endl; */
#endif

        if (isZero((right_child_center - center).norm()) || (right_child_center - center).dot(direction) >= 0) {
            ChildInfo right_node;
            right_node.node_index = cur_index;
            right_node.is_left = false;
#ifdef DEBUG
            std::cout << "find right child of Cluster [" << right_node.node_index << "]" << std::endl;
#endif
            candidates.push_back(right_node);
        } else {
#ifdef DEBUG
            std::cout << "right child of Cluster [" << cur_index << "] far from the center, not match" << std::endl;
#endif
        }
    } else {
#ifdef DEBUG
            std::cout << "right child of Cluster [" << cur_index << "] discarded, not match" << std::endl;
#endif
    }

    return;
}

// sort otu size by decending order
bool otu_comp(const OTUSize& a, const OTUSize& b) {
    return a.otu_size > b.otu_size;
}

// filter the zero column from matrix
void KmLDM::filterZeroVarColumn(MatrixXi& x_i, std::vector<int>& feature_list, MatrixXi& x_filter) {
    // convert integer to double
    MatrixXd x_d = x_i.cast <double> ();
    int n = x_d.rows();
    std::vector<OTUSize> feature_info;

    for (int i = 0 ;i < P; ++i) {
        double var_i = computeVar(x_d.col(i));
        int zero_i = zeroNum(x_d.col(i));

        if (isZero(var_i) || (zero_i * 1.0 / n > zero_ratio)) {
            // std::cout << i << " " << zero_i << " " << var_i << std::endl;
            continue;
        }

        double average_i = x_d.col(i).sum() * 1.0 / (n - zero_i);

        OTUSize otu_i;
        otu_i.otu_index = i;
        otu_i.otu_size = average_i;
        otu_i.otu_zero = zero_i;
        feature_info.push_back(otu_i);
    }

    // select some big otu and ensure p < n
    std::sort (feature_info.begin(), feature_info.end(), otu_comp);
    // std::cout << "big otu: " << feature_info[0].otu_index << " " << feature_info[0].otu_size << " " << feature_info[0].otu_zero << std::endl;
    int p = feature_info.size() < ceil(n / 2) ? feature_info.size() : ceil(n / 2);

    // std::cout << "small otu: " << feature_info[p - 1].otu_index << " " << feature_info[p - 1].otu_size << " " << feature_info[p - 1].otu_zero << std::endl;

    // if there are some otus meet the condition
    // otherwise discard the cluster
    if (p > 0) {
        x_filter = MatrixXi::Zero(n, p);
        for (int i = 0; i < p; ++ i) {
            feature_list.push_back(feature_info[i].otu_index);
            x_filter.col(i) = x_i.col(feature_list[i]);
        }

        VectorXi samples_size_list = x_filter.rowwise().sum();
        for (int i = 0; i < n; ++i) {
            // if there is bad sample, add some random counts to avoid
            if (samples_size_list[i] < (1 - zero_ratio) * p) {
                // generate random number to decide adding 1 count to which otu
                VectorXd random_ratio = VectorXd::Random(p);
                // set the prob of add 1 count to zero_ratio
                for (int j = 0; j < random_ratio.size(); ++j) {
                    if (random_ratio[j] > zero_ratio || random_ratio[j] < -zero_ratio) {
                        x_filter(i, j) += 1;
                    }
                }
            }
        }

    }
    return;
}

// read matrix from otu_file
bool KmLDM::loadOTUs(char* otu_file) {
    std::ifstream otu_in;
    bool flag = true;
    try {
        otu_in.open(otu_file);
        int i, j;
        for (i = 0; i < N; ++i) {
            for (j = 0; j < P; ++j) {
                otu_in >> X(i, j);
            }
        }
        otu_in.close();
    } catch (const char* msg) {
        std::cout << msg << std::endl;
        flag = false;
    }

    return flag;
}

// read matrix from meta_file
bool KmLDM::loadMetas(char* meta_file) {
    std::ifstream meta_in;
    bool flag = true;

    try {
        meta_in.open(meta_file);
        int i, j;
        for (i = 0; i < N; ++i) {
            for (j = 0; j < Q; ++j) {
                meta_in >> M(i, j);
            }
        }
        meta_in.close();

        // add small noise to meta data, to avoid the zero variance
        for (int i = 0; i < Q; ++i) {
            double var_i = computeVar(M.col(i));
            // if (isZero(var_i)) {
                // std::cout << i << "-th meta data's variance is zero, add small noise ~" << std::endl;
                M.col(i) += (VectorXd::Random(N, 1) * SMALL_NOISE);
            // }
        }

	// scale M with mean 0 and std 1
	// MatrixXd M_scale = MatrixXd::Zero(N, Q);
	// scale(M, M_scale);
	// M = M_scale;
    } catch (const char* msg) {
        std::cout << msg << std::endl;
        flag = false;
    }

    return flag;
}

// get otus' number
int KmLDM::getOTUsNumber() {
    return P;
}

// read otu table and meta table, load into otus and efs respectively
void KmLDM::loadData(char* otu_file, char* meta_file, char* shape) {
#ifdef VERBOSE
    std::cout << "======================         Load Dataset         =====================" << std::endl;
    std::cout << "The OTU File is: " << otu_file << std::endl;
    std::cout << "The Meta File is: " << meta_file << std::endl;
#endif
    // read N, P, Q from shape file
    std::ifstream shape_in;
    shape_in.open(shape);
    shape_in >> N >> P >> Q;
    // the minimum number of samples of one cluster is setted to P
    // N_min = P;
    // set the minimum cluster'size >= 2 * noise size
    /*if (N_min < 2 * noise_size) {
        N_min = noise_size * 2;
    } */
    shape_in.close();
#ifdef VERBOSE
    std::cout << "N = " << N << " P = " << P << " Q = " << Q << std::endl;
#endif
    X = MatrixXi::Zero(N, P);
    M = MatrixXd::Zero(N, Q);

#ifdef VERBOSE
    std::cout << "OTU Counts Reading..." << std::endl;
#endif

    // read otu table (N * P) from otu file
    if (loadOTUs(otu_file)) {
#ifdef VERBOSE
        std::cout << "OTU Counts Loaded ~" << std::endl;
#endif
    } else {
        std::cout << "OTU Counts Read ERROR !" << std::endl;
        exit(0);
    }

#ifdef VERBOSE
    std::cout << "Meta Values Reading..." << std::endl;
#endif

    // read meta table (N * Q) from meta file
    if (loadMetas(meta_file)) {
#ifdef VERBOSE
        std::cout << "Meta Values Loaded ~" << std::endl;
#endif
    } else {
        std::cout << "Meta Values Read ERROR !" << std::endl;
        exit(0);
    }

    return;
}

// output result into file
void exportCluster(int number, mLDMResult& result, std::vector<int>& otus, VectorXd& meta_mean, MatrixXd& meta_cov, std::vector<int> sample_list, int sample_num, char* result_dir_name) {
    // the cluster data files
    std::string buffer = std::to_string(number);
    std::string result_dir = result_dir_name;
    std::string cluster_theta_file = result_dir + "/Theta_" + buffer;
    std::string cluster_b_file = result_dir + "/B_" + buffer;
    std::string cluster_otu_list_file = result_dir + "/OTU_Index_" + buffer;
    std::string cluster_meta_mean_file = result_dir + "/Meta_Mean_" + buffer;
    std::string cluster_meta_cov_file = result_dir + "/Meta_Cov_" + buffer;
    std::string cluster_sample_list_file = result_dir + "/Sample_Index_" + buffer;
    std::string cluster_sample_num_file = result_dir + "/Sample_Num_" + buffer;
    std::string cluster_otu_otu_file = result_dir + "/OTU_OTU_Association_" + buffer;
    std::string cluster_meta_otu_file = result_dir + "/Meta_OTU_Association_" + buffer;
    //std::cout << "===========================================================================" << "\n";
    //std::cout << "Cluster [" << number << "] (samples num: " << sample_num << ")"<< "\n";
    //std::cout << "===========================================================================" << "\n";
    // fout << "otus list: " << "\n"; 

    // export sample number
    std::ofstream num_out(cluster_sample_num_file);
    num_out << sample_num;
    num_out.close();

    int i,j;
    int p = otus.size();
    int q = result.B.rows();
    // export otu list
    std::ofstream otu_out(cluster_otu_list_file);
    for (i = 0; i < p; ++ i) {
        otu_out << otus[i];
        if (i != p - 1) {
            otu_out << " ";
        }
    }
    otu_out.close();

    // export sample list 
    //std::cout << "export sample_list" << std::endl;
    std::ofstream sample_out(cluster_sample_list_file);
    for (i = 0; i < sample_num; ++ i) {
        sample_out << sample_list[i];
        if (i != sample_num - 1) {
            sample_out << " ";
        }
        // std::cout << sample_list[i] << " ";
    }
    sample_out.close();

    // fout << "B" << "\n";
    // export B 
    std::ofstream b_out(cluster_b_file);
    for (i = 0; i < q; ++i) {
        for (j = 0; j < p; ++j) {
            b_out << result.B(i, j);
            if (j != p - 1) {
                b_out << " ";
            }
        }
        if (i != q - 1) {
            b_out << "\n";
        }
    }
    b_out.close();

    // export meta-otu associations
    std::ofstream meta_otu_out(cluster_meta_otu_file);
    for (i = 0; i < q; ++i) {
        for (j = 0; j < p; ++j) {
            meta_otu_out << result.B(i, j);
            if (j != p - 1) {
                meta_otu_out << " ";
            }
        }
        if (i != q - 1) {
            meta_otu_out << "\n";
        }
    }
    meta_otu_out.close();

    std::ofstream theta_out(cluster_theta_file);
    for (i = 0; i < p; ++i) {
        for (j = 0; j < p; ++j) {
            theta_out << result.Theta(i, j);
            if (j != p - 1) {
                theta_out << " ";
            }
        }
        if (i != p - 1) {
            theta_out << "\n";
        }
    }
    theta_out.close();

    // fout << "Theta" << "\n";
    // export Theta 
    // Theta(i,j) = - Theta(i,j) / sqrt(Theta(i,i) * Theta(j,j))
    std::ofstream otu_otu_out(cluster_otu_otu_file);
    for (i = 0; i < p; ++i) {
        for (j = 0; j < p; ++j) {
            if (i != j) {
                otu_otu_out << - result.Theta(i, j) / sqrt(result.Theta(i,i) * result.Theta(j,j));
            } else {
                otu_otu_out << 0.0;
            }

            if (j != p - 1) {
                otu_otu_out << " ";
            }
        }
        if (i != p - 1) {
            otu_otu_out << "\n";
        }
    }
    otu_otu_out.close();

    // fout << "Mean of EFs" << "\n";
    // export ef mean 
    std::ofstream mean_out(cluster_meta_mean_file);
    for (i = 0; i < q; ++i) {
        mean_out << meta_mean[i];
        if (i != q - 1) {
            mean_out << " ";
        }
    }
    mean_out.close();

    // fout << "Cov of EFs" << "\n";
    // export ef cov
    std::ofstream cov_out(cluster_meta_cov_file);
    for (i = 0; i < q; ++i) {
        for (j = 0; j < q; ++j) {
            cov_out << meta_cov(i, j);
            if (j != q - 1) {
                cov_out << " ";
            }
        }
        if (i != q - 1) {
            cov_out << "\n";
        }
    }
    cov_out.close();

    return;
}

// export hierarchical results 
void KmLDM::exportHierarchicalNodes(int cluster_number, int node_index, bool is_left, std::string result_dir_name) {

    SplitResult cur_node = split_results_list[node_index];
    mLDMResult cur_result;
    std::vector<int> cur_otus;
    std::vector<int> cur_sample_list;
    VectorXd cur_meta_mean;
    MatrixXd cur_meta_cov;
    int cur_sample_num;
    std::vector<int> cur_old_list;
    std::vector<ChildInfo> cur_merge_list;

    if (is_left) {
        cur_result = cur_node.child1_mldm_res;
        cur_otus = cur_node.child1_otus_list;
        cur_meta_mean = cur_node.gmm2_res.cluster1_mean;
        cur_meta_cov = cur_node.gmm2_res.cluster1_cov;
        cur_sample_list = cur_node.gmm2_res.cluster1_index_list;
        cur_sample_num = cur_sample_list.size();
        cur_old_list = cur_node.child1_old_list;
        cur_merge_list = cur_node.child1_merge_list;
    } else {
        cur_result = cur_node.child2_mldm_res;
        cur_otus = cur_node.child2_otus_list;
        cur_meta_mean = cur_node.gmm2_res.cluster2_mean;
        cur_meta_cov = cur_node.gmm2_res.cluster2_cov;
        cur_sample_list = cur_node.gmm2_res.cluster2_index_list;
        cur_sample_num = cur_sample_list.size();
        cur_old_list = cur_node.child2_old_list;
        cur_merge_list = cur_node.child2_merge_list;
    }

    // data dir 
    std::string cur_dir = result_dir_name + "/merged/";
    safe_mkdir(cur_dir);

    // export current node 
    exportCluster(cluster_number, cur_result, cur_otus, cur_meta_mean, cur_meta_cov, cur_sample_list, cur_sample_num, const_cast<char*>(cur_dir.c_str()));

    int layer = cur_merge_list.size();
    for (int i = layer - 1; i > -1; --i) {
        // mkdir 
        std::string origin_dir = result_dir_name + "/original/";
        safe_mkdir(origin_dir);
        origin_dir = origin_dir + "/merged/";
        safe_mkdir(origin_dir);
        std::string add_dir = result_dir_name + "/added/"; 
        safe_mkdir(add_dir);

        // export original node 
        int origin_index = cur_old_list[i];
        OldCluster origin_node = old_results_list[origin_index];
        VectorXd origin_meta_mean = origin_node.meta_mean;
        MatrixXd origin_meta_cov = origin_node.meta_cov;
        std::vector<int> origin_sample_index = origin_node.sample_index;
        int origin_sample_num = origin_sample_index.size();
        std::vector<int> origin_otus = origin_node.otu_index;
        mLDMResult origin_result = origin_node.mldm_res;
        exportCluster(1, origin_result, origin_otus, origin_meta_mean, origin_meta_cov, origin_sample_index, origin_sample_num, const_cast<char*>(origin_dir.c_str()));
        
        // export added node 
        ChildInfo another_node = cur_merge_list[i];
        exportHierarchicalNodes(1, another_node.node_index, another_node.is_left, add_dir);

        // update dir 
        result_dir_name = result_dir_name + "/original/";
    }

    return;
}

// export the estimated association networks into file
void KmLDM::exportResult(char* result_dir) {
    // push parent nodes into stack, to find leaf nodes
    std::vector<int> stack;
    stack.push_back(0);
    int index = -1;
    SplitResult cur_node;
    mLDMResult leaf_result;
    std::vector<int> leaf_otus;
    std::vector<int> leaf_sample_list;
    VectorXd leaf_meta_mean;
    MatrixXd leaf_meta_cov;
    int leaf_sample_num;
    int cluster_number = 1;
    std::string result_dir_name = result_dir;
    std::string hierarchical_dir = result_dir_name + "/hierarchical_results/";
    // mkdir 
    safe_mkdir(hierarchical_dir);

#ifdef VERBOSE
    std::cout << "======================        Export Result        ======================" << std::endl;
#endif

    try{
        // std::ofstream fout;
        // fout.open(result_file);

        while(! stack.empty()) {
            index = stack.back();
            // std::cout << "pop cluster index: " << index << std::endl;
            stack.pop_back();
            // get current node
            cur_node = split_results_list[index];
            // if child1 is not leaf, push into stack
            if (cur_node.child1_index != -1) {
                stack.push_back(cur_node.child1_index);
            }
            // if child2 is not leaf, push into stack
            if (cur_node.child2_index != -1) {
                stack.push_back(cur_node.child2_index);
            }
            // if child1 is leaf
            if (cur_node.child1_index == -1 && ! cur_node.child1_discard) {
                // std::cout << "export child1 " << std::endl;
                leaf_result = cur_node.child1_mldm_res;
                leaf_otus = cur_node.child1_otus_list;
                leaf_meta_mean = cur_node.gmm2_res.cluster1_mean;
                leaf_meta_cov = cur_node.gmm2_res.cluster1_cov;
                leaf_sample_list = cur_node.gmm2_res.cluster1_index_list;
                leaf_sample_num = leaf_sample_list.size();

                exportCluster(cluster_number, leaf_result, leaf_otus, leaf_meta_mean, leaf_meta_cov, leaf_sample_list, leaf_sample_num, result_dir);
                // export hierachical results  
                std::string number_str = std::to_string(cluster_number);
                std::string leaf_result_dir = hierarchical_dir + "/cluster-" + number_str + "/";
                // mkdir 
                safe_mkdir(leaf_result_dir);

                exportHierarchicalNodes(1, index, true, leaf_result_dir);

                cluster_number += 1;
            }
            // if child2 is leaf
            if (cur_node.child2_index == -1 && ! cur_node.child2_discard) {
                // std::cout << "export child2" << std::endl;
                leaf_result = cur_node.child2_mldm_res;
                leaf_otus = cur_node.child2_otus_list;
                leaf_meta_mean = cur_node.gmm2_res.cluster2_mean;
                leaf_meta_cov = cur_node.gmm2_res.cluster2_cov;
                leaf_sample_list = cur_node.gmm2_res.cluster2_index_list;
                leaf_sample_num = leaf_sample_list.size();

                exportCluster(cluster_number, leaf_result, leaf_otus, leaf_meta_mean, leaf_meta_cov, leaf_sample_list, leaf_sample_num, result_dir);
                // export hierachical results  
                std::string number_str = std::to_string(cluster_number);
                std::string leaf_result_dir = hierarchical_dir + "cluster-" + number_str + "/";
                // mkdir
                safe_mkdir(leaf_result_dir);

                exportHierarchicalNodes(1, index, false, leaf_result_dir);
                

                cluster_number += 1;
            }
        }

        // fout.close();

        // export cluster number
        std::string cluster_number_file = result_dir;
        cluster_number_file += "/Cluster_Number";
        std::ofstream number_out(cluster_number_file);
        number_out << cluster_number - 1;
        number_out.close();

#ifdef VERBOSE
        std::cout << "export success ~" << std::endl;
#endif

    } catch (const char* msg) {
        std::cout << msg << std::endl;
        exit(0);
    }

}

