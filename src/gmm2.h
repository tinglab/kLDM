#ifndef GMM2_H
#define GMM2_H

#include "utils.h"

// record the result of 2-Gaussian Mixtures Model
struct GMM2Result {
    // store the real index of datapoints belongs to cluster1
    std::vector<int> cluster1_index_list;
    // store the real index of datapoints belongs to cluster2
    std::vector<int> cluster2_index_list;
    // store the mean vector of the cluster1
    VectorXd cluster1_mean;
    // store the mean vector of the cluster2
    VectorXd cluster2_mean;
    // store the covariance matrix of the cluster1
    MatrixXd cluster1_cov;
    // store the covariance matrix of the cluster2
    MatrixXd cluster2_cov;
};

// K-Means
void initializationCenters(MatrixXd& efs, GMM2Result& result, double& delta_kmeans,
    int& Q, int& max_iterations_kmeans, int& repeat_times_kmeans);
// 2 Gaussian Mixture Model
void GMM2(MatrixXd& efs, GMM2Result& result, double& delta_2gmm, double& delta_kmeans,
    int& Q, int& max_iterations_2gmm, int& max_iterations_kmeans, int& repeat_times_kmeans);

#endif
