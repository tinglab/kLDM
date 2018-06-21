#include "gmm2.h"

//#define DEBUG

// K-Means
void initializationCenters(MatrixXd& efs, GMM2Result& result, double& delta_kmeans, int& Q,
                           int& max_iterations_kmeans, int& repeat_times_kmeans) {
    // initialize two centers
    int samples_num = efs.rows();
    // record the best mean
    VectorXd cluster1_mean = VectorXd::Zero(Q);
    VectorXd cluster2_mean = VectorXd::Zero(Q);
    double cost_best = 1e+20;
    srand(unsigned(time(0)));

    // set parameters
    double cost_now;
    double cost_new;
    int cluster1_size;
    VectorXd center1_sum;
    VectorXd center2_sum;
    int t = 0;
    int step;
    int index1, index2;

    while (t < repeat_times_kmeans) {
        t += 1;

        cost_now = 0.0;
        cost_new = 1.0;
        cluster1_size = 0;
        step = 0;
        center1_sum = VectorXd::Zero(Q);
        center2_sum = VectorXd::Zero(Q);

        index1 = floor(rand() % samples_num);
        index2 = floor(rand() % samples_num);
        while (index1 == index2) {
            index2 = floor(rand() % samples_num);
        }
	    // std::cout << "index: " << index1 << " " << index2 << std::endl;
        cluster1_mean = efs.row(index1).transpose();
        cluster2_mean = efs.row(index2).transpose();

        // when the variation <= delta_kmeans or the number of iterations > max_iterations_kmeans, the kmeans stop
        while (fabs(cost_new - cost_now) > delta_kmeans) {
            cost_now = cost_new; cost_new = 0.0;
            cluster1_size = 0;
            center1_sum = VectorXd::Zero(Q);
            center2_sum = VectorXd::Zero(Q);

            for (int i = 0; i < samples_num; ++ i) {
                double distance1 = (efs.row(i).transpose() - cluster1_mean).squaredNorm();
                double distance2 = (efs.row(i).transpose() - cluster2_mean).squaredNorm();

                if (distance2 > distance1) {
                    center1_sum += efs.row(i).transpose();
                    cost_new += distance1;
                    cluster1_size += 1;
                } else {
                    center2_sum += efs.row(i).transpose();
                    cost_new += distance2;
                }
            }

            // update two centers
            if (cluster1_size < 1) {
                cluster1_mean = efs.row(floor(rand() % samples_num)).transpose();
                cluster2_mean = center2_sum / (samples_num - cluster1_size);
            } else if (cluster1_size == samples_num) {
                cluster1_mean = center1_sum / cluster1_size;
                cluster2_mean = efs.row(floor(rand() % samples_num));
            } else {
                cluster1_mean = center1_sum / cluster1_size;
                cluster2_mean = center2_sum / (samples_num - cluster1_size);
            }
            // std::cout << "cluster1 size: " << cluster1_size << std::endl;

            step ++;
            if (step > max_iterations_kmeans) {
                break;
            }
        }
	// std::cout << "new kmeans cost: " << cost_new << std::endl;
        if (cost_new < cost_best) {
	    // std::cout << "find better centers" << std::endl;
            cost_best = cost_new;
            result.cluster1_mean = cluster1_mean;
            result.cluster2_mean = cluster2_mean;
        }
    }
    //std::cout << "kmeans steps: " << step << std::endl;
    //std::cout << "kmeans cost: " << cost_new << std::endl;
}

// 2 Gaussian Mixture Model
void GMM2(MatrixXd& efs, GMM2Result& result, double& delta_2gmm, double& delta_kmeans,
    int& Q, int& max_iterations_2gmm, int& max_iterations_kmeans, int& repeat_times_kmeans) {
    // initialize the mean and covariance of 2 centers via K-Means
    initializationCenters(efs, result, delta_kmeans, Q, max_iterations_kmeans, repeat_times_kmeans);
    // set covariance matrixs to identity matrix
    result.cluster1_cov = MatrixXd::Identity(Q,Q);
    result.cluster2_cov = MatrixXd::Identity(Q,Q);
    // get the inverse of the covariance matrix
    MatrixXd cluster1_icov = result.cluster1_cov.inverse();
    MatrixXd cluster2_icov = result.cluster2_cov.inverse();
    int samples_num = efs.rows();
    // record the propability belongs to per cluster for every point
    VectorXd points_cluster1_prob = VectorXd::Zero(samples_num);
    VectorXd points_cluster2_prob = VectorXd::Zero(samples_num);
    // record the cluster propability
    double cluster1_prob = 0.5;
    double obj_now = 0;
    double obj_new = 1;
    int step = 0;

    while (fabs(obj_new - obj_now) > delta_2gmm) {
        // std::cout << "mean1: " << result.cluster1_mean << std::endl;
        // std::cout << "mean2: " << result.cluster2_mean << std::endl;
        // std::cout << "cov1: " << result.cluster1_cov << std::endl;
        // std::cout << "cov2: " << result.cluster2_cov << std::endl;

        obj_now = obj_new; obj_new = 0;

        // - 1/2 * ln(cov_1^-1)
        double icov1_det_ln = - 0.5 * computeLnDet(result.cluster1_cov);
        // - 1/2 * ln(cov_2^-1)
        double icov2_det_ln = - 0.5 * computeLnDet(result.cluster2_cov);
        double likeli_per = 0;
        // std::cout << "ln det icov " << icov1_det_ln << " " << icov2_det_ln << std::endl;
        // std::cout << "cluster1 prob: " << cluster1_prob << std::endl;
        for (int i = 0; i < samples_num; ++ i) {
            // (x - u_1) and (x - u_2)
            VectorXd distance_vector1 = efs.row(i).transpose() - result.cluster1_mean;
            VectorXd distance_vector2 = efs.row(i).transpose() - result.cluster2_mean;
            // ln(p(x | z_1)) + ln(p(z_1)) and ln(p(x | z_2)) + ln(p(z_2))
            // -0.5 * ln(cov_1^-1) - 0.5 * (x - u_1)^T * cov_1^-1 * (x - u_1) + ln(p(z_1))
            points_cluster1_prob(i) = - 0.5 * (distance_vector1.transpose() * cluster1_icov * distance_vector1).array()(0, 0) + icov1_det_ln + log(cluster1_prob);
            // -0.5 * ln(cov_2^-1) - 0.5 * (x - u_2)^T * cov_2^-1 * (x - u_2) + ln(p(z_2))
            points_cluster2_prob(i) = - 0.5 * (distance_vector2.transpose() * cluster2_icov * distance_vector2).array()(0, 0) + icov2_det_ln + log(1 - cluster1_prob);
            // compute the denominator
            if (points_cluster1_prob(i) > points_cluster2_prob(i)) {
                likeli_per = points_cluster1_prob(i) + log(exp(points_cluster2_prob(i) - points_cluster1_prob(i)) + 1);
                points_cluster1_prob(i) = exp(points_cluster1_prob(i) - likeli_per);
                points_cluster2_prob(i) = 1 - points_cluster1_prob(i);
            } else {
                likeli_per = points_cluster2_prob(i) + log(exp(points_cluster1_prob(i) - points_cluster2_prob(i)) + 1);
                points_cluster2_prob(i) = exp(points_cluster2_prob(i) - likeli_per);
                points_cluster1_prob(i) = 1 - points_cluster2_prob(i);
            }
            // update the obj value
            obj_new += likeli_per - log(2 * PI) * Q/2;
        }

        // update p(z_1) and p(z_2)
        double cluster1_prob_sum = points_cluster1_prob.array().sum();
        double cluster2_prob_sum = points_cluster2_prob.array().sum();
        cluster1_prob = cluster1_prob_sum / samples_num;

        // update mean1 and mean2
        MatrixXd temp1 = MatrixXd::Zero(samples_num, Q);
        computeVecMatMult(points_cluster1_prob, efs, temp1);
        result.cluster1_mean = temp1.colwise().sum().transpose() / cluster1_prob_sum;

        computeVecMatMult(points_cluster2_prob, efs, temp1);
        result.cluster2_mean = temp1.colwise().sum().transpose() / cluster2_prob_sum;

        // update cov1 and cov2
        MatrixXd temp2 = MatrixXd::Zero(Q, samples_num);
        MatrixXd mean1_matrix = efs.transpose().colwise() - result.cluster1_mean;
        computeVecMatMultTranspose(points_cluster1_prob, mean1_matrix, temp2);
        result.cluster1_cov = temp2 * mean1_matrix.transpose() / cluster1_prob_sum;

        MatrixXd mean2_matrix = efs.transpose().colwise() - result.cluster2_mean;
        computeVecMatMultTranspose(points_cluster2_prob, mean2_matrix, temp2);
        result.cluster2_cov = temp2 * mean2_matrix.transpose() / cluster2_prob_sum;
        // update icov1 and icov2
        cluster1_icov = result.cluster1_cov.inverse();
        cluster2_icov = result.cluster2_cov.inverse();

        obj_new /= samples_num;
        step ++;
        // std::cout << "delta gmm: " << obj_new - obj_now << std::endl;
        if (step > max_iterations_2gmm) {
            break;
        }
    }
    // std::cout << "cluster1 prob: " << cluster1_prob << std::endl;
    // std::cout << "after_clustering" << std::endl;
    // record the index belongs to two clusters
    // compute the final mean and cov (not use the probability)
    int cluster1_size = 0, cluster2_size = 0;
    for (int i = 0; i < samples_num; ++ i) {
        if (points_cluster1_prob(i) > 0.5) {
            result.cluster1_index_list.push_back(i);
            cluster1_size += 1;
        } else {
            result.cluster2_index_list.push_back(i);
            cluster2_size += 1;
        }
    }
    std::cout << "c1 size: " << cluster1_size << std::endl;
    std::cout << "c2 size: " << cluster2_size << std::endl;

    MatrixXd cluster1_efs = MatrixXd::Zero(cluster1_size, Q);
    MatrixXd cluster2_efs = MatrixXd::Zero(cluster2_size, Q);
    for (int i = 0; i < cluster1_size; ++i) {
        cluster1_efs.row(i) = efs.row(result.cluster1_index_list[i]);
    }
    for (int i = 0; i < cluster2_size; ++i) {
        cluster2_efs.row(i) = efs.row(result.cluster2_index_list[i]);
    }
    // update mean1 and cov1
    result.cluster1_mean = cluster1_efs.colwise().mean().transpose();
    MatrixXd cluster1_meta_cov = MatrixXd::Zero(Q, Q);
    computeCov(cluster1_efs, cluster1_meta_cov);
    result.cluster1_cov = cluster1_meta_cov;
    // update mean2 and cov2
    result.cluster2_mean = cluster2_efs.colwise().mean().transpose();
    MatrixXd cluster2_meta_cov = MatrixXd::Zero(Q, Q);
    computeCov(cluster2_efs, cluster2_meta_cov);
    result.cluster2_cov = cluster2_meta_cov;

#ifdef DEBUG
    std::cout << "GMM2 result: " << result.cluster1_index_list.size() << " and " << result.cluster2_index_list.size() << std::endl;
#endif

    return;
}

