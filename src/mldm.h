#ifndef MLDM_H
#define MLDM_H

#include "utils.h"
#include "proximalQusiNewton.h"

#include "./lbfgs/lbfgs.h"

#include "./QUIC/QUIC.h"

// record the estimated association networks of the mldm
struct mLDMResult {
    // lambda1 and lambda2
    double lambda1;
    double lambda2;
    // store the associations among OTUs
    MatrixXd Theta;
    // store the associations between OTUs and metas
    MatrixXd B;
    // the basis value
    VectorXd B0;
    // the EBIC value
    double EBIC;
    // the obejctive value
    double obj;
    // the result exist
    bool exists;
};

class mLDM {
private:
    // sample num
    int N;
    // otu num
    int P;
    // meta num
    int Q;
    // otu table
    MatrixXi X;
    // meta table
    MatrixXd M;
    // delta of two objective function values
    double threshold;
    // max iterations of mLDM
    int max_iterations;

    // lbfgs
    // max iterations Z
    int max_iteration_Z;
    // num of vector to appproximate the hessian
    int approx_num_Z;
    // max times of line search
    int max_linesearch_Z;
    // eplison of Z to check convergence
    double eplison_Z;

    // proximal qusi-newton
    // num of vector to appproximate the hessian
    int approx_num_B;
    // max times of line search
    int max_linesearch_B;
    // max iterations of qusi-newton
    int max_iteration_B;
    // stop threshold of qusi-newton
    double threshold_B;
    // threshold of strong wolfe condition
    double delta1_threshold_B;
    double delta2_threshold_B;
    // check the singularity of the matrix
    double sy_threshold_B;
    // max times of coordinate descent of B
    int max_iteration_B_coor;
    // stop threshold of coordinate
    double threshold_B_coor;

    // quic
    int max_iteration_Theta;
    double threshold_Theta;
    int msg_Theta;


    double ratio1;
    double ratio2;

    int model_selection_num;
    // if show debug info or not
    bool verbose;

public:
    mLDMResult best_result;

public:
    // get the elements within lower triangle of matrix a
    void getStrictLowerValues(MatrixXd& a, std::vector<double> values_list);

    void setLambdaRatio(double r1, double r2);
    void setModelSelectionNum(int model_num);
    // void enableRemoveSample(bool flag);

    // set the parameters of lbfgs
    void setParametersLBFGS(int ap_num_z, int max_line_z, int max_iters_z, double eplison);

    // set the parameters of QUIC
    void setParametersQUIC(int max_iters, double thres, int msg);

    // set the parameters of qusi-newton
    void setParametersQusiNewton(int ap_num_B, int max_line_B, int max_iters_B, double thres_B, double delta1_thres_B, double delta2_thres_B, double sy_thres_B, int max_iter_B_coor, double thres_B_coor);

    // lognormal-dirichlet-multinomial model
    void LognormalDirichletMultinomial(MatrixXd& X, MatrixXd& M, MatrixXd& B, VectorXd& B0, MatrixXd& Theta, MatrixXd& Z_init, bool verbose, bool loopFlag, mLDMResult& per_solution);

    // initialize parameters
    mLDM(MatrixXi x, MatrixXd m, int n, int p, int q, int max_iters, double thres, bool verb/*, int noise_size*/);

    // estimate the association networks via metagenomic lognormal dirichlet multinomial model
    void work();

private:
    // compute the value of obejctive function
    double computeObjf(MatrixXd& X, MatrixXd& M, MatrixXd& B, MatrixXd& Theta, MatrixXd& Z, double& lambda1, double& lambda2);
    // Compute the EBIC
    double computeEBIC(double& objNew, MatrixXd& B, MatrixXd& Theta, double& lambda1, double& lambda2);
};

// define the objective function for proximal qusi-newton
class PQNPoint: public Point {
    int index;
    int P;
    int Q;
    int N;
    MatrixXd BZ;
    MatrixXd X;
    MatrixXd M;

public:
    PQNPoint(int& p, int& q, int& n, MatrixXd& xm, MatrixXd& mm);
    virtual bool Evaluate(const VectorXd& w, double& obj, VectorXd& grad) const;
    void setIndex(int& i);
    void setBZ(MatrixXd& bz);
    void setX(MatrixXd& xm);
    void setM(MatrixXd& mm);
    void setN(int& N);
};

#endif
