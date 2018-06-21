#include "mldm.h"

// compute alphai = exp(u_i + z_i)
void getAlpha(const VectorXd& mui, const VectorXd& zi, VectorXd& alpha) {
    alpha = (mui + zi).array().exp();
    return;
}

// compute the lgamma part of objective funciton for zi
void computeAlphaVector(VectorXd& ax, VectorXd& a, double& obj) {
    ArrayXd a_sum(2);
    a_sum[0] = ax.sum();
    a_sum[1] = a.sum();
    a_sum = a_sum.lgamma();

    obj = - (ax.array().lgamma() - a.array().lgamma()).sum() + (a_sum[0] - a_sum[1]);
    return;
}

// compute the digamma part of objective funciton for zi
void computeAlphaDerVector(VectorXd& ax, VectorXd& a, VectorXd& der) {
    ArrayXd a_sum(2);
    a_sum[0] = ax.sum();
    a_sum[1] = a.sum();
    a_sum = a_sum.digamma();

    der = - (ax.array().digamma() - a.array().digamma()).array() + (a_sum[0] - a_sum[1]);
    return;
}

void computeAlphaMatrix(MatrixXd& AlphaMX, MatrixXd& AlphaM, double& obj) {
    int n = AlphaMX.cols();
    int p = AlphaMX.rows();
    double sum1 = 0;
    double sum2 = 0;

    #pragma omp parallel for reduction(+:sum1,sum2)
    for (int i = 0; i < n; ++i) {
        VectorXd alphax_i(p);
        alphax_i = AlphaMX.col(i);
        VectorXd alpha_i(p);
        alpha_i = AlphaM.col(i);
        sum1 += - (alphax_i.array().lgamma() - alpha_i.array().lgamma()).sum();
        ArrayXd two_sum(2);
        two_sum << alphax_i.sum(), alpha_i.sum();
        two_sum = two_sum.lgamma();
        // sum2 += (lgamma(alphax_i.sum()) - lgamma(alpha_i.sum()));
        sum2 += (two_sum[0] - two_sum[1]);
    }
    obj = sum1 + sum2;

    return;
}

void computeAlphaDerMatrix(MatrixXd& AlphaMX, MatrixXd& AlphaM, MatrixXd& der) {
    int n = AlphaMX.cols();
    int p = AlphaMX.rows();

    #pragma omp parallel for
    for(int i = 0 ; i < n; ++i) {
        VectorXd alphax_i(p);
        alphax_i = AlphaMX.col(i);
        VectorXd alpha_i(p);
        alpha_i = AlphaM.col(i);
        ArrayXd two_sum(2);
        two_sum << alphax_i.sum(), alpha_i.sum();
        two_sum = two_sum.digamma();
        der.col(i) = - (alphax_i.array().digamma() - alpha_i.array().digamma()).array() + (two_sum[0] - two_sum[1]);
    }

    return;
}

struct ZiData {
    VectorXd Mui;
    VectorXd Xi;
    VectorXd B0;
    MatrixXd Theta;
    int N;
};

// compute obj and gradient
static lbfgsfloatval_t lbfgs_evaluate(void* instance, const lbfgsfloatval_t* parameters, lbfgsfloatval_t* gradient, const int P, const lbfgsfloatval_t step) {
    ZiData* data = (ZiData*) instance;
    VectorXd z = Eigen::Map<const VectorXd>(parameters, P);
    VectorXd alphai;
    getAlpha(data->Mui, z, alphai);
    VectorXd alphaix = alphai + data->Xi;
    double obj;
    computeAlphaVector(alphaix, alphai, obj);
    obj += ((z - data->B0).transpose() * data->Theta * (z - data->B0))(0, 0) / 2;
    // obj 
    obj = obj / data->N;
    // grad
    if (gradient != NULL) {
        VectorXd der;
        computeAlphaDerVector(alphaix, alphai, der);
        VectorXd grad = (der.array() * alphai.array() + (data->Theta * (z - data->B0)).array()) / data->N;
        for (int i = 0; i < P; ++i) {
            gradient[i] = grad[i];
        }

    }
    return obj;
}

// implement the obj and grad functions of B for proximal qusi-newton method
PQNPoint::PQNPoint(int& p, int& q, int& n, MatrixXd& xm, MatrixXd& mm) {
    P = p;
    Q = q;
    N = n;
    X = xm;
    M = mm;
}

// set the index
void PQNPoint::setIndex(int& i) {
    this->index = i;
}
// set the matrix BZ
void PQNPoint::setBZ(MatrixXd& bz) {
    this->BZ = bz;
}

void PQNPoint::setX(MatrixXd& xm) {
    this->X = xm;
}

void PQNPoint::setM(MatrixXd& mm) {
    this->M = mm;
}

void PQNPoint::setN(int& N) {
    this->N = N;
}

bool PQNPoint::Evaluate(const VectorXd& w, double& obj, VectorXd& grad) const {
    MatrixXd AlphaM(P, N); 
    MatrixXd AlphaMX(P, N);
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        AlphaM.col(i) = (BZ.col(i).array() + M(i, index) * w.array()).array().exp();
        AlphaMX.col(i) = AlphaM.col(i) + X.row(i).transpose();
    }

    computeAlphaMatrix(AlphaMX, AlphaM, obj);
    obj /= N;

    MatrixXd der(P, N);
    computeAlphaDerMatrix(AlphaMX, AlphaM, der);
    VectorXd mi(N);
    mi = M.col(index);
    #pragma omp parallel for
    for (int i = 0; i < P; ++i) {
        VectorXd temp(N);
        temp = der.row(i).array() * AlphaM.row(i).array();
        grad[i] = (temp.transpose() * mi)(0,0) / N;
    }

    return true;
}

// init dataset
mLDM::mLDM(MatrixXi x, MatrixXd m, int n, int p, int q, int max_iters, double thres, bool verb/*, int noise_size*/) {
    // sampel num
    N = n;
    // otu num
    P = p;
    // meta num
    Q = q;
    // otu table
    X = x;
    // meta table
    M = m;
    // max iterations and delta of two objectives
    max_iterations = max_iters;
    threshold = thres;
    // set the ratio to choose the lambda list
    ratio1 = 0.6;
    ratio2 = 0.9;
    model_selection_num = 4;
    // lbfgs parameters
    approx_num_Z = 20;
    max_linesearch_Z = 30;
    eplison_Z = 1e-4;
    max_iteration_Z = 200;
    // proximal qusi-newton method
    approx_num_B = 20;
    max_linesearch_B = 30;
    max_iteration_B = 200;
    threshold_B = 1e-4;
    delta1_threshold_B = 1e-4;
    delta2_threshold_B = 0.9;
    sy_threshold_B = 1e-6;
    max_iteration_B_coor = 20;
    threshold_B_coor = 1e-6;
    // debug info
    verbose = verb;
    // quic parameters
    max_iteration_Theta = 200;
    threshold_Theta = 1e-4;
    msg_Theta = 0;
}

// set the quantiles when choose lambda1 and lambda2
void mLDM::setLambdaRatio(double r1, double r2) {
    ratio1 = r1;
    ratio2 = r2;
}

// set the num of model selection
void mLDM::setModelSelectionNum(int model_num) {
    model_selection_num = model_num;
}

// set the parameters of lbfgs
void mLDM::setParametersLBFGS(int ap_num_z, int max_line_z, int max_iters_z, double eplison) {
    approx_num_Z = ap_num_z;
    max_linesearch_Z = max_line_z;
    eplison_Z = eplison;
    max_iteration_Z = max_iters_z;
}

// set the parameters of qusi-newton
void mLDM::setParametersQusiNewton(int ap_num_B, int max_line_B, int max_iters_B, double thres_B, double delta1_thres_B, double delta2_thres_B, double sy_thres_B, int max_iter_B_coor, double thres_B_coor) {
    // num of vector to appproximate the hessian
    approx_num_B = ap_num_B;
    // max times of line search
    max_linesearch_B = max_line_B;
    // max iterations of qusi-newton
    max_iteration_B = max_iters_B;
    // stop threshold of qusi-newton
    threshold_B = thres_B;
    // threshold of strong wolfe condition
    delta1_threshold_B = delta1_thres_B;
    delta2_threshold_B = delta2_thres_B;
    // check the singularity of the matrix
    sy_threshold_B = sy_thres_B;
    // max times of coordinate descent of B
    max_iteration_B_coor = max_iter_B_coor;
    // stop threshold of coordinate
    threshold_B_coor = thres_B_coor;
}

// get the elements within lower triangle of matrix a
void getSymmStrictLowerValues(MatrixXd& a, std::vector<double>& values_list) {
    int p = a.cols();
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < i; ++j) {
            values_list.push_back(a(j, i));
        }
    }

    return;
}

// push values of matrix into vector
void getMatrixValuesList(MatrixXd& a, std::vector<double>& values_list) {
    int n = a.rows();
    int p = a.cols();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            values_list.push_back(a(i, j));
        }
    }

    return;
}

// Compute the current objective function for f
double mLDM::computeObjf(MatrixXd& X, MatrixXd& M, MatrixXd& B, MatrixXd& Theta, MatrixXd& Z, double& lambda1, double& lambda2) {
    VectorXd B0 = Z.colwise().mean();
    MatrixXd Zcenter(N, P);
    Zcenter = (Z.transpose().colwise() - B0).transpose();
    MatrixXd AlphaM(P, N); 
    MatrixXd AlphaMX(P, N);
    AlphaM = ((M * B).transpose() + Z.transpose()).array().exp();
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        AlphaM.col(i) = (B.transpose() * M.row(i).transpose() + Z.row(i).transpose()).array().exp();
        AlphaMX.col(i) = AlphaM.col(i) + X.row(i).transpose();
    }
    MatrixXd S(P, P);
    S = Zcenter.transpose() * Zcenter / N;

    double obj;
    computeAlphaMatrix(AlphaMX, AlphaM, obj);
    obj /= N;

    obj = obj - computeLnDet(Theta) / 2 + (S * Theta).diagonal().sum() / 2 + lambda1 * Theta.array().abs().sum() / 2 + lambda2 * B.array().abs().sum();
    return obj;
}

// Compute the edges of Theta
int computeEdgesTheta(MatrixXd& Theta) {
    int edges = 0;
    int p = Theta.rows();
    for(int i = 0; i < p; ++i) {
        for(int j = 0; j < i; ++j) {
            if(Theta(i,j) != 0)
                edges += 1;
        }
    }

    return edges;
}

// Compute the edges of B
int computeEdgesB(MatrixXd& B) {
    int edges = 0;
    int q = B.rows();
    int p = B.cols();
    for(int i = 0; i < q; ++i) {
        for(int j = 0; j < p; ++j) {
            if(B(i,j) != 0)
                edges += 1;
        }
    }

    return edges;
}

// Compute the EBIC
double mLDM::computeEBIC(double& objNew, MatrixXd& B, MatrixXd& Theta, double& lambda1, double& lambda2) {
    double g = 0.5;
    int E1 = computeEdgesTheta(Theta);
    int E2 = computeEdgesB(B);
    double EBIC = 2 * N * (objNew - lambda2 * B.array().abs().sum() - lambda1 * Theta.array().abs().sum() / 2)
    + (E1 + E2) * log(N) + 4 * g * E1 * log(P) + 2 * g * E2 * log(P*Q);

    return EBIC;
}

// Compute vary of edges of Theta
int computeEdgesVaryTheta(MatrixXd& Theta, MatrixXd& ThetaOld) {
    int edges = 0;
    int p = Theta.rows();
    bool per1, per2;

    for(int i = 0; i < p; ++i) {
        for(int j = 0; j < i; ++j) {
            (Theta(i,j) != 0) ? per1 = true : per1 = false;
            (ThetaOld(i,j) != 0) ? per2 = true : per2 = false;
            if(! (per1&& per2)) {
                edges += 1;
            }
        }
    }

    return edges;
}

// Compute vary of edges of B
int computeEdgesVaryB(MatrixXd& B, MatrixXd& BOld) {
    int edges = 0;
    int q = B.rows();
    int p = B.cols();
    bool per1, per2;

    for(int i = 0; i < q; ++i) {
        for(int j = 0; j < p; ++j) {
            (B(i,j) != 0) ? per1 = true : per1 = false;
            (BOld(i,j) != 0) ? per2 = true : per2 = false;
            if(! (per1 && per2)) {
                edges += 1;
            }
        }
    }

    return edges;
}

// compute the derivative of matrix B
MatrixXd derObjB(MatrixXd& B, MatrixXd& ZB0, MatrixXd& X, MatrixXd& M, int& n) {
    MatrixXd AlphaM = ((M * B).transpose() + ZB0).array().exp();
    MatrixXd AlphaMX = AlphaM + X.transpose();
    MatrixXd der(X.cols(), X.rows());
    computeAlphaDerMatrix(AlphaMX, AlphaM, der);
    MatrixXd derB = (der.array() * AlphaM.array()).matrix() * M;

    return derB.transpose() / n;
}

// lognormal-dirichlet-multinomial model
void mLDM::LognormalDirichletMultinomial(MatrixXd& X, MatrixXd& M, MatrixXd& B, VectorXd& B0, MatrixXd& Theta, MatrixXd& Z,  bool verbose, bool loopFlag, mLDMResult& solution) {
    // Record values of the objective function between two iterations
    double objOldB = 0.0;
    double objNewB = 0.0;
    double objOldTheta = 0.0;
    double objNewTheta = 0.0;
    double delta = 0.0;
    double deltaB = 0.0;
    double deltaTheta = 0.0;
    double EBIC = 0.0;
    double objNew = 0.0;
    double objOld = 0.0;

    // The round of iterations
    int round = 0;
    int roundB = 0;
    int roundTheta = 0;
    // check if parameters are all finite
    bool isFinite = true;

    // Record old values for B and Theta
    MatrixXd BOld = B;
    MatrixXd ThetaOld = Theta;
    VectorXd B0Old = B0;
    MatrixXd ZOld = Z;
    // record element delta max 
    double B_delta_max = 0.0;
    double Theta_delta_max = 0.0;

    MatrixXd derZ(N, P);
    std::cout << "Origin parameters: " << N << " " << P << " " << Q << std::endl;
    while(true) {
        round += 1;
        //if (verbose) {
        //    std::cout << "round " << round << std::endl;
        //}
        if (N <= P) {
            std::cout << "The number of samples <= the number of OTUs, stop estimating !!!" << std::endl;
            isFinite = false;
            break;
        }

        double z_start, z_end;
        z_start = omp_get_wtime();
        double init_sum = 0;
        double solve_sum = 0;
        double post_sum = 0;
        // record old Z
        ZOld = Z;
        // Estimate Z parallel
        #pragma omp parallel for
        for(int i = 0; i < N; ++i) {
            // if z is not finite, use LBFGSZ_CERES to re-estimate
            ZiData zi_data;
            VectorXd zi = Z.row(i);
            VectorXd xi = X.row(i);
            zi_data.Xi = xi;
            VectorXd mui = B.transpose() * M.row(i).transpose();
            zi_data.Mui = mui; 
            zi_data.N = N;
            zi_data.Theta = Theta;
            zi_data.B0 = B0;
            lbfgsfloatval_t fx;
            lbfgsfloatval_t* zi_x = lbfgs_malloc(P);
            lbfgs_parameter_t lbfgs_param;
            lbfgs_parameter_init(&lbfgs_param);
            lbfgs_param.max_iterations = max_iteration_Z;
            lbfgs_param.epsilon = eplison_Z;
            lbfgs_param.max_linesearch = max_linesearch_Z;
            lbfgs_param.m = approx_num_Z;
            for (int iv = 0; iv < P; ++iv) zi_x[iv] = zi[iv];
            int ret = lbfgs(P, zi_x, &fx, lbfgs_evaluate, NULL, &zi_data, &lbfgs_param);
            for (int iv = 0; iv < P; ++iv) zi[iv] = zi_x[iv];
            lbfgs_free(zi_x);
            // if fail 
            if (! zi.allFinite()) {
                std::cout << "zi is not all finite !" << std::endl;
                zi = Z.row(i);
                lbfgsfloatval_t fx;
                lbfgsfloatval_t* zi_x = lbfgs_malloc(P);
                lbfgs_parameter_t lbfgs_param2;
                lbfgs_parameter_init(&lbfgs_param2);
                lbfgs_param2.max_iterations = max_iteration_Z;
                lbfgs_param2.epsilon = eplison_Z;
                lbfgs_param2.max_linesearch = max_linesearch_Z;
                lbfgs_param2.m = approx_num_Z;
                lbfgs_param2.linesearch = LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE;
                for (int iv = 0; iv < P; ++iv) zi_x[iv] = zi[iv];
                int ret = lbfgs(P, zi_x, &fx, lbfgs_evaluate, NULL, &zi_data, &lbfgs_param2);
                for (int iv = 0; iv < P; ++iv) zi[iv] = zi_x[iv];
                lbfgs_free(zi_x);
                // update the Z.row(i) and its derivative
                if (! zi.allFinite()) {
                    isFinite = false;
                }
            }

            Z.row(i) = zi.transpose();
            VectorXd derZi(P);
            double* zi_v = new double[P];
            double zi_grad[P];
            for (int iv = 0; iv < P; ++iv) zi_v[iv] = zi[iv];
            double rv = lbfgs_evaluate(&zi_data, zi_v, zi_grad, P, 0);
            derZi = Eigen::Map<VectorXd>(zi_grad, P);
            derZ.row(i) = derZi;
        }
        // std::cout << "der Z norm " << derZ.norm() << std::endl;
        z_end = omp_get_wtime();
        // std::cout << "z estimate time: " << (z_end - z_start) * 1.0 << std::endl;
        // check if elements of Z are all finite
        if (! isFinite) {
            std::cout << "Inf or Nan found in latent variable Z" << std::endl;
            break;
        }

        // Estimate B0 
        // record old B0 
        B0Old = B0;
        B0 = Z.colwise().mean();

        // Estimate B
        if(loopFlag) {
            roundB += 1;
            objOldB = objNewB;
            // std::cout << "round B: " << roundB << std::endl;
            double b_start, b_end;
            b_start = omp_get_wtime();
            double b_solve_sum = 0;
            
            // Optimize the B[i,] respectively
            MatrixXd ZB0 = Z.transpose();
            // record old B
            BOld = B;
            // std::cout << "start B" << std::endl;
            for(int i = 0 ; i < Q; ++i) {
                // init problem and parameters
                PQNPoint pqn_point(P, Q, N, X, M);
                VectorXd bv = B.row(i).transpose();
                MatrixXd Bp = B;
                Bp.row(i).setZero();
                MatrixXd BZ = (M * Bp).transpose() + ZB0;
                // set the parameter BZ
                pqn_point.setBZ(BZ);
                pqn_point.setIndex(i);
                // proximal qusi-newton to estimate B.row(i)
                PQNResult pqn_result;
                double solve_s, solve_e;
                solve_s = omp_get_wtime();
                proximalQusiNewton(bv, solution.lambda2, approx_num_B, max_linesearch_B, max_iteration_B,
                                                   threshold_B, delta1_threshold_B, delta2_threshold_B, sy_threshold_B,
                                                   max_iteration_B_coor, threshold_B_coor, verbose, pqn_result, pqn_point);
                solve_e = omp_get_wtime();
                b_solve_sum += (solve_e - solve_s);
                // update B.row(i)
                if (! pqn_result.w.allFinite()) {
                    isFinite = false;
                }
                B.row(i) = pqn_result.w.transpose();
                /*std::cout << "i: " << i << std::endl;
                std::cout << "Bi: " << pqn_result.w << std::endl;
                std::cout << "Bi obj: " << pqn_result.obj << std::endl;*/
            }
            b_end = omp_get_wtime();
           // std::cout << "B estimate time: " << (b_end - b_start) * 1.0 << std::endl;
            // std::cout << "B solve time: " << b_solve_sum << std::endl; 
            // check if elements of B are all finite
            if (! isFinite) {
                isFinite = false;
                std::cout << "Inf or Nan found in OTU-Meta association matrix B" << std::endl;
                break;
            }

            objNewB = computeObjf(X, M, B, Theta, Z, solution.lambda1, solution.lambda2);
            // judge if it is nan or inf 
            if (! std::isfinite(objNewB) || ! std::isnormal(objNewB)) {
                isFinite = false;
                std::cout << "value of obj function is nan or inf when estimate B !!!" << std::endl;
                break;
            }
            // std::cout << "B obj: " << objNewB << std::endl;
            deltaB = fabs(objNewB - objOldB);
            // compute element delta 
            B_delta_max = (B - BOld).array().abs().maxCoeff();
            // std::cout << "B delta max: " << B_delta_max << std::endl;

            // B is finished
            if(deltaB < threshold || round >= max_iterations || B_delta_max < threshold) {
                loopFlag = false;
                objNew = objNewB;
                objOld = objNewTheta;

                if(verbose) {
                    MatrixXd ZB0 = Z.transpose();
                    MatrixXd derB = derObjB(B, ZB0, X, M, N);

                    std::cout << "der Z norm " << derZ.norm() << std::endl;

                    std::cout << "Max B " << B.maxCoeff() << "min B " << B.minCoeff() << std::endl;
                    std::cout << "Round B is:" << roundB << std::endl;
                    std::cout << "Norm2 B is:" << sqrt(B.squaredNorm()) << std::endl;
                    std::cout << "Norm2 derB is:" << sqrt(derB.squaredNorm()) << std::endl;
                }
                roundB = 0;
            } else {
                // std::cout << "delta B: " << deltaB << std::endl;
                continue;
            }
        } else {// Estimate Theta
            roundTheta += 1;
            objOldTheta = objNewTheta;

            // Estimate Theta via graphical lasso
            MatrixXd Zcenter(N, P);
            Zcenter = (Z.transpose().colwise() - B0).transpose();
            MatrixXd S(P, P);
            S = Zcenter.transpose() * Zcenter / N;
            // QUIC to estimate Theta
            double theta_start, theta_end;
            theta_start = omp_get_wtime();
            // record Theta old
            ThetaOld = Theta;
            // std::cout << "start Theta" << std::endl;
            QUICC(S, solution.lambda1, max_iteration_Theta, threshold_Theta, msg_Theta, Theta);
            theta_end = omp_get_wtime();

            // std::cout << "Theta estimate time: " << (theta_end - theta_start) * 1.0 << std::endl;

            // update the Theta
            if (! Theta.allFinite()) {
                isFinite = false;
                std::cout << "Inf or Nan found in OTU-OTU association matrix Theta" << std::endl;
                break;
            }

            objNewTheta = computeObjf(X, M, B, Theta, Z, solution.lambda1, solution.lambda2);
            if (! std::isfinite(objNewTheta) || ! std::isnormal(objNewTheta)) {
                isFinite = false;
                std::cout << "value of obj function is nan or inf when estimate Theta !!!" << std::endl;
                break;
            }
            // std::cout << "Theta obj: " << objNewTheta << std::endl;
            deltaTheta = fabs(objNewTheta - objOldTheta);
            // compute element delta max 
            Theta_delta_max = (Theta - ThetaOld).array().abs().maxCoeff();
            // std::cout << "Theta delta max: " << Theta_delta_max << std::endl;

            if(deltaTheta < threshold || round >= max_iterations || Theta_delta_max < threshold) {
                loopFlag = true;
                objNew = objNewTheta;
                objOld = objNewB;
                if(verbose) {
                    std::cout << "der Z norm " << derZ.norm() << std::endl;

                    std::cout << "Round Theta is:" << roundTheta << std::endl;
                }
                roundTheta = 0;
            } else {
                // std::cout << "delta theta: " << deltaTheta << std::endl;
                continue;
            }
        }

        if(verbose) {
            int edgesTheta = computeEdgesTheta(Theta);
            int edgesB = computeEdgesB(B);
            std::cout << "edges of Theta:" << edgesTheta << std::endl;
            std::cout << "edges of B:" << edgesB << std::endl;
            int varyTheta = computeEdgesVaryTheta(Theta, ThetaOld);
            int varyB = computeEdgesVaryB(B, BOld);
            std::cout << "vary of edges Theta:" << varyTheta << std::endl;
            std::cout << "vary of edges B:" << varyB << std::endl;
            std::cout << "objNew " << objNew << std::endl;
            std::cout << "objOld " << objOld << std::endl;
        }

        // Judge termination
        delta = fabs(objNew - objOld);
        if(verbose) {
            std::cout << "delta is : " << delta  << std::endl;
            std::cout << "lambda1 :" << solution.lambda1 << std::endl;
            std::cout << "lambda2 :" << solution.lambda2 << std::endl;
        }

        if(delta < threshold || round >= max_iterations) {
            if(verbose) {
                std::cout << "Rounds of iterations : " << round << std::endl;
                std::cout << "Iteration Success ~" << std::endl;
            }
            break;
        }
    }

    // Compute the EBIC, return results
    // save the result
    if (! isFinite) {
        solution.exists = true;
        objNew = computeObjf(X, M, BOld, ThetaOld, ZOld, solution.lambda1, solution.lambda2);
        EBIC = computeEBIC(objNew, BOld, ThetaOld, solution.lambda1, solution.lambda2);
        solution.exists = true;
        solution.EBIC = EBIC;
        solution.B = BOld;
        solution.Theta = ThetaOld;
        solution.B0 = B0Old;
        solution.obj = objNew;

    } else {
        EBIC = computeEBIC(objNew, B, Theta, solution.lambda1, solution.lambda2);
        solution.exists = true;
        solution.EBIC = EBIC;
        solution.B = B;
        solution.Theta = Theta;
        solution.B0 = B0;
        solution.obj = objNew;
    }

    return;
}


// estimate the association netoworks via mLDM model
// otus:  N * P matrix otu table
// efs:   N * Q matrix meta table
void mLDM::work() {
    // initialize the parameters
    // get the relative abundance matrix
    MatrixXd X_relative(N, P);
    MatrixXd X_d = X.cast <double> ();
    toRelative(X_d, X_relative);
    // normalize the meta data
    MatrixXd M_scale(N, Q);
    scale(M, M_scale, true);

    // correlation among otus
    MatrixXd cor_otus(P, P);
    corr(X_d, "spearman", cor_otus);
    cor_otus = cor_otus.array().abs();
    // correlation between otus and meta data
    MatrixXd cor_otus_meta(P, Q);
    corr(X_d, M, "spearman", cor_otus_meta);
    cor_otus_meta = cor_otus_meta.array().abs();

    if(!cor_otus.allFinite() || ! cor_otus_meta.allFinite()) {
        std::cout << "Initialization Error! NaN found in the spearman correlation matrix among X, or between X and M" << std::endl;
        exit(0);
    }

    // For Z
    MatrixXd Z_init = (X_d.array() + 1).array().log() + 1.0;
    if(!Z_init.allFinite()) {
        std::cout << "Initialization Error! NaN found in the Z" << std::endl;
        exit(0);
    }
    // For B
    MatrixXd B_init = cor_otus_meta.transpose();
    if(!B_init.allFinite()) {
        std::cout << "Initialization Error! NaN found in the B" << std::endl;
        exit(0);
    }
    // For Theta
    MatrixXd Theta_cov = cor_otus;
    double diagValue = 1;
    while(Theta_cov.determinant() < P) {
        Theta_cov.diagonal().array() += diagValue;
    }
    MatrixXd Theta_init = Theta_cov.inverse();
    if(!Theta_init.allFinite()) {
        std::cout << "Initialization Error! NaN found in the Theta" << std::endl;
        exit(0);
    }
    // For B0
    VectorXd B0_init = Z_init.colwise().mean();
    if(!B0_init.allFinite()) {
        std::cout << "Initialization Error! NaN found in the B0" << std::endl;
        exit(0);
    }
    // set the lambda1 and lambda2 list
    std::vector<double> cor1_list;
    getSymmStrictLowerValues(cor_otus, cor1_list);
    std::vector<double> cor2_list;
    getMatrixValuesList(cor_otus_meta, cor2_list);

    std::vector<double> ratios_list;
    ratios_list.push_back(ratio1);
    ratios_list.push_back(ratio2);
    // select the quantiles of correlations as the lambda values
    std::vector<double> lambda1_bound;
    std::vector<double> lambda2_bound;
    quantile(cor1_list, ratios_list, lambda1_bound);
    quantile(cor2_list, ratios_list, lambda2_bound);

    // get the lambda1 list and lambda2 list
    std::vector<double> lambda1_list;
    std::vector<double> lambda2_list;
    if (isZero(std::fabs(lambda1_bound[0] - lambda1_bound[1]))) {
        lambda1_list.push_back(lambda1_bound[0]);
    } else {
        seq(lambda1_bound[0], lambda1_bound[1], model_selection_num, lambda1_list);
    }

    if (isZero(std::fabs(lambda2_bound[0] - lambda2_bound[1]))) {
        lambda2_list.push_back(lambda2_bound[0]);
    } else {
        seq(lambda2_bound[0], lambda2_bound[1], model_selection_num, lambda2_list);
    }

    int length1 = lambda1_list.size();
    int length2 = lambda2_list.size();

    // Record the minimum of EBIC
    double EBIC_min = 1e+20;
    // LoopFlag to control the order of optimaztion for B or Theta
    // flase optimize the Theta first
    // true optimize the B first
    bool loopFlag = false;
    int count = 0;
    // If the optimal lambda1 is found
    bool exist_best = false;
    // is the first solution
    bool is_first = true;
    // If is the optimal lambda1
    int best_i = -1;
    // save the index of minimun EBIC
    int min_index = -1;
    std::vector<bool> exist_per;
    for (int i = 0; i < length1; ++i) {
        exist_per.push_back(false);
    }
    // save the results of every mLDM model with various combinations of lambda1 and lambda2
    std::vector<mLDMResult> LDM_result_all;

    // select the best result from various combinations of lambda1 and lambda2
    for (int j = 0; j < length2; ++j) {
        for (int i = 0; i < length1; ++i) {
            j > 0 ? loopFlag = true : loopFlag = false;

            if (! exist_best || exist_per[i]) {
                double lambda1 = lambda1_list[i];
                double lambda2 = lambda2_list[j];

                if (verbose) {
                    std::cout << "lambda1 " << lambda1 << " lambda2 " << lambda2 << std::endl;
                }

                mLDMResult per_solution;
                per_solution.lambda1 = lambda1;
                per_solution.lambda2 = lambda2;
                per_solution.EBIC = 1e+30;

                LognormalDirichletMultinomial(X_d, M_scale, B_init, B0_init, Theta_init, Z_init, verbose, loopFlag, per_solution);

                if(per_solution.exists && per_solution.EBIC > 0 && per_solution.EBIC < EBIC_min) {
                    EBIC_min = per_solution.EBIC;
                    best_i = i;
                    // save the position of best result
                    min_index = LDM_result_all.size();

                    if (is_first) {
                        B_init = per_solution.B;
                        B0_init = per_solution.B0;
                        Theta_init = per_solution.Theta;

                        is_first = false;
                    }
                }
                LDM_result_all.push_back(per_solution);
            }

        }

        if(best_i != -1) {
            exist_per[best_i] = true;
            exist_best = true;
        }
    }

    // extract the best result
    if (min_index != -1) {
        best_result = LDM_result_all[min_index];
        best_result.exists = true;
    } else {
        best_result.lambda1 = 0.0;
        best_result.lambda2 = 0.0;
        best_result.EBIC = -1;
        best_result.exists = false;
    }

    return;
}
