#include "proximalQusiNewton.h"

MatrixXd cBind(MatrixXd a, MatrixXd b) {
    int row = a.rows();
    int col = a.cols() + b.cols();
    MatrixXd c(row, col);
    
    c.block(0,0,row, a.cols()) = a;
    c.block(0,a.cols(),row, b.cols()) = b;
    
    return c;
}

MatrixXd rBind(MatrixXd a, MatrixXd b) {
    int col = a.cols();
    int row = a.rows() + b.rows();
    MatrixXd c(row, col);
    
    c.block(0,0,a.rows(), col) = a;
    c.block(a.rows(), 0,b.rows(), col) = b;
    
    return c;
}

void computeBtQ(int& d, double& gammat, MatrixXd& St, MatrixXd& Yt, int& approx_num, int& approxCount, int& approxIndex, MatrixXd& Q, MatrixXd& Qh) {
    if(approxCount != 0) {
        MatrixXd Stp(d, approxCount);
        MatrixXd Ytp(d, approxCount);
        
        int basis = 0;
        if(approxCount == approx_num)
            basis = approxIndex + 1;
        for(int i = 0; i < approxCount; ++i) {
            Stp.col(i) = St.col((i + basis) % approx_num);
            Ytp.col(i) = Yt.col((i + basis) % approx_num);
        }
        
        // cout << "Stp " << Stp << endl;
        // cout << "Ytp " << Ytp << endl;
        
        Q = cBind(gammat*Stp, Ytp);
        // cout << "Qin" << Q << endl;
        // cout << "gammat" << gammat << endl;
        MatrixXd SY = Stp.transpose() * Ytp;
        MatrixXd Dt = SY.diagonal().asDiagonal();
        
        MatrixXd Lt = SY.triangularView<Eigen::StrictlyLower>();
        // cout << "SY" << SY << endl;
        // cout << "Dt" << Dt << endl;
        MatrixXd Rt = rBind(cBind(gammat*Stp.transpose() * Stp, Lt), cBind(Lt.transpose(), - Dt));
        Qh = Rt.inverse() * Q.transpose();
        // cout << "Qhin" << Qh << endl;
    }
    
    return;
}

void getActiveSet(VectorXd& w, VectorXd& gt, double& lambda, int& d, std::vector<int>& index) {
    double threshold = 1e-6;
    double a, b, subgt;
    bool addFlag = false;
    
    for(int i = 0; i < d; ++i) {
        a = w[i];
        b = gt[i];
        addFlag = true;
        if(fabs(a) < threshold) {
            subgt = fmax(0.0, fmax(b - lambda, -(b + lambda)));
            
            if(fabs(subgt) < threshold) {
                addFlag = false;
            }
        }
        
        if(addFlag) {
            index.push_back(i);
        }
    }
    
    return;
}

double softthreshold(double a, double b) {
    return a > 0 ? fmax(a - b, 0.0) : -1 * fmax(-a - b, 0.0);
}

// Compute the direction via coordinate descent
VectorXd coordinateDescent(int& d, VectorXd& w, double& gammat, MatrixXd& Q, MatrixXd& Qh, VectorXd& gt, double& lambda,
                           int& max_iteration, double& threshold) {
    VectorXd wt = w;
    VectorXd wtOld;
    std::vector<int> activeSet;
    getActiveSet(w, gt, lambda, d, activeSet);
    int activeSize = activeSet.size();
    
    int round = 0;
    double delta = 1;
    int i;
    double a, b;
    MatrixXd Bt = gammat * MatrixXd::Identity(d,d) - Q * Qh;
    
    while(delta > threshold && round < max_iteration) {
        round += 1;
        wtOld = wt;
        
        for(int index = 0 ; index < activeSize; ++index) {
            i = activeSet[index];
            a = Bt(i,i);
            b = gt[i] + (Bt.row(i) * (wt - w)) - a * wt[i];
            
            wt[i] = softthreshold(-b, lambda) / a;
        }
        
        delta = (wt - wtOld).array().abs().sum() / d;
        
    }
    
    VectorXd direction = wt - w;
    return direction.normalized();
}

void filterDir(VectorXd& w) {
    int d = w.rows();
    for(int i = 0 ; i < d; ++i) {
        if(fabs(w[i]) < 1e-10)
            w[i] = 0;
    }
    return;
}

// Find new w via linesearch based on strong wolfe condition
void linesearch(VectorXd& w, VectorXd& direction, double& lambda, int& max_linesearch, double& f0, VectorXd& g0, double& delta1, double& delta2, Point& pqn_point, LineSearchResult& linesearch_result) {
    double beta = 0.5;
    double alpha = 2;
    VectorXd wt = w;
    int k = 0;
    
    double f1 = f0;
    double dg1 = 0;
    VectorXd g1 = g0;
    
    if(!direction.allFinite() || !g0.allFinite()) {
        std::cout << "linesearch direction failed!" << std::endl;
        linesearch_result.exist = false;
        linesearch_result.w = w;
        linesearch_result.obj = f0;
        linesearch_result.grad = g0;
    }
    
    double d1 = g0.transpose() * direction;
    double dg0 = d1 + lambda * w.array().abs().sum();
    d1 += lambda * ((w + direction).array().abs().sum() - w.array().abs().sum());
    d1 *= delta1;
    
    double part;
    
    //double line_s, line_e;
    //double line_sum;
    VectorXd f1_grad(w.size());
    double f1_func;
    bool exist = false;
    while(k <= max_linesearch) {
        alpha *= beta;
        wt = w + alpha * direction;
        // line_s = omp_get_wtime();
        pqn_point.Evaluate(wt, f1_func, f1_grad);
        // line_e = omp_get_wtime();
        // line_sum += (line_e - line_s);
        f1 = f1_func + lambda * wt.array().abs().sum();
        part = alpha * d1;
        k += 1;
        
        if(f1 <= f0 + part) {
            dg1 = f1_grad.transpose() * direction + alpha * lambda * (w + direction).array().abs().sum();
            if(fabs(dg1 / dg0) <= delta2) {
               /* std::cout << "wt: " << wt << std::endl;
                std::cout << "f1: " << f1 << std::endl;
                std::cout << "f1 grad: " << f1_grad << std::endl;*/
                g1 = f1_grad;
                exist = true;
                break;
            }
        }
    }

    // std::cout << "line k times : " << k << std::endl;
    /*std::cout << "line func evaluate time: " << line_sum << std::endl; 
    std::cout << "func evaluate time: " << line_sum / k << std::endl;
    */
    filterDir(wt);
    // 如果达到最大搜索次数还没有满足条件,也是返回true?需要修改?
    // save the final linesearch result
    linesearch_result.exist = exist;
    linesearch_result.w = wt;
    linesearch_result.obj = f1;
    linesearch_result.grad = g1;
    
    return;
}

// priximal qusi-newton method
void proximalQusiNewton(VectorXd& w, double lambda, int approx_num, int max_linesearch,
                        int max_iteration, double threshold, double delta1_threshold,
                        double delta2_threshold, double sy_threshold, int max_iteration_coor,
                        double threshold_coor, bool verbose, PQNResult& solution, Point& pqn_point) {
    int approxCount = 0;
    int approxIndex = -1;
    double gammat = 1;
    double gammat0 = 1;
    
    int d = w.rows();
    MatrixXd St(d, approx_num);
    St.setZero();
    MatrixXd Yt(d, approx_num);
    Yt.setZero();
    MatrixXd Q(d, 2);
    Q.setZero();
    MatrixXd Q0(d, 2);
    Q0.setZero();
    MatrixXd Qh(2, d);
    Qh.setZero();
    MatrixXd Qh0(2, d);
    Qh0.setZero();
    VectorXd direction(d);
    VectorXd wt = w;
    
    int round = 0;
    // When SYdot < sy_threshold, change to steepest descent
    bool steepestActive = false;
    int steepestActiveHeight = 1;
    int steepestCount = 0;
    
    // Record obj values, the gradient for B
    double objNew = 0;
    double obj_func;
    VectorXd gt0(d);
    pqn_point.Evaluate(w, obj_func, gt0);
    double objOld = obj_func + lambda * w.array().abs().sum();
    double delta = 0;
    // cout << "gt0" << gt0 << endl;
    VectorXd gt1(d);
    VectorXd StPer(d);
    VectorXd YtPer(d);
    double SYDot = 0;
    bool exist;
    LineSearchResult linesearch_result;
    
    double coor_s, coor_e;
    double dir_s, dir_e;
    double post_s, post_e;
    double coor_sum, dir_sum, post_sum;
    while(true) {
        round += 1;
        // coor_s = omp_get_wtime();
        // Obtain direction
        if(!steepestActive) {
            computeBtQ(d, gammat, St, Yt, approx_num, approxCount, approxIndex, Q, Qh);
            direction = coordinateDescent(d, w, gammat, Q, Qh, gt0, lambda, max_iteration_coor, threshold_coor);
        } else { // when the hessian nears singular, change to steepest descent
            direction = coordinateDescent(d, w, gammat0, Q0, Qh0, gt0, lambda, max_iteration_coor, threshold_coor);
            steepestCount += 1;
            if(steepestCount >= steepestActiveHeight) {
                steepestActive = false;
            }
        }
        // coor_e = omp_get_wtime();
        // coor_sum += (coor_e - coor_s);
        // dir_s = omp_get_wtime();
        // line search find new w
        linesearch(w, direction, lambda, max_linesearch, objOld, gt0, delta1_threshold,
                          delta2_threshold, pqn_point, linesearch_result);
        // dir_e = omp_get_wtime();
        // dir_sum += (dir_e - dir_s);
        // post_s = omp_get_wtime();
        exist = linesearch_result.exist;
        if(!exist) {
            break;
        }
        
        wt = linesearch_result.w;
        gt1 = linesearch_result.grad;
        objNew = linesearch_result.obj;
        // cout << "wt " << wt << endl;
        
        delta = objNew - objOld;
        StPer = wt - w;
        YtPer = gt1 - gt0;
        SYDot = StPer.transpose() * YtPer;
        if(SYDot > sy_threshold) {
            approxCount > approx_num ? approxCount = approx_num: approxCount += 1;
            
            approxIndex = (approxIndex + 1) % approx_num;
            St.col(approxIndex) = StPer;
            Yt.col(approxIndex) = YtPer;
            gammat = SYDot / StPer.squaredNorm();
        } else {
            steepestActive = true;
            steepestCount = 0;
            if(verbose) {
                //cout << "steepest descent active!" << endl;
            }
        }
        
        objOld = objNew;
        gt0 = gt1;
        w = wt;
        // post_e = omp_get_wtime();
        // post_sum += (post_e - post_s);
        if(fabs(delta) < threshold || round > max_iteration) {
            break;
        }
    }

   // std::cout << "round pqn: " << round << std::endl;
   /* std::cout << "coor time : " << coor_sum << std::endl;
    std::cout << "dir time : " << dir_sum << std::endl;
    std::cout << "post time : " << post_sum << std::endl; */
    // save the final results
    solution.w = w;
    solution.obj = objOld;
    solution.gradient = gt0;

    return;
}
