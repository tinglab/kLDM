#ifndef PQN_H
#define PQN_H

#include "utils.h"

class Point {
public:
    virtual bool Evaluate(const VectorXd&x, double& obj, VectorXd& grad) const = 0;
    // virtual double value(VectorXd& x) = 0;
    // virtual void gradient(VectorXd& x, VectorXd& grad) = 0;
};

struct LineSearchResult {
    VectorXd w;
    VectorXd grad;
    double obj;
    bool exist;
};

struct PQNResult {
    VectorXd w;
    double obj;
    VectorXd gradient;
};

void proximalQusiNewton(VectorXd& w, double lambda, int approx_num, int max_linesearch,
                        int max_iteration, double threshold, double delta1_threshold,
                        double delta2_threshold, double sy_threshold, int max_iteration_coor,
                        double threshold_coor, bool verbose, PQNResult& solution, Point& pqn_point);

#endif
