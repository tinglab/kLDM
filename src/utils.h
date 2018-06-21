#ifndef UTILS_H
#define UTILS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/SpecialFunctions>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cstring>
#include <math.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include <map>
#include <omp.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using Eigen::SelfAdjointEigenSolver;
using Eigen::ArrayXd;

#define PI 3.14159265358979323846

struct ValueOrder {
    double value;
    int index;
};

// compute the log det of precision matrix
double computeLnDet(MatrixXd& cov);

// some functions for matrix calculations
void computeVecMatMult(VectorXd& a, MatrixXd& b, MatrixXd& c);
void computeVecMatMultTranspose(VectorXd& a, MatrixXd& b, MatrixXd& c);

// compute the variance of variable
double computeVar(VectorXd v);

// test if the variable is all zero
bool isZero(double a);

// count the zero num of variable
int zeroNum(VectorXd v);

// compute the covariance matrix
void computeCov(MatrixXd& a, MatrixXd& cov_a);
void computeCov(MatrixXd& a, MatrixXd& b, MatrixXd& cov_ab);

// convert every column into vaiable with mean = 0 and variance = 1
void scale(MatrixXd& a, MatrixXd& a_scale, bool skip_zero = false);
// convert count matrix into relative abundance matrix
// the sum of every row is 1.0
void toRelative(MatrixXd& a, MatrixXd& a_relative);

// compute the correlation among columns of the matrix
// method can be "pearson" or "spearman"
void corr(MatrixXd& a, std::string method, MatrixXd& result);

// compute the correlation between matrix a and b
void corr(MatrixXd& a, MatrixXd& b, std::string method, MatrixXd& result);

// compute the quantile of the variable
void quantile(std::vector<double> temp_a, std::vector<double>& ratios, std::vector<double>& a_ratios);

// convert the element of matrix a into the index of element
void convertMatrixValue2Index(MatrixXd& a, MatrixXd& a_index);

// construct n values list with the same interval from a to b
void seq(double a, double b, int num, std::vector<double>& sequences);

// compute the union set of two integer sets
void getUnionList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& merged_list);
// compute the difference set of two integer sets
// exist in list_a but not in list_b
void getDiffList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& diff_list);
// compute the intersection set of two integer sets
// exist both in list_a and in list_b
void getIntersectList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& intersect_list);

// remove the i-th row of matrix X, then update the X
void removeMatrixRow(MatrixXd& X, int i);

#endif
