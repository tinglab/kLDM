#include "utils.h"

// compute the log det of precision matrix
double computeLnDet(MatrixXd& cov) {
    SelfAdjointEigenSolver<MatrixXd> es(cov);
    // std::cout << "eigen values " << es.eigenvalues() << std::endl;
    return es.eigenvalues().array().log().sum();
}

// point multiplication for vector A * matrix B
// the i-th column of C equals the vector A mutiple the i-th column of matrix B
void computeVecMatMult(VectorXd& a, MatrixXd& b, MatrixXd& c) {
    for (int i = 0; i < b.cols(); ++ i) {
        c.col(i) = a.array() * b.col(i).array();
    }
}

// point multiplication for vector A * matrix B
// the i-th row of C equals the transpose of the vector A mutiple the i-th row of matrix B
void computeVecMatMultTranspose(VectorXd& a, MatrixXd& b, MatrixXd& c) {
    for (int i = 0; i < b.rows(); ++ i) {
        c.row(i) = a.transpose().array() * b.row(i).array();
    }
}

// compute the variance of variable
double computeVar(VectorXd v) {
    if (v.size() <= 1) {
        return 0;
    }
    double v_mean = v.mean();
    double var = (v.array() - v_mean).matrix().squaredNorm() / (v.size() - 1);

    return var;
}

// test if the variable is zero
bool isZero(double a) {
    if (fabs(a) < 1e-10) {
        return true;
    } else {
        return false;
    }
}

// test if the variable is all zero
int zeroNum(VectorXd v) {
    return (v.array().abs() < 1e-10).count();
}

// compute the covariance matrix
// x -- N * P
// cov_x -- P * P
void computeCov(MatrixXd& x, MatrixXd& cov_x) {
    if (x.rows() == 1) {
        cov_x = MatrixXd::Identity(x.cols(), x.cols());
    } else {
        MatrixXd x_center = x.rowwise() - x.colwise().mean();
        cov_x = x_center.transpose() * x_center / (x.rows() - 1);
    }
}

// compute the covariance between two matrix
// a -- N * P
// b -- N * Q
// cov_ab -- P * Q
void computeCov(MatrixXd& a, MatrixXd& b, MatrixXd& cov_ab) {
    int n = a.rows();
    int p = a.cols();
    int q = b.cols();
    if (n != b.rows()) {
        std::cout << "the sample size of matrix a and b is not equal ! (Function: corr(Matrix a, Matrix b, ...))" << std::endl;

        exit(0);
    }

    MatrixXd a_center = a.rowwise() - a.colwise().mean();
    MatrixXd b_center = b.rowwise() - b.colwise().mean();
    cov_ab = a_center.transpose() * b_center / (n - 1);
}

// convert count into relative ratio
void toRelative(MatrixXd& a, MatrixXd& a_relative) {
    int n = a.rows();

    VectorXd xRowSums = a.rowwise().sum();
    for(int i = 0; i < n; ++i) {
        if (isZero(xRowSums[i])) {
            std::cout << "row with all zero elements found when convert the matrix into relative abundance; please check the matrix ! (Function: toRelative(Matrix a, Matrix a_relative))" << std::endl;

            exit(0);
        }
        a_relative.row(i) = a.row(i) * 1.0/ xRowSums[i];
    }

    return;
}

// Scale matrix with mean=0 and std=1
// skip_zero
// when find a column with same values, if divide the std or not
void scale(MatrixXd& a, MatrixXd& a_scale, bool skip_zero) {
    int q = a.cols();
    int n = a.rows();

    for(int i = 0; i < q; ++i) {
        VectorXd perCol = a.col(i);
        double perVar = computeVar(perCol);
        if (isZero(perVar)) {
            std::cout << "Find zero variance at column " << i << " of matrix a! (Function: scale(Matrix a, Matrix a_scale))" << std::endl;
            // std::cout << perCol << std::endl;
            if (! skip_zero) {
                exit(0);
            }
            a_scale.col(i) = VectorXd::Zero(n);
        } else {
            a_scale.col(i) = (perCol.array() - perCol.mean())/ sqrt(perVar);
        }
    }

    return;
}

// compute the quantile of the variable
void quantile(std::vector<double> temp_a, std::vector<double>& ratios, std::vector<double>& a_ratios) {
    int n = temp_a.size();
    // sort the variable a
    std::sort(temp_a.begin(), temp_a.end());
    for (int i = 0; i < ratios.size(); ++i) {
        double position = (n + 1) * ratios[i];
        int index = floor(position);

        if (position >= n) {
            a_ratios.push_back(temp_a[n - 1]);
        } else if (position <= 1) {
            a_ratios.push_back(temp_a[0]);
        } else {
            a_ratios.push_back(temp_a[index - 1] * (1 - (position - index)) + temp_a[index] * (position - index));
        }
    }
}

bool order_comp(const ValueOrder& a, const ValueOrder& b) {
    return a.value < b.value;
}

// convert the element of matrix a into the index of element
void convertMatrixValue2Index(MatrixXd& a, MatrixXd& a_index) {
    int p = a.cols();
    int n = a.rows();
    for (int i = 0; i < p; ++i) {
        std::vector<ValueOrder> per_column;
        for (int j = 0; j < n; ++j) {
            ValueOrder per_pair;
            per_pair.value = a(j, i);
            per_pair.index = j;
            per_column.push_back(per_pair);
        }
        // sort value and order at the same time
        sort(per_column.begin(), per_column.end(), order_comp);
        // convert value vector into index vector
        // use the fractional rank for the same value
        double batch_value = per_column[0].value;
        int batch_size = 1;
        int batch_index = 1;
        int k = 1;
        while (k < n) {
            if (per_column[k].value != batch_value) {

                for (int g = 0; g < batch_size; ++g) {
                    a_index(per_column[k - g - 1].index, i) = batch_index * 1.0 / batch_size;
                }

                batch_index = k + 1;
                batch_value = per_column[k].value;
                batch_size = 1;
            } else {
                batch_index += k + 1;
                batch_size += 1;
            }

            k += 1;
        }
        for (int g = 0; g < batch_size; ++g) {
            a_index(per_column[n - 1 - g].index, i) = batch_index * 1.0 / batch_size;
        }
    }
}

// method can be "pearson" or "spearman"
void corr(MatrixXd& a, std::string method, MatrixXd& result) {
    if (method == "spearman") {
        int p = a.cols();
        int n = a.rows();
        MatrixXd a_index(n, p);
        convertMatrixValue2Index(a, a_index);
        // compute the pearson correlation for a_index
        corr(a_index, "pearson", result);
    } else if (method == "pearson") {
        int p = a.cols();
        computeCov(a, result);
        VectorXd a_std(p);
        for (int i = 0; i < p; ++i) {
            a_std[i] = sqrt(result(i, i));
            if (isZero(a_std[i])) {
                std::cout << "Find zero variance ! (Function: corr(Matrix a, ...))" << std::endl;

                exit(0);
            }
        }
        for (int i = 0; i < p; ++ i) {
            for (int j = 0; j < i; ++j) {
                result(j, i) = result(j, i) / (a_std[i] * a_std[j]);
                result(i, j) = result(j, i);
            }

            result(i, i) = 1.0;
        }

    } else {
        std::cout << "corr don't support the method: " << method << " (Function: corr(Matrix a, ...))"<< std::endl;
        exit(0);
    }

    return;
}

// compute the correlation between matrix a and b
void corr(MatrixXd& a, MatrixXd& b, std::string method, MatrixXd& result) {
    int n = a.rows();
    int p = a.cols();
    int q = b.cols();
    if (n != b.rows()) {
        std::cout << "the sample size of matrix a and b is not equal ! (Function: corr(Matrix a, Matrix b, ...))" << std::endl;

        exit(0);
    }

    if (method == "spearman") {
        // convert matrix a and b into their index matrix
        MatrixXd a_index(n, p);
        convertMatrixValue2Index(a, a_index);
        MatrixXd b_index(n, q);
        convertMatrixValue2Index(b, b_index);
        // compute the pearson correlation between the index matrix a and b
        corr(a_index, b_index, "pearson", result);
    } else if (method == "pearson") {
        computeCov(a, b, result);
        VectorXd a_std(p);
        VectorXd b_std(q);
        for (int i = 0; i < p; ++i) {
            a_std[i] = sqrt(computeVar(a.col(i)));
            if (isZero(a_std[i] * a_std[i])) {
                std::cout << "Find zero variance in matrix a! (Function: corr(Matrix a, Matrix b, 'pearson', ...))" << std::endl;

                exit(0);
            }
        }
        for (int i = 0; i < q; ++i) {
            b_std[i] = sqrt(computeVar(b.col(i)));
            if (isZero(b_std[i] * b_std[i])) {
                std::cout << "Find zero variance in matrix b! (Function: corr(Matrix a, Matrix b, 'pearson', ...))" << std::endl;

                exit(0);
            }
        }
        for (int i = 0; i < p; ++ i) {
            for (int j = 0; j < q; ++j) {
                result(i, j) /= (a_std[i] * b_std[j]);
            }
        }

    } else {
        std::cout << "corr don't support the method: " << method << " (Function: corr(Matrix a, Matrix b, ...))"<< std::endl;
        exit(0);
    }
}

// construct n values list with the same interval from a to b
void seq(double a, double b, int num, std::vector<double>& sequences) {
    if (num == 1) {
        sequences.push_back(a);

        return;
    }
    /*if (num < 2) {
        std::cout << "the length of sequence cann't be less than 2 (Function: seq(a, b, num, ...))" << std::endl;
        exit(0);
    }*/

    double delta = (b - a) / (num - 1);
    for (int i = 0; i < num; ++i) {
        sequences.push_back(delta * i + a);
    }

    return;
}

// compute the union of two vector <int>
void getUnionList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& merged_list) {
    std::map<int, bool> value_map;
    for (int i = 0; i < list_a.size(); ++i) {
        value_map.insert(std::map<int, bool>::value_type (list_a[i], true));
    }
    for (int i = 0; i < list_b.size(); ++i) {
        value_map.insert(std::map<int, bool>::value_type (list_b[i], true));
    }

    for (std::map<int, bool>::iterator iter = value_map.begin(); iter != value_map.end(); iter++) {
        merged_list.push_back(iter->first);
    }
}

// compute the difference set of two integer sets
// exist in list_a but not in list_b
void getDiffList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& diff_list) {
    std::map<int, bool> value_map;
    for (int i = 0; i < list_b.size(); ++i) {
        value_map.insert(std::map<int, bool>::value_type (list_b[i], true));
    }
    for (int i = 0; i < list_a.size(); ++i) {
        std::map<int, bool>::iterator it;
        it = value_map.find(list_a[i]);
        if (it == value_map.end()) {
            diff_list.push_back(list_a[i]);
        }
    }
}

// compute the intersection set of two integer sets
// exist both in list_a and in list_b
void getIntersectList(std::vector<int>& list_a, std::vector<int>& list_b, std::vector<int>& intersect_list) {
    std::map<int, bool> value_map;
    for (int i = 0; i < list_b.size(); ++i) {
        value_map.insert(std::map<int, bool>::value_type (list_b[i], true));
    }
    for (int i = 0; i < list_a.size(); ++i) {
        std::map<int, bool>::iterator it;
        it = value_map.find(list_a[i]);
        if (it != value_map.end()) {
            intersect_list.push_back(list_a[i]);
        }
    }
}


// remove the i-th row of matrix X, the first row is 0-th row
void removeMatrixRow(MatrixXd& X, int i) {
    int n = X.rows();
    if (n <= 1) {
        std::cout << "The matrix X has only one row !!! in function: removeMatrixRow(X, i)" << std::endl;
        exit(0);
    }
    if (i >= n || i < 0) {
        std::cout << "The index i exceed the row number of matrix X !!! in function: removeMatrixRow(X, i)" << std::endl;
        exit(0);
    }
    MatrixXd X_new = MatrixXd::Zero(n - 1, X.cols());

    if (i > 0) {
        X_new.topRows(i) = X.topRows(i);
    }
    if (i < n - 1) {
        X_new.bottomRows(n - i - 1) = X.bottomRows(n - i - 1);
    }

    X = X_new;

    return;
}

void safe_mkdir(std::string dir_name) {
    int flag = mkdir(dir_name.c_str(), 0775);
    if (flag == 0) {
        std::cout << "create directory success ~" << std::endl;
    } else {
        std::cout << "create directory fail !!!" << std::endl;
        exit(0);
    }
    return;
}


