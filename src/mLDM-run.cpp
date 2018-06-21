#include "utils.h"
#include "mldm.h"
#include "mldm-openmp.cpp"
#include <fstream>

#define SMALL_NOISE 1e-4
#define MAX_CHARS 1000

int main(int argc, char* argv[]) {
    // google::InitGoogleLogging(argv[0]);
    // FLAGS_logtostderr = 1;
    // LOG(INFO) << "test glog";
    
    char otu_file[MAX_CHARS];
    char meta_file[MAX_CHARS];
    char shape_file[MAX_CHARS];
    char result_file[MAX_CHARS];

    if (argc < 2) {
        std::cout << "Please specify the path of otu table !" << std::endl;
        return 0;
    }
    if (argc < 3) {
        std::cout << "Please specify the path of meta table !" << std::endl;
        return 0;
    }
    if (argc < 4) {
        std::cout << "Please specify the path of matrix shape !" << std::endl;
        return 0;
    }
    if (argc < 5) {
        std::cout << "Please specify the path of result file !" << std::endl;
        return 0;
    }

    strcpy(otu_file, argv[1]);
    strcpy(meta_file, argv[2]);
    strcpy(shape_file, argv[3]);
    strcpy(result_file, argv[4]);
    
    int N, P, Q;
    // read shape
    std::ifstream shape_in;
    shape_in.open(shape_file);
    shape_in >> N >> P >> Q;
    shape_in.close();
    MatrixXi X = MatrixXi::Zero(N, P);
    MatrixXd M = MatrixXd::Zero(N, Q);
    // read otu
    std::ifstream otu_in;
    otu_in.open(otu_file);
    int i,j;
    for (i = 0; i < N; ++i) {
        for (j = 0; j < P; ++j) {
            otu_in >> X(i, j);
        }
    }
    otu_in.close();
    // read meta
    std::ifstream meta_in;
    meta_in.open(meta_file);
    for (i = 0; i < N; ++i) {
        for (j = 0; j < Q; ++j) {
            meta_in >> M(i, j);
        }
    }
    for (i = 0; i < Q; ++i) {
        M.col(i) += (VectorXd::Random(N, 1) * SMALL_NOISE);
    }
    meta_in.close();
    // mldm inference
    int max_iters = 1500;
    double thres = 1e-4;
    bool verb = true;
    double m_start, m_end;

    // convert integer to double
    MatrixXd x_d = X.cast <double> ();
    std::vector<int> otu_index;
    double zero_ratio = 0.9;

    for (int i = 0 ;i < P; ++i) {
        double var_i = computeVar(x_d.col(i));
        int zero_i = zeroNum(x_d.col(i));

        if (isZero(var_i) || (zero_i * 1.0 / N > zero_ratio)) {
            continue;
        }

        otu_index.push_back(i);
    }

    int p = otu_index.size();

    // filter zero variance
    MatrixXi X_new = MatrixXi::Zero(N, p);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < p; ++j) {
            X_new(i, j) = X(i, otu_index[j]);
        }
    }

    m_start = omp_get_wtime();
    mLDM mldm_worker(X_new, M, N, p, Q, max_iters, thres, verb);
    mldm_worker.setModelSelectionNum(4);
    mldm_worker.work();

    m_end = omp_get_wtime();
    if (mldm_worker.best_result.exists) {
        MatrixXd B = mldm_worker.best_result.B;
        MatrixXd Theta = mldm_worker.best_result.Theta;
        std::string result_dir = result_file;
        std::string cluster_otu_otu_file = result_dir + "/Theta";
        std::string cluster_meta_otu_file = result_dir + "/B";
        std::string cluster_otu_list_file = result_dir + "/OTU_Index";
        // export meta-otu associations
        std::ofstream meta_otu_out(cluster_meta_otu_file);
        int q = B.rows();
        int p = B.cols();
        for (i = 0; i < q; ++i) {
            for (j = 0; j < p; ++j) {
                meta_otu_out << B(i, j);
                if (j != p - 1) {
                    meta_otu_out << " ";
                }
            }
            if (i != q - 1) {
                meta_otu_out << "\n";
            }
        }
        meta_otu_out.close();
        // otu otu association
        std::ofstream otu_otu_out(cluster_otu_otu_file);
        for (i = 0; i < p; ++i) {
            for (j = 0; j < p; ++j) {
                if (i != j) {
                    otu_otu_out << - Theta(i, j) / sqrt(Theta(i,i) * Theta(j,j));
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
        std::cout << "Export Result Success ~" << std::endl;

        // export otu list
        std::ofstream otu_out(cluster_otu_list_file);
        for (i = 0; i < p; ++ i) {
            otu_out << otu_index[i];
            if (i != p - 1) {
                otu_out << " ";
            }
        }
        otu_out.close();

        double obj = mldm_worker.best_result.obj;
        double EBIC = mldm_worker.best_result.EBIC;
    } else {
        std::cout << "Error when estimate parameters " << std::endl;
        return 0;
    }

    return 0;
}
