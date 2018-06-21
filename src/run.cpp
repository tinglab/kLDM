#include "kmldm.h"

#define MAX_CHARS 1000

int main(int argc, char* argv[]) {
    char otu_file[MAX_CHARS];
    char meta_file[MAX_CHARS];
    char shape_file[MAX_CHARS];
    char result_dir[MAX_CHARS];

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
        std::cout << "Please specify the directory's path of results !" << std::endl;
        return 0;
    }
    strcpy(otu_file, argv[1]);
    strcpy(meta_file, argv[2]);
    strcpy(shape_file, argv[3]);
    strcpy(result_dir, argv[4]);
   
    double k_begin = omp_get_wtime();
    KmLDM solver;
    // glog 输出信息
    // FLAGS_logtostderr = 1;
    // solver.enableRemoveSample(true);
    solver.loadData(otu_file, meta_file, shape_file);
    solver.setClusterThreshold(4 * solver.getOTUsNumber());
    // solver.setClusterThreshold(2 * solver.getOTUsNumber());
    solver.work();

    // save result into file
    solver.exportResult(result_dir);
    double k_end = omp_get_wtime();
    std::cout << "kmldm run time: " << k_end - k_begin << " seconds"<< std::endl;

    return 0;
}
