#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#define _ERROR_ANALYSIS_
static double* res_errs;

#include "../../../include/linbcg.h"


int main(int argc, char *argv[])
{
    res_errs = new double[5000];  // max iterations is 5000;

    std::string floats_file[7];
    floats_file[0] = "..\\data\\sample1spd.txt";    // 25 x 25
    floats_file[1] = "..\\data\\sample2spd.txt";    // 50 x 50
    floats_file[2] = "..\\data\\sample3spd.txt";    // 50 x 50
    floats_file[3] = "..\\data\\sample4spd.txt";    // 25 x 25
    floats_file[4] = "..\\data\\diff2banded10.txt"; // 10 x 10
    floats_file[5] = "..\\data\\diff2banded25.txt"; // 25 x 25
    floats_file[6] = "..\\data\\diff2banded50.txt"; // 50 x 50

    std::string SD_output_file[7];
    SD_output_file[0] = "results\\SD_sample01_results.txt";
    SD_output_file[1] = "results\\SD_sample02_results.txt";
    SD_output_file[2] = "results\\SD_sample03_results.txt";
    SD_output_file[3] = "results\\SD_sample04_results.txt";
    SD_output_file[4] = "results\\SD_diff2banded10_results.txt";
    SD_output_file[5] = "results\\SD_diff2banded25_results.txt";
    SD_output_file[6] = "results\\SD_diff2banded50_results.txt";

    std::string CGD_output_file[7];
    CGD_output_file[0] = "results\\CGD_sample01_results.txt";
    CGD_output_file[1] = "results\\CGD_sample02_results.txt";
    CGD_output_file[2] = "results\\CGD_sample03_results.txt";
    CGD_output_file[3] = "results\\CGD_sample04_results.txt";
    CGD_output_file[4] = "results\\CGD_diff2banded10_results.txt";
    CGD_output_file[5] = "results\\CGD_diff2banded25_results.txt";
    CGD_output_file[6] = "results\\CGD_diff2banded50_results.txt";

    int sizes[7] = { 25, 50, 50, 25, 10, 25, 50 };

    for (int m=0; m<7; m++) {
        // get data from file;
        std::ifstream inFile(floats_file[m]);
        std::string line;
        std::vector< std::vector<float> > all_floats;
        while ( getline(inFile, line) ) {
            std::istringstream is(line);
            all_floats.push_back(
                std::vector<float>( std::istream_iterator<float>(is),
                                    std::istream_iterator<float>() ) 
                );
        }
        inFile.close();

        // put data to array;
        int N = sizes[m];
        double* randfloats = new double[(N+1)*(N)];  // the last row is b vector;
        for (int i=0; i < N+1; i++) {
            for (int j=0; j < N; j++) {
                randfloats[i*N+j] = all_floats[i][j];
            }
        }
        double* bvec = new double[N];
        for (int j=0; j<N; j++) {
            bvec[j] = all_floats[N][j];
        }
        all_floats.clear();

        // put array to linbcg class;
        Linbcg* descentMethods = new Linbcg(N, N, randfloats);
        VecDoub* b = new VecDoub(N, bvec);
        VecDoub* x = new VecDoub(N, bvec);  // initialize x0 to b
        int iter=0, itmax=5000;
        double reserr;
        double steperr;
		int norm=2;  // l2 norm;
        double restol = 1.0e-08;
        double steptol = 1.0e-08;
        std::ofstream test01file;


        // <--- steepest descent SD <---
        descentMethods->GDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

        // write results to file;
        test01file.open( SD_output_file[m] );
        for (int k=0; k<N; k++) {
            test01file << std::to_string((*b)[k]) << " ";
        }
        test01file << "\n";
        for (int k=0; k<N; k++) {
            test01file << std::to_string((*x)[k]) << " ";
        }
        test01file << "\n";
        test01file << "iterations: " << std::to_string(iter) << "\n";
        test01file << "error: " << std::to_string(reserr) << "\n";
        for (int k=0; k<iter; k++) {
            test01file << std::to_string(res_errs[k]) << " ";
        }
        test01file << "\n";
        test01file.close();
        // >--- steepest descent SD >---


        // <--- conjugate gradient descent CGD <---
        for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
        descentMethods->CGDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

        // write results to file;
        test01file.open( CGD_output_file[m] );
        for (int k=0; k<N; k++) {
            test01file << std::to_string((*b)[k]) << " ";
        }
        test01file << "\n";
        for (int k=0; k<N; k++) {
            test01file << std::to_string((*x)[k]) << " ";
        }
        test01file << "\n";
        test01file << "iterations: " << std::to_string(iter) << "\n";
        test01file << "error: " << std::to_string(reserr) << "\n";
        for (int k=0; k<iter; k++) {
            test01file << std::to_string(res_errs[k]) << " ";
        }
        test01file << "\n";
        test01file.close();
        // >--- conjugate gradient descent CGD >---
        

        // ~ clean up;
        delete[] randfloats;
        delete[] bvec;
        delete descentMethods;
        delete b;
        delete x;
    }

    // ~ clean up residual errors;
    delete[] res_errs;
}