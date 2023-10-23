#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../../../include/richardson.h"

extern void Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax);

extern void JOR_Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax);

extern void GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err, 
    const Int norm, const Doub tol, const Int itmax);

extern void SOR_GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax);




int main(int argc, char *argv[])
{
    std::string floats_file[4];
    floats_file[0] = "..\\data\\sample1.txt";  // 25 x 25
    floats_file[1] = "..\\data\\sample2.txt";  // 50 x 50
    floats_file[2] = "..\\data\\sample3.txt";  // 50 x 50
    floats_file[3] = "..\\data\\sample4.txt";  // 25 x 25

    std::string output_file[4];
    output_file[0] = "results\\sample01_results.txt";
    output_file[1] = "results\\sample02_results.txt";
    output_file[2] = "results\\sample03_results.txt";
    output_file[3] = "results\\sample04_results.txt";

    int sizes[4] = { 25, 50, 50, 25 };

    for (int m=0; m<4; m++) {
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
        MatDoub* A = new MatDoub(N, N, randfloats);
        VecDoub* b = new VecDoub(N, bvec);
        VecDoub* x = new VecDoub(N);
        int iter=0, itmax=5000;
        double err;
		int norm=2;  // l2 norm;
        double tol = 10e-7;
        GaussSeidel(*A, *b, *x, &iter, &err, 2, tol, itmax);

        // write results to file;
        std::ofstream test05file;
        test05file.open( output_file[m] );
        for (int k=0; k<N; k++) {
            test05file << std::to_string((*b)[k]) << " ";
        }
        test05file << "\n";
        for (int k=0; k<N; k++) {
            test05file << std::to_string((*x)[k]) << " ";
        }
        test05file << "\n";
        test05file << "iterations: " << std::to_string(iter) << "\n";
        test05file << "error: " << std::to_string(err) << "\n";
        test05file.close();

        delete[] randfloats;
        delete[] bvec;
        delete b;
        delete x;
    }
}