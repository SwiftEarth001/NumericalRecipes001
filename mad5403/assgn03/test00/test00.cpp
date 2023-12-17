#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../../../include/linbcg.h"


int main(int argc, char *argv[])
{
    int N = 4;
    double* Atest = new double[16]{5, 7, 6, 5, 7, 10, 8, 7, 6, 8, 10, 9, 5, 7, 9, 10};
    double* bvec = new double[4]{57, 79, 88, 86};


    //---------------------------------------------
    // task 01 : without preconditioning
    //---------------------------------------------

    // put array to linbcg class;
    Linbcg* descentMethods = new Linbcg(N, N, Atest);
    VecDoub* b = new VecDoub(N, bvec);
    VecDoub* x = new VecDoub(N, bvec);
    int iter=0, itmax=5000;
    double reserr;
    double steperr;
    int norm=2;  // l2 norm
    double restol = 1.0e-08;
    double steptol = 1.0e-08;
    std::string output_file;  // filename to record results to
    std::ofstream test00file;  // output stream to write results to

    // <--- steepest descent SD <---
    descentMethods->GDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write steepest descent SD results to file;
    output_file = "results\\NOP_SD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- steepest descent SD >---

    // <--- conjugate gradient descent CGD <---
    for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
    descentMethods->CGDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write conjugate gradient descent CGD results to file;
    output_file = "results\\NOP_CGD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- conjugate gradient CGD >---
    


    //---------------------------------------------
    // task 02 : with Jacobi preconditioning
    //---------------------------------------------

    double* Ptest1 = new double[16]{5, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10, 0, 0, 0, 0, 10};
    MatDoub* PJacobi = new MatDoub(4, 4, Ptest1);
    descentMethods->assignPreconditioner(PJacobi);

    // <--- steepest descent SD <---
    for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
    descentMethods->GDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write steepest descent SD results to file;
    output_file = "results\\PJacobi_SD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- steepest descent SD >---

    // <--- conjugate gradient descent CGD <---
    for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
    descentMethods->CGDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write conjugate gradient descent CGD results to file;
    output_file = "results\\PJacobi_CGD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- conjugate gradient CGD >---

    delete Ptest1;
    delete PJacobi;



    //---------------------------------------------
    // task 03 : with symmetric 
    //      Gauss-Seidel preconditioning
    //---------------------------------------------

    double* Ptest2 = new double[16]{5, 7, 6, 5, 7, 19.8, 16.4, 14, 6, 16.4, 23.6, 20.6, 5, 14, 20.6, 28};
    MatDoub* PSGS = new MatDoub(4, 4, Ptest2);
    descentMethods->assignPreconditioner(PSGS);

    // <--- steepest descent SD <---
    for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
    descentMethods->GDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write steepest descent SD results to file;
    output_file = "results\\PSGS_SD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- steepest descent SD >---

    // <--- conjugate gradient descent CGD <---
    for (int k=0; k<N; k++) { (*x)[k] = (*b)[k]; }
    descentMethods->CGDsolve(*b, *x, &iter, &reserr, &steperr, 2, restol, steptol, itmax);

    // write conjugate gradient descent CGD results to file;
    output_file = "results\\PSGS_CGD_results.txt";
    test00file.open( output_file );
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*b)[k]) << " ";
    }
    test00file << "\n";
    for (int k=0; k<N; k++) {
        test00file << std::to_string((*x)[k]) << " ";
    }
    test00file << "\n";
    test00file << "iterations: " << std::to_string(iter) << "\n";
    test00file << "residual error: " << std::to_string(reserr) << "\n";
    test00file.close();
    // >--- conjugate gradient CGD >---

    delete Ptest2;
    delete PSGS;


    //---------------------------------------------
    // ~ clean up fin ~
    //---------------------------------------------

    delete[] Atest;
    delete[] bvec;
    delete descentMethods;
    delete b;
    delete x;
}