#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>
#include <chrono>

#define _ERROR_ANALYSIS_
#define _DEBUG_PRINT_
#include "../../include/interp_poly.h"

/*--- test cases ---*/
static void test_function01(VecDoub_I &x, VecDoub_IO &y, double p) 
{
    for (int k=0; k<x.size(); k++) { y[k] = std::pow(x[k]-2, p); }
}
static void test_function02(VecDoub_I &x, VecDoub_IO &y) 
{
    for (int k=0; k<x.size(); k++) { y[k] = 1 / (1+x[k]*x[k]); }
}

/*--- randomly pregenerated grid points ---*/
static double xgrid[20] = { -4.0, -2.76636, -2.72608, -2.27393, -2.0523, -1.98908, -1.37064, 
                            -0.60345, -0.270115, -0.256095, 0.7301, 0.78089, 0.919678, 
                            2.16769, 2.2435, 2.80859, 3.04971, 3.34071, 3.90351, 4.0 };

/*--- recording test results ---*/
static void inline record_results(
    std::string file,  std::ios_base::openmode mode, 
    int deg, double time, VecDoub &v
    ) 
{
    std::ofstream testfile;
    testfile.open(file, mode);
    for (int k=0; k<v.size(); k++) {
        testfile << std::to_string(v[k]) << " ";
    }
    testfile << "\n";
    testfile << "polydeg: " << deg << " :: ";
    testfile << "time: " << time << " (ms)" << std::endl;
    testfile.close();
}
static void inline record_weights(std::string file,  std::ios_base::openmode mode, VecDoub &w) 
{
    std::ofstream testfile;
    testfile.open(file, mode);
    for (int k=0; k<w.size(); k++) {
        testfile << std::to_string(w[k]) << " ";
    }
    testfile << "\n";
    testfile.close();
}


/**
 * main :: test
 *  A suite of 5 test cases; 
 *  test case 01 -
 *  test case 02 -
 *  test case 03 -
 *  test case 04 -
 *  test case 05 -
 */
int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    VecDoub xval(20, xgrid);
    VecDoub yval1(20), yval2(20);
    test_function01(xval, yval1, 9);
    test_function02(xval, yval2);

    int polydeg[7] = {3, 4, 5, 8, 9, 12, 19}; // list of polynomial degree
    double times[7]; // list of test times

    // interpolated values of x,y;
    //  a uniform grid on [-4,4] with spacing h=0.2;
    VecDoub xtest(39), ytest(39);
    xtest[0]=-3.8; xtest[39]=3.8;
    for (int i=0; i<39; i++) {
        xtest[i] = xtest[0] + i*0.2;
    }

    // polynomial interpolation object;
    Poly_interp* polyinterp;
/*-----------------------------------------------------------------------/
/    main :: test case 01 - Lagrange barycentric form 02 interpolation
/-----------------------------------------------------------------------*/

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval1, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Lagrange_interp(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/LagrangeForm02_function01_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval2, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Lagrange_interp(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/LagrangeForm02_function02_results.txt", std::ios_base::app, deg, times[k], ytest);
    }


/*-----------------------------------------------------------------------/
/    main :: test case 01 - Lagrange barycentric form 01 interpolation
/-----------------------------------------------------------------------*/

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval1, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Lagrange_interp_alt(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/LagrangeForm01_function01_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval2, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Lagrange_interp_alt(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/LagrangeForm01_function02_results.txt", std::ios_base::app, deg, times[k], ytest);
    }


/*-----------------------------------------------------------------------/
/    main :: test case 02 - Newton Horner Rule interpolation
/-----------------------------------------------------------------------*/

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval1, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Newton_Horner_rule(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/NewtonHorner_function01_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval2, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Newton_Horner_rule(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/NewtonHorner_function02_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

/*-----------------------------------------------------------------------/
/    main :: test case 02 - Newton Neville scheme interpolation
/-----------------------------------------------------------------------*/

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval1, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Neville_scheme(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/NevilleScheme_function01_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

    for (int k=0; k<7; k++) {
        int deg = polydeg[k];
        polyinterp = new Poly_interp(xval, yval2, deg);
        auto t1 = high_resolution_clock::now();
        for (int i=0; i<39; i++) {
            int lowidx = polyinterp->locate(xtest[i]);
            ytest[i] = polyinterp->Neville_scheme(lowidx, xtest[i]);
        }
        auto t2 = high_resolution_clock::now();
        duration<double, std::milli> ms_double = t2 - t1; // get the number of milliseconds as a double
        times[k] = ms_double.count();
        delete polyinterp;
        record_results("./results/NevilleScheme_function02_results.txt", std::ios_base::app, deg, times[k], ytest);
    }

/*-----------------------------------------------------------------------/
/    main :: test case 03 - uniform mesh interpolation
/-----------------------------------------------------------------------*/

    VecDoub xv1=VecDoub(14), yv1=VecDoub(14), bv1=VecDoub(14);
    VecDoub xv2=VecDoub(13), yv2=VecDoub(13), bv2=VecDoub(13);
    Poly_interp::uniform_mesh(-4.0, 4.0, xv1);
    Poly_interp::uniform_mesh(-4.0, 4.0, xv2);
    test_function01(xv1, yv1, 9);
    test_function02(xv2, yv2);

    polyinterp = new Poly_interp(xv1, yv1, 13);
    (void) polyinterp->Lagrange_interp(0, 0.0);
    polyinterp->get_inverse_barycweights(bv1);
    delete polyinterp;

    polyinterp = new Poly_interp(xv2, yv2, 12);
    (void) polyinterp->Lagrange_interp(0, 0.0);
    polyinterp->get_inverse_barycweights(bv2);
    delete polyinterp;

    record_weights("./results/UniformMeshWeights_function01.txt", std::ios_base::app, bv1);
    record_weights("./results/UniformMeshWeights_function02.txt", std::ios_base::app, bv2);

/*-----------------------------------------------------------------------/
/    main :: test case 04 - Chebyshev interpolation
/-----------------------------------------------------------------------*/

    Poly_interp::Chebyshev_mesh(xv1);
    Poly_interp::Chebyshev_mesh(xv2);
    test_function01(xv1, yv1, 9);
    test_function02(xv2, yv2);

    polyinterp = new Poly_interp(xv1, yv1, 13);
    (void) polyinterp->Lagrange_interp(0, 0.0);
    polyinterp->get_inverse_barycweights(bv1);
    delete polyinterp;

    polyinterp = new Poly_interp(xv2, yv2, 12);
    (void) polyinterp->Lagrange_interp(0, 0.0);
    polyinterp->get_inverse_barycweights(bv2);
    delete polyinterp;

    record_weights("./results/ChebyshevMeshWeights_function01.txt", std::ios_base::app, bv1);
    record_weights("./results/ChebyshevMeshWeights_function02.txt", std::ios_base::app, bv2);

/*-----------------------------------------------------------------------/
/    main :: test case 05 - Vandermonde interpolation
/-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------/
/    main :: ~ clear memory ~
/-----------------------------------------------------------------------*/
    
}
