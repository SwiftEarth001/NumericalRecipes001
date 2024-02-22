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
static void test_function01(VecDoub_I &x, VecDoub_IO &y) 
{
    for (int k=0; k<x.size(); k++) { y[k] = std::pow(x[k]-2, 9); }
}
static void test_function02(VecDoub_I &x, VecDoub_IO &y) 
{
    for (int k=0; k<x.size(); k++) { y[k] = 1 / (1+x[k]*x[k]); }
}

/*--- randomly pregenerated grid points ---*/
static double xgrid[25] = { -4.95755616, -3.80701176, -3.76215901, -3.70701338, -2.72490868, -2.66172848, -2.46249911, -2.42192397, -1.52695176, -0.50615107, -0.34544951,  0.08156321,  0.25824084,  0.41845151,  0.88614905, 1.65196149,  1.75338246,  2.10957078,  2.30142646,  2.43887168, 2.49140943,  2.52375555,  2.88355883,  3.71524474,  4.45083864  };
static double f1[4][25] =  { { -38204144.849, -7508865.961, -7002726.187, -6422143.685, -1173655.861, -1039735.252, -701783.715, -646400.142, -84449.228, -3900.006, -2148.003, -351.988, -147.534, -61.907, -2.639, -0.0, -0.0, 0.0, 0.0, 0.001, 0.002, 0.003, 0.328, 128.509, 3190.307 },
                             { 49419263.878, 11637619.554, 10937659.933, 10127765.487, 2235578.181, 2007327.821, 1415362.396, 1315626.707, 215495.731, 14005.56, 8242.353, 1651.287, 762.338, 352.292, 21.324, 0.002, 0.0, 0.0, 0.001, 0.012, 0.031, 0.051, 3.343, 674.297, 11715.485 }, 
                             { -56823703.918, -16032506.947, -15185502.41, -14196939.533, -3785179.071, -3444778.616, -2537344.855, -2380188.741, -488797.682, -44707.792, -28113.513, -6885.967, -3501.46, -1782.01, -153.152, -0.045, -0.004, 0.0, 0.016, 0.226, 0.498, 0.778, 30.268, 3144.96, 38241.554 }, 
                             { 57170350.98, 19326213.414, 18447688.925, 17413412.25, 5607781.081, 5172641.529, 3980149.582, 3767889.562, 970124.913, 124874.573, 83904.851, 25125.545, 14072.108, 7887.25, 962.483, 0.896, 0.113, 0.001, 0.378, 3.601, 7.097, 10.404, 239.797, 12834.741, 109224.194 } };
static double f2[2][25] = { { 0.039, 0.065, 0.066, 0.068, 0.119, 0.124, 0.142, 0.146, 0.3, 0.796, 0.893, 0.993, 0.937, 0.851, 0.56, 0.268, 0.245, 0.183, 0.159, 0.144, 0.139, 0.136, 0.107, 0.068, 0.048 },
                            { 6486.481, 1827.697, 1727.876, 1611.254, 386.843, 347.962, 245.751, 228.33, 33.897, 1.597, 0.866, -0.165, -0.588, -1.156, -5.649, -45.942, -58.213, -125.332, -182.488, -235.478, -258.821, -274.116, -500.399, -1628.244, -3854.913 } };

/*--- recording test results ---*/
static void inline record_results(
    std::string file,  std::string comments, 
    std::ios_base::openmode mode, VecDoub &v
    ) 
{
    std::ofstream testfile;
    testfile.open(file, mode);
    for (int k=0; k<v.size(); k++) {
        testfile << std::to_string(v[k]) << " ";
    }
    testfile << "\n";
    testfile << comments << std::endl;
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
 *  A suite of 2 test cases and a base case to compare against; 
 *  test case 01 - Hermite-Newton Interpolation
 *  test case 02 - Cubic Spline Interpolation
 *  base case 03 - Newton-Neville Scheme & Lagrange Interpolations
 */
int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    typedef std::chrono::steady_clock::time_point tictoc;
    tictoc t1, t2;

    // data - ordinates and abscissa;
    VecDoub xval(25, xgrid);
    // y values for test function 01
    MatDoub ydat(4, 25, 0.0);
    for (int r=0; r<4; r++) {
        for (int j=0; j<25; j++) {
            ydat[r][j] = f1[r][j];
        }
    }
    // y values for test function 02
    MatDoub ydat20(2, 25, 0.0);
    for (int r=0; r<2; r++) {
        for (int j=0; j<25; j++) {
            ydat20[r][j] = f2[r][j];
        }
    }

    // for the test cases
    struct Hermite_interp* HM;
    struct Spline_interp* CS;
    int lowidx;
    std::string comments;
    duration<double, std::milli> ms_double;
    // for the base cases
    NewtonPolyInterp* NewtPoly;
    LagrangePolyInterp* LagrPoly;
    // for test case 01, 
    //  a uniform grid on [-4,4] with spacing h=0.2;
    VecDoub xtest(41), ytest(41);
    Poly_Base_Interp::uniform_mesh(-4.0, 4.0, xtest);

    int test=4;
    switch(test) {
    
/*----------------------------------------------------------------------------/
/    main :: test case 01 - Hermite-Newton interpolation
/----------------------------------------------------------------------------*/
    case 1:
    {
        int numpoints[4] = {4, 3, 5, 10};
        int difforder[4] = {2, 3, 1, 0};
        // test_function01
        for (int k=0; k<4; k++) {
            HM = new Hermite_interp(xval, ydat, numpoints[k], difforder[k]);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = HM->locate(xtest[i]);
                ytest[i] = HM->Newton_Hermite_rule(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "numpoints: " + std::to_string(numpoints[k]) + " :: " +
                        "difforder: " + std::to_string(difforder[k]) + " :: " +
                        "polydeg: " + std::to_string(HM->get_degree()) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/testcase01_function01_results.txt", 
                            comments, std::ios_base::app, ytest );
            delete HM;
        }
        int order[3] = {2, 3, 2};
        VecInt od(3, order);
        HM = new Hermite_interp(xval, ydat, od, 3);
        t1 = high_resolution_clock::now();
        for (int i=0; i<xtest.size(); i++) {
            lowidx = HM->locate(xtest[i]);
            ytest[i] = HM->Newton_Hermite_rule(lowidx, xtest[i]);
        }
        t2 = high_resolution_clock::now();
        ms_double = t2 - t1; // get the number of milliseconds as a double
        comments =  "numpoints: " + std::to_string(3) + " :: " + 
                    "orders: 2, 3, 2 :: " +
                    "polydeg: " + std::to_string(HM->get_degree()) + " :: " +
                    "time: " + std::to_string(ms_double.count()) + " (ms)";
        record_results( "./results/testcase01_function01_results.txt", 
                        comments, std::ios_base::app, ytest );
        delete HM;
    }
    break;

/*----------------------------------------------------------------------------/
/    main :: test case 02 - Cubic Spline interpolation
/----------------------------------------------------------------------------*/
    case 2:
    {
        int numpoints02[4] = {5, 11, 21, 31};
        VecDoub xpoints, ypoints;
        // test_function02
        for (int k=0; k<4; k++) {
            xpoints = VecDoub(numpoints02[k]);
            ypoints = VecDoub(numpoints02[k]);
            Poly_Base_Interp::uniform_mesh(-4.0, 4.0, xpoints);
            test_function02(xpoints, ypoints);
            CS = new Spline_interp(xpoints, ypoints, Spline_interp::EdgeConstraints::Natural);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = CS->locate(xtest[i]);
                ytest[i] = CS->interpolate(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "numpoints: " + std::to_string(numpoints02[k]) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/testcase02_function02_results.txt", 
                            comments, std::ios_base::app, ytest );
            delete CS;
        }
    }
    break;

/*----------------------------------------------------------------------------/
/    main :: base case 03 - Newton & Lagrange schemes - testfunction01
/----------------------------------------------------------------------------*/
    case 3:
    {
        int polydeg[2] = {9, 11};
        VecDoub yval(25, 0.0);
        for (int j=0; j<25; j++) { yval[j]=ydat[0][j]; }
        
        for (int deg : polydeg) {
            NewtPoly = new NewtonPolyInterp(xval, yval, deg);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = NewtPoly->locate(xtest[i]);
                ytest[i] = NewtPoly->Neville_scheme(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "polydeg: " + std::to_string(deg) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/basecase03_function01_NewtonNeville.txt", 
                            comments, std::ios_base::app, ytest );
            delete NewtPoly;
        }

        for (int deg : polydeg) {
            NewtPoly = new NewtonPolyInterp(xval, yval, deg);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = NewtPoly->locate(xtest[i]);
                ytest[i] = NewtPoly->Newton_Horner_rule(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "polydeg: " + std::to_string(deg) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/basecase03_function01_NewtonHorner.txt", 
                            comments, std::ios_base::app, ytest );
            delete NewtPoly;
        }

        for (int deg : polydeg) {
            LagrPoly = new LagrangePolyInterp(xval, yval, deg);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = LagrPoly->locate(xtest[i]);
                ytest[i] = LagrPoly->Lagrange_interp(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "polydeg: " + std::to_string(deg) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/basecase03_function01_Lagrange.txt", 
                            comments, std::ios_base::app, ytest );
            delete LagrPoly;
        }
    }
    break;

/*----------------------------------------------------------------------------/
/    main :: base case 04 - uniform Lagrange scheme - testfunction02
/----------------------------------------------------------------------------*/
    case 4:
    {
        int numpoints04[4] = {5, 11, 21, 31};
        VecDoub xpoints, ypoints;
        // test_function02
        for (int k=0; k<4; k++) {
            xpoints = VecDoub(numpoints04[k]);
            ypoints = VecDoub(numpoints04[k]);
            Poly_Base_Interp::uniform_mesh(-4.0, 4.0, xpoints);
            test_function02(xpoints, ypoints);
            LagrPoly = new LagrangePolyInterp(xpoints, ypoints, numpoints04[k]-1);
            t1 = high_resolution_clock::now();
            for (int i=0; i<xtest.size(); i++) {
                lowidx = LagrPoly->locate(xtest[i]);
                ytest[i] = LagrPoly->Lagrange_interp(lowidx, xtest[i]);
            }
            t2 = high_resolution_clock::now();
            ms_double = t2 - t1; // get the number of milliseconds as a double
            comments =  "numpoints: " + std::to_string(numpoints04[k]) + " :: " +
                        "time: " + std::to_string(ms_double.count()) + " (ms)";
            record_results( "./results/basecase04_function02_Lagrange.txt", 
                            comments, std::ios_base::app, ytest );
            delete LagrPoly;
        }
    }
    break;

/*----------------------------------------------------------------------------/
/    main :: default case
/----------------------------------------------------------------------------*/
    default:
        break;
    }  // fin switch;

/*----------------------------------------------------------------------------/
/    main :: ~ clear memory ~
/----------------------------------------------------------------------------*/
    
}
