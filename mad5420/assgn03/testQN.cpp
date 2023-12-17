#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>

#define _ERROR_ANALYSIS_

/*--- global variables for error analysis ---*/
#ifdef _ERROR_ANALYSIS_
#ifdef _ERROR_ANALYSIS_RESIDUALS_
double* res_errs=nullptr;
#endif
double* step_errs;
double* fvals;
#endif

#include "../../include/quasinewton.h"

/*--- test functions ---*/
double Beck(VecDoub_I &x);
double TwoVarRosenbrock(VecDoub_I &x);
double NVarRosenbrock(VecDoub_I &x);
double LJpotential(VecDoub_I &x);

/**
 * --------------------//
 * static definitions
 * --------------------//
 */
static int MAXITER=1000;  // maximum descent iterations
static std::ofstream testfile;  // output file stream
static int TESTCASE;
/*--- clean up function ---*/
static void clear_errs()
{
    for (int i=0; i<MAXITER; i++) {
        step_errs[i]=0.0;
        fvals[i]=0.0;
    }
}

/*--- step errors write up ---*/
static void inline record_steperrs(std::string file,  std::ios_base::openmode mode)
{
    testfile.open(file, mode);
    for (int k=0; k<MAXITER; k++) {
        testfile << std::to_string(step_errs[k]) << " ";
    }
    testfile << "\n";
    testfile.close();
}

/*--- results write up ---*/
static void inline record_results(
    std::string file,  std::ios_base::openmode mode, 
    std::string testdetails, int iter, double time, VecDoub &x) 
{
    testfile.open(file, mode);
    testfile << testdetails << " :: ";
    testfile << "iterations: " << std::to_string(iter) << " :: ";
    testfile << "time: " << time << " (ms)" << " :: ";
    for (int k=0; k<x.size(); k++) {
        testfile << std::to_string(x[k]) << " ";
    }
    testfile << "\n";
    testfile.close();
}

/**
 * ----------------------------------------------//
 * main
 * Entry point to test cases.
 * All results are APPENDED to file.
 * AFFECT variable [ static int TESTCASE ] to
 *  fix the test case.
 * ----------------------------------------------//
 */
int main(int argc, char *argv[])
{
    int iter;
    double fp;
    line_search_params LSparams;
    LSparams.rho=0.350;  // rho
    LSparams.sigma=0.850;  // sigma
    LSparams.omega=0.0;  // omega
    LSparams.alpha0=1.0;
    LSparams.maxiter=MAXITER;

    static int TESTCASE=4;
    step_errs = new double[MAXITER];
    fvals = new double[MAXITER];
    objPhi testfunc;
    VecDoub p;
    MatDoub* P;
    int dimn[9] = {2,4,8,16,32,64,128,256,512};  // for test case 4

    switch(TESTCASE) {
    case 1:
/*------------------------------/
/    test 01
/------------------------------*/
    testfunc = Beck;

    LSparams.IES = WolfePowellS;
    p.assign(2,0.0);
    p[0]=5.0; p[1]=0.0;
    iter=0;
    BFGSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
    record_results("test01/BFGSresults.txt", std::ios_base::app, "test01 :: BFGS :: WolfePowell", iter, 0.0, p);
    record_steperrs("test01/BFGSsteperrs.txt", std::ios_base::app);
    clear_errs();

    LSparams.IES=WolfePowellS;
    p.assign(2,0.0);
    p[0]=1.0; p[1]=1.0;
    iter=0;
    TruncatedCGDupdate(p, testfunc, MAXITER, iter, fp, &LSparams, nullptr);
    record_results("test01/CGDresults.txt", std::ios_base::app, "test01 :: CGD :: NoP :: WolfePowellS", iter, 0.0, p);
    record_steperrs("test01/CGDsteperrs.txt", std::ios_base::app);
    clear_errs();

    LSparams.IES=WolfePowellS;
    p.assign(2,0.0);
    p[0]=5.0; p[1]=5.0;
    iter=0;
    GMSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
    record_results("test01/GMSresults.txt", std::ios_base::app, "test01 :: GMS :: WolfePowell", iter, 0.0, p);
    record_steperrs("test01/GMSsteperrs.txt", std::ios_base::app);
    clear_errs();

    break;
    case 2:
/*------------------------------/
/    test 02
/------------------------------*/
    testfunc = TwoVarRosenbrock;
    LSparams.IES = WolfePowellS;

    p.assign(2,7.0);
    iter=0;
    BFGSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
    record_results("test02/BFGSresults.txt", std::ios_base::app, "test02 :: BFGS :: WolfePowell", iter, 0.0, p);
    record_steperrs("test02/BFGSsteperrs.txt", std::ios_base::app);
    clear_errs();

    p.assign(2,7.0);
    iter=0;
    GMSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
    record_results("test02/GMSresults.txt", std::ios_base::app, "test02 :: GMS :: WolfePowell", iter, 0.0, p);
    record_steperrs("test02/GMSsteperrs.txt", std::ios_base::app);
    clear_errs();

    p.assign(2,2.0);
    iter=0;
    TruncatedCGDupdate(p, testfunc, MAXITER, iter, fp, &LSparams, nullptr);
    record_results("test02/CGDresults.txt", std::ios_base::app, "test02 :: CGD :: NoP :: WolfePowell", iter, 0.0, p);
    record_steperrs("test02/CGDsteperrs.txt", std::ios_base::app);
    clear_errs();

    p.assign(2,2.0);
    iter=0;
    P = new MatDoub(2,2);
    (*P)[0][0]=12000.;   (*P)[0][1]=0.0;
    (*P)[1][0]=0.0;     (*P)[1][1]=400.;
    TruncatedCGDupdate(p, testfunc, MAXITER, iter, fp, &LSparams, P);
    record_results("test02/CGDresults.txt", std::ios_base::app, "test02 :: CGD :: P :: WolfePowell", iter, 0.0, p);
    record_steperrs("test02/CGDsteperrs.txt", std::ios_base::app);
    clear_errs();
    delete P;

    break;
    case 3:
/*------------------------------/
/    test 03
/------------------------------*/
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    testfunc = LJpotential;

    LSparams.IES = WolfePowellS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        auto tic = high_resolution_clock::now();
        BFGSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test03/BFGSresults.txt", std::ios_base::app, 
            "test03 :: BFGS :: " + to_string(dimn[k]) + " :: WolfePowell", 
            iter, ms_double.count(), p);
        record_steperrs("test03/BFGSsteperrs.txt", std::ios_base::app);
        clear_errs();
    }

    LSparams.IES=StrongWolfeS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        P = new MatDoub(dimn[k],dimn[k],0.0);
        for (int i=0; i<dimn[k]; i++) {
            (*P)[i][i]=1/(dimn[k]*dimn[k]);
        }
        auto tic = high_resolution_clock::now();
        TruncatedCGDupdate(p, testfunc, MAXITER, iter, fp, &LSparams, nullptr);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test03/CGDresults.txt", std::ios_base::app, 
            "test03 :: CGD :: " + to_string(dimn[k]) + " :: StrongWolfe", 
            iter, ms_double.count(), p);
        record_steperrs("test03/CGDsteperrs.txt", std::ios_base::app);
        clear_errs();
        delete P;
    }

    LSparams.IES=StrongWolfeS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        auto tic = high_resolution_clock::now();
        GMSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test03/GMSresults.txt", std::ios_base::app, 
            "test03 :: GMS :: " + to_string(dimn[k]) + " :: StrongWolfe", 
            iter, ms_double.count(), p);
        record_steperrs("test03/GMSsteperrs.txt", std::ios_base::app);
        clear_errs();
    }

    break;
    case 4:
/*------------------------------/
/    test 04
/------------------------------*/
    testfunc = NVarRosenbrock;

    LSparams.IES = WolfePowellS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        auto tic = high_resolution_clock::now();
        BFGSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test04/BFGSresults.txt", std::ios_base::app, 
            "test04 :: BFGS :: " + to_string(dimn[k]) + " :: WolfePowell", 
            iter, ms_double.count(), p);
        record_steperrs("test04/BFGSsteperrs.txt", std::ios_base::app);
        clear_errs();
    }

    LSparams.IES=StrongWolfeS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        auto tic = high_resolution_clock::now();
        TruncatedCGDupdate(p, testfunc, MAXITER, iter, fp, &LSparams, nullptr);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test04/CGDresults.txt", std::ios_base::app, 
            "test04 :: CGD :: " + to_string(dimn[k]) + " :: StrongWolfe", 
            iter, ms_double.count(), p);
        record_steperrs("test04/CGDsteperrs.txt", std::ios_base::app);
        clear_errs();
    }

    LSparams.IES=StrongWolfeS;
    for (int k=0; k<9; k++) {
        p.assign(dimn[k],15.0);
        iter=0;
        auto tic = high_resolution_clock::now();
        GMSupdate(p, testfunc, MAXITER, iter, fp, &LSparams);
        auto toc = high_resolution_clock::now();
        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;
        record_results("test04/GMSresults.txt", std::ios_base::app, 
            "test04 :: GMS :: " + to_string(dimn[k]) + " :: StrongWolfe", 
            iter, ms_double.count(), p);
        record_steperrs("test04/GMSsteperrs.txt", std::ios_base::app);
        clear_errs();
    }    
    
    default:
/*------------------------------/
/    default case
/------------------------------*/
    break;
    }  // end switch

    // ~ clear memory ~ //
    delete[] step_errs;
    delete[] fvals;
    // P should already be cleared;
}


/*------------------------------/
/    test functions
/------------------------------*/

double Beck(VecDoub_I &x)
{
    using namespace std;
    return ( sqrt(x[0]*x[0]+1) + sqrt(x[1]*x[1]+1) );
}

double TwoVarRosenbrock(VecDoub_I &x)
{
    using namespace std;
    return ( 100*pow(x[1]-x[0]*x[0], 2) + pow(1-x[0],2) );
}

double NVarRosenbrock(VecDoub_I &x)
{
    using namespace std;
    int n=x.size(), N;
    N = 2*(n/2);
    double sum=0.0;
    for (int i=0; i<N; i+=2) {
        sum += 100*pow(x[i+1]-x[i]*x[i],2) + pow(1-x[i],2);
    }
    return (sum);
}

double LJpotential(VecDoub_I &x)
{
    using namespace std;
    double r=VECNORM(x,2), sigma=7.0, welldepth=4.0;
    return ( 4 * welldepth * ( pow(sigma/r,12) - pow(sigma/r,6) ) );
}