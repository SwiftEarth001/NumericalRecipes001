#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>

// global vars for error analysis;
double* res_errs;
double* step_errs;
double* fvals;
static int MAXITER=250;

#define _ERROR_ANALYSIS_
#include "../../../include/projdesc.h"

Doub LJPotential(VecDoub_I &x)
{
    using namespace std;
    double r=VECNORM(x,2), sigma=7.0, welldepth=4.0;
    return ( 4 * welldepth * ( pow(sigma/r,12) - pow(sigma/r,6) ) );
}

Doub ineq_constraint00(VecDoub_I &x)
{
    using namespace std;
    double r=VECNORM(x,2);
    return (r-17);
}

static void inline clear_errvals() 
{
    for (int k=0; k<MAXITER; k++) {
        res_errs[k]=0.0;
        step_errs[k]=0.0;
        fvals[k]=0.0;
    }
}

int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    double fp;
    double stpmax;
    int iter;

    res_errs = new double[MAXITER];
    step_errs = new double[MAXITER];
    fvals = new double[MAXITER];
    clear_errvals();

    double dimsize[8] = { 2, 3, 5, 10, 25, 50, 75, 100 };

    std::string output_file[8];
    output_file[0] = "results\\size02_results.txt";
    output_file[1] = "results\\size03_results.txt";
    output_file[2] = "results\\size05_results.txt";
    output_file[3] = "results\\size10_results.txt";
    output_file[4] = "results\\size25_results.txt";
    output_file[5] = "results\\size50_results.txt";
    output_file[6] = "results\\size75_results.txt";
    output_file[7] = "results\\size100_results.txt";
    std::ofstream test00file;


/*------------------------------/
/    test :: ...
/------------------------------*/
    objPhi test00 = LJPotential;
    cnstrPhi cnstrnt00[1];
    cnstrnt00[0]=ineq_constraint00;
    int n;  // stores the dimension size
    VecDoub* x0;

    for (int m=0; m<8; m++) {
        n=dimsize[m];
        stpmax = 1.0/sqrt(n);
        x0 = new VecDoub( n, 11.0/sqrt((double)(n)) );  // starting at a radius r=sqrt(11)

        auto tic = high_resolution_clock::now();
        ProjGradDesc(*x0, test00, cnstrnt00, 0, 1, stpmax, iter, fp);
        auto toc = high_resolution_clock::now();

        /* Getting number of milliseconds as a double. */
        duration<double, std::milli> ms_double = toc-tic;

        // write results to file;
        test00file.open( output_file[m] );
        for (int k=0; k<n; k++) {
            test00file << std::to_string((*x0)[k]) << " ";
        }
        test00file << "\n";
        test00file << "iterations: " << std::to_string(iter) << "\n";
        test00file << "time: " << ms_double.count() << " ms\n";
        for (int k=0; k<MAXITER; k++) {
            test00file << std::to_string(res_errs[k]) << " ";
        }
        test00file << "\n";
        for (int k=0; k<MAXITER; k++) {
            test00file << std::to_string(step_errs[k]) << " ";
        }
        test00file << "\n";
        for (int k=0; k<MAXITER; k++) {
            test00file << std::to_string(fvals[k]) << " ";
        }
        test00file << "\n";
        test00file.close();

        // clear vector x0;
        delete x0;
        // clear residuals, step-errors, and function values;
        clear_errvals();
    }

/*------------------------------/
/    ~ clear memory ~
/------------------------------*/
    delete[] step_errs;
    delete[] res_errs;
    delete[] fvals;
}