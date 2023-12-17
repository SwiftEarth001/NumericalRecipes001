#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

// global vars for error analysis;
double* res_errs;
double* step_errs;
double* fvals;
static int MAXITER=250;

#define _ERROR_ANALYSIS_
#define _DEBUG_PRINT_
#include "../../../include/projdesc.h"

/*--- test cases ---*/
// test01 functions;
double objective_function001(VecDoub_I &x);
double ineq_constraint11(VecDoub_I &x);
double ineq_constraint12(VecDoub_I &x);

// test02 functions;
double objective_function002(VecDoub_I &x);
double ineq_constraint23(VecDoub_I &x);
double ineq_constraint25(VecDoub_I &x);
double ineq_constraint26(VecDoub_I &x);
double ineq_constraint27(VecDoub_I &x);

// test03 functions;
double objective_function003(VecDoub_I &x);
double ineq_constraint31(VecDoub_I &x);
double ineq_constraint32(VecDoub_I &x);
double ineq_constraint33(VecDoub_I &x);

// test04 functions;
double objective_function004(VecDoub_I &x);
double ineq_constraint41(VecDoub_I &x);

// test05 functions;
double objective_function005(VecDoub_I &x);
double ineq_constraint51(VecDoub_I &x);
double ineq_constraint52(VecDoub_I &x);
double ineq_constraint53(VecDoub_I &x);
double ineq_constraint54(VecDoub_I &x);
double ineq_constraint55(VecDoub_I &x);
double ineq_constraint56(VecDoub_I &x);


/*--- static definitions ---*/
// results files
static std::string output_file[5] = {
    "results\\test01_results.txt",
    "results\\test02_results.txt",
    "results\\test03_results.txt",
    "results\\test04_results.txt",
    "results\\test05_results.txt" 
};
static std::ofstream testfile;

// for writing to results files;
static void inline record_results(int test, int iter, VecDoub &x)
{
    testfile.open( output_file[test-1] );
    for (int k=0; k<x.size(); k++) {
        testfile << std::to_string(x[k]) << " ";
    }
    testfile << "\n";
    testfile << "iterations: " << std::to_string(iter) << "\n";
    for (int k=0; k<MAXITER; k++) {
        testfile << std::to_string(res_errs[k]) << " ";
    }
    testfile << "\n";
    for (int k=0; k<MAXITER; k++) {
        testfile << std::to_string(step_errs[k]) << " ";
    }
    testfile << "\n";
    for (int k=0; k<MAXITER; k++) {
        testfile << std::to_string(fvals[k]) << " ";
    }
    testfile << "\n";
    testfile.close();
}

// for clearing recorded errors and residuals;
static void inline clear_errvals() 
{
    for (int k=0; k<MAXITER; k++) {
        res_errs[k]=0.0;
        step_errs[k]=0.0;
        fvals[k]=0.0;
    }
}


/*--- main method ---*/
int main(int argc, char *argv[])
{
    int iter;
    double stpmax=1.0;
    double fp;

    res_errs = new double[MAXITER];
    step_errs = new double[MAXITER];
    fvals = new double[MAXITER];
    clear_errvals();

/*------------------------------/
/    test case 01
/------------------------------*/
    objPhi test01 = objective_function001;
    cnstrPhi cnstrnt01[2];
    cnstrnt01[0]=ineq_constraint11;
    cnstrnt01[1]=ineq_constraint12;
    VecDoub x01(2,0.0);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "starting test 01...\n");
#endif

    ProjGradDesc(x01, test01, cnstrnt01, 0, 2, stpmax, iter, fp);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "finished test 01...\n\n");
#endif

    record_results(1, iter, x01);
    clear_errvals();

/*------------------------------/
/    test case 02
/------------------------------*/
    objPhi test02 = objective_function002;
    cnstrPhi cnstrnt02[4];
    cnstrnt02[0]=ineq_constraint23;
    cnstrnt02[1]=ineq_constraint25;
    cnstrnt02[2]=ineq_constraint26;
    cnstrnt02[3]=ineq_constraint27;
    VecDoub x02(2,1.0);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "starting test 02...\n");
#endif

    ProjGradDesc(x02, test02, cnstrnt02, 0, 4, stpmax, iter, fp);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "finished test 02...\n\n");
#endif

    record_results(2, iter, x02);
    clear_errvals();

/*------------------------------/
/    test case 03
/------------------------------*/
    objPhi test03 = objective_function003;
    cnstrPhi cnstrnt03[3];
    cnstrnt03[0]=ineq_constraint31;
    cnstrnt03[1]=ineq_constraint32;
    cnstrnt03[2]=ineq_constraint33;
    VecDoub x03(2,1.0);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "starting test 03...\n");
#endif

    ProjGradDesc(x03, test03, cnstrnt03, 0, 3, stpmax, iter, fp);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "finished test 03...\n\n");
#endif

    record_results(3, iter, x03);
    clear_errvals();

/*------------------------------/
/    test case 04
/------------------------------*/
    objPhi test04 = objective_function004;
    cnstrPhi cnstrnt04[1];
    cnstrnt04[0]=ineq_constraint41;
    VecDoub x04(2,2.0);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "starting test 04...\n");
#endif

    ProjGradDesc(x04, test04, cnstrnt04, 1, 0, stpmax, iter, fp);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "finished test 04...\n\n");
#endif

    record_results(4, iter, x04);
    clear_errvals();

/*------------------------------/
/    test case 05
/------------------------------*/
    objPhi test05 = objective_function005;
    cnstrPhi cnstrnt05[6];
    cnstrnt05[0]=ineq_constraint51;
    cnstrnt05[1]=ineq_constraint52;
    cnstrnt05[2]=ineq_constraint53;
    cnstrnt05[3]=ineq_constraint54;
    cnstrnt05[4]=ineq_constraint55;
    cnstrnt05[5]=ineq_constraint56;
    VecDoub x05(5,0.5);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "starting test 05...\n");
#endif

    ProjGradDesc(x05, test05, cnstrnt05, 1, 0, stpmax, iter, fp);

#ifdef _DEBUG_PRINT_
    fprintf(stderr, "finished test 05..\n\n");
#endif

    record_results(5, iter, x05);
    clear_errvals();

/*------------------------------/
/    ~ clear memory ~
/------------------------------*/
    delete[] step_errs;
    delete[] res_errs;
    delete[] fvals;
}



/*------------------------------/
/    test 01 functions
/------------------------------*/

double objective_function001(VecDoub_I &x)
{
    using namespace std;
    return ( pow(x[0]+1,2) + pow(x[1]+1,2) );
}

double ineq_constraint11(VecDoub_I &x)
{
    using namespace std;
    return ( 2 - pow(x[0],2) - pow(x[1],2) );
}

double ineq_constraint12(VecDoub_I &x)
{
    using namespace std;
    return ( 1-x[1] );
}


/*------------------------------/
/    test 02 functions
/------------------------------*/

double objective_function002(VecDoub_I &x)
{
    using namespace std;
    return ( -1000*x[0] - 1000*x[1] + pow(x[0],2) + pow(x[1],2) );
}

double ineq_constraint23(VecDoub_I &x)
{
    using namespace std;
    return ( 3*x[0]+x[1]-3 );
}

double ineq_constraint25(VecDoub_I &x)
{
    using namespace std;
    return ( x[0] );
}

double ineq_constraint26(VecDoub_I &x)
{
    using namespace std;
    return ( x[1] );
}

double ineq_constraint27(VecDoub_I &x)
{
    using namespace std;
    return ( -x[0]-x[1]+10 );
}


/*------------------------------/
/    test 03 functions
/------------------------------*/

double objective_function003(VecDoub_I &x)
{
    using namespace std;
    return ( -pow(x[0],2) - pow(x[1],2) );
}

double ineq_constraint31(VecDoub_I &x)
{
    using namespace std;
    return ( 8.0-x[0] );
}

double ineq_constraint32(VecDoub_I &x)
{
    using namespace std;
    return ( 8.0-x[1] );
}

double ineq_constraint33(VecDoub_I &x)
{
    using namespace std;
    return ( x[0]+x[1]-1 );
}


/*------------------------------/
/    test 04 functions
/------------------------------*/

double objective_function004(VecDoub_I &x)
{
    using namespace std;
    return ( -2*x[0]+x[1] );
}

double ineq_constraint41(VecDoub_I &x)
{
    using namespace std;
    return ( x[1]-pow(x[0],2) );
}


/*------------------------------/
/    test 05 functions
/------------------------------*/

double objective_function005(VecDoub_I &x)
{
    using namespace std;
    double r=VECNORM(x,2);
    return ( exp(-0.5*r) / (1 + exp(-0.5*r)) );
}

double ineq_constraint51(VecDoub_I &x)
{
    using namespace std;
    double r=VECNORM(x,1);
    return ( 1-r );
}

double ineq_constraint52(VecDoub_I &x)
{
    using namespace std;
    return ( x[0] );
}

double ineq_constraint53(VecDoub_I &x)
{
    using namespace std;
    return ( x[1] );
}

double ineq_constraint54(VecDoub_I &x)
{
    using namespace std;
    return ( x[2] );
}

double ineq_constraint55(VecDoub_I &x)
{
    using namespace std;
    return ( x[3] );
}

double ineq_constraint56(VecDoub_I &x)
{
    using namespace std;
    return ( x[4] );
}

