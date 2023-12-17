#ifndef _PROJDESC_H_
#define _PROJDESC_H_

/*------------------------------------------------------/
// we implement a suite of projected descent methods 
//   using the Numerical Recipes interface;
/------------------------------------------------------*/

//#include <float.h>  // already included in nr3math.h
//#include <cmath>  // already included in nr3.h

#include "nr3math.h"
#include "linsearch.h"
#include "qrdcmp.h"


/*--------------------------------------------------------/
// 1. static (and global) scope parameters
/--------------------------------------------------------*/

// starting parameters;
static int ITMAX=250;  // maximum descent iterations
static const double EPS32=(double)(numeric_limits<float>::epsilon());  // EPS32 is float32 precision
static double TOLX=1.0e-13, GTOL=1.0e-13, TOLY=1.0e-13;  // tolerances for descent-step, residual, and Newton-step errors
static const double TOLZERO=1.0e-13;  // tolerance for zero value
static int NMAX=50;  // maximum Newton iterations

#ifdef _ERROR_ANALYSIS_
extern double* res_errs;
extern double* step_errs;
extern double* fvals;
#endif



/*--------------------------------------------------------/
// 2. active constraint methods
/--------------------------------------------------------*/

static inline void evalConstrnts(VecDoub_I &x, VecDoub_O &cn, const cnstrPhi cnstrnt[], const int mm) 
{
    for (int i=0; i<mm; i++) { cn[i] = cnstrnt[i](x); }
}

static inline bool isFeasible(VecDoub_I &cn, const int mm01, const int mm)
{
    for (int i=0; i<mm01; i++) {
        if (abs(cn[i]) > TOLZERO) { return false; }
    }
    for (int i=mm01; i<mm; i++) {
        if (cn[i] < -TOLZERO) { return false; }
    }
    return true;
}

static inline void getActiveIndx(VecDoub_I &cn, int* indx, const int mm01, const int mm, int &ma) 
{
    ma=mm01;
    for (int i=0; i<mm; i++) {
        if (i<mm01) { indx[i]=1; }
        else if (cn[i] < TOLZERO) { 
            indx[i]=1;
            ma++;
        }
        else { indx[i]=0; }
    }
}

static inline void getActiveCn(VecDoub_I &cn, VecDoub_O &ca, int* indx, const int mm) 
{
    int ma=0;
    for (int i=0; i<mm; i++) {
        if (indx[i] == 1) { 
            ca[ma]=cn[i];
            ma++;
        }
    }
}

static void getActiveMat(VecDoub_I &x, MatDoub_IO &A, cnstrPhi cnstrnt[], int indx[], 
    const int mm01, const int mm) 
{
    int n=A.nrows(), midx=0;  // midx indexes the active constraints
    VecDoub cg(n);
    CnstrFuncd* cfuncd;
    for (int i=0; i<mm; i++) {
        if ( (i >= mm01) && (indx[i] != 1) ) { continue; }
        else {
            cfuncd = new CnstrFuncd(cnstrnt[i]);
            cfuncd->df(x,cg);
            for (int j=0; j<n; j++) { A[j][midx]=cg[j]; }
            midx++;
            delete cfuncd;
        }
    }
}



/*--------------------------------------------------------/
// 2. Projected Gradient and Projected Newton methods
/--------------------------------------------------------*/

static void CProj(MatDoub_I &A, MatDoub_O &S) 
{
    int n=A.nrows();
    int ma=A.ncols();
    QRdcmp qrA(A);
    qrA.factorise();
    for (int i=ma, k=0; i<n; i++, k++) {
        for (int j=0; j<n; j++) {
            S[j][k]=qrA.qt[i][j];
        }
    }
}

static void NewtonStep(MatDoub_I &A, VecDoub_I &c, VecDoub_O &y)
{
    int n=A.nrows(), ma=A.ncols();
    for (int i=0; i<n; i++) {
        y[i]=0.0;
        for (int j=0; j<ma; j++) {
            y[i] += A[i][j]*c[j];
        }
    }
    QRdcmp qrA(A);
    qrA.factorise();
    VecDoub z(n);
    qrA.qtmult(y,z);
    qrA.rsolve(z,y);
    qrA.rtsolve(y,z);
    qrA.qmult(z,y);
}

/**
 * Projected Gradient Descent
 *
 * @param x is an initial seeded value, and contains the final value upon
 * completion of the algorithm.  @param func is the objection function. 
 * @param cnstrnt[] is an array of constraint functions, the first 
 * @param mm01 of which contain equality contraints, and the latter 
 * @param mm02 of which contain inequality constraints.  @param stpmax 
 * specifies the maximum line search length. @param iter stores the 
 * iteration count.  @param fp stores the function value across iterations. 
 * @return void
 *
 * @exception singular reduced QR factorisation - does not occur by the
 * rank deficient nature of R, but by a zero in the diagonal of R.
 */
void ProjGradDesc(VecDoub_IO &x, objPhi &func, cnstrPhi cnstrnt[], 
    const int mm01, const int mm02, 
    const double stpmax, int &iter, double &fp)
{
    int n=x.size();  // total dimension of space
    // mm is total number of constraints, ma stores number of active constraints;
    int mm=mm01+mm02, ma;  // ma is at least mm01
    int* idx = new int[mm];  // records indices of active constraints
	MatDoub A(0,0);  // column-oriented active constraint gradients
    MatDoub S(0,0);  // complement image of A
    QRdcmp* qrA;  // holds QR factorisation of A

    // d is direction, y is for projected Newton iterations;
    VecDoub xnew(n), y(n), d(n), xi(n);  // xi is step difference xnew-xold
    VecDoub g(n), gr(0);  // g is gradient, g is reduced gradient
    double gdnorm;  // stores the norm of descent gradient
    VecDoub cn(mm), ca(0);  // cn stores values of constraints, ca of active constraints
    ObjFuncd funcd(func);  // objective function
    double alpha;  // line search length
    double fnew;  // new value of function
    int Niter, NNiter;  // Newton iterations
    double yerr;  // for calculating yerr in Newton iterations
    double gerr=DBL_MAX, xerr=DBL_MAX;  // for testing tolerance of errors
    double temp, den;  // for calculating xerr and gerr

    iter=-1;  // first descent iteration will start at 0
    while ( ((gerr > GTOL) && (xerr > TOLX)) && (iter++ < ITMAX) ) {
        evalConstrnts(x, cn, cnstrnt, mm);
        fp=funcd(x);
        getActiveIndx(cn, idx, mm01, mm, ma);
        funcd.df(x,g);
        if (ma==0) {
            gdnorm = MAX(VECNORM(g,2),1.0);
            for (int i=0; i<n; i++) { d[i]=-g[i]/gdnorm; }
        } else {
            A.resize(n,ma);
            S.resize(n,n-ma);
            getActiveMat(x, A, cnstrnt, idx, mm01, mm);
            CProj(A,S);
            gr.resize(n-ma);
            for (int i=0; i<n-ma; i++) {
                gr[i]=0.0;
                for (int j=0; j<n; j++) {
                    gr[i] += S[j][i]*g[j];
                }
            }
            gdnorm = MAX(VECNORM(gr,2),1.0);
            for (int i=0; i<n; i++) {
                d[i]=0.0;
                for (int j=0; j<n-ma; j++) {
                    d[i] -= S[i][j]*gr[j]/gdnorm;
                }
            }
        }
        alpha=stpmax;
        fnew=DBL_MAX;
        Niter=0;
        while ( (fnew > fp) && (Niter++ < ITMAX) ) {
            for (int i=0; i<n; i++) {
                xnew[i] = x[i] + alpha*d[i];
            }
            fnew = funcd(xnew);
            if (fnew > fp) { 
                alpha = 0.9*alpha;  // perhaps call line search here
                continue; 
            } else {
                NNiter=0;
                evalConstrnts(xnew, cn, cnstrnt, mm);
                while (!isFeasible(cn, mm01, mm) && (NNiter++ < NMAX)) {
                    getActiveIndx(cn, idx, mm01, mm, ma);
                    ca.resize(ma);
                    getActiveCn(cn, ca, idx, mm);
                    A.resize(n,ma);
                    getActiveMat(xnew, A, cnstrnt, idx, mm01, mm);
                    NewtonStep(A, ca, y);
#ifdef _DEBUG_PRINT_
                    fprintf(stderr, "NNiter: %d :: xnew :: ", NNiter);
                    for (int i=0; i<n; i++) {
                        fprintf(stderr, "%.3f ", xnew[i]);
                    }
                    fprintf(stderr, ":: ca :: ");
                    for (int i=0; i<ma; i++) {
                        fprintf(stderr, "%.3f ", ca[i]);
                    }
                    fprintf(stderr, ":: y :: ");
                    for (int i=0; i<n; i++) {
                        fprintf(stderr, "%.3f ", y[i]);
                    }
                    fprintf(stderr, "\n");
#endif
                    // calculate max relative change in xnew;
                    temp=0.0, yerr=0.0;
                    for (int i=0; i<n; i++) {
                        temp = abs(y[i]) / MAX(abs(xnew[i]),1.0);
                        if (temp > yerr) { yerr=temp; }
                    }
                    for (int i=0; i<n; i++) {
                        xnew[i] -= y[i];
                    }
                    if (yerr < TOLY) { break; }
                    evalConstrnts(xnew, cn, cnstrnt, mm);
                }
                fnew = funcd(xnew);
                if (fnew > fp) { 
                    alpha = 0.9*alpha;  // perhaps call line search here
                    continue; 
                }
            }
        }
        // calculate max relative change in x;
        for (int i=0; i<n; i++) { xi[i] = xnew[i]-x[i]; }
        temp=0.0, xerr=0.0;
        for (int i=0; i<n; i++) {
            temp = abs(xi[i]) / MAX(abs(x[i]),1.0);
            if (temp > xerr) { xerr=temp; }
        }
        // calculate max relative change in function;
        den = MAX(abs(fp),1.0);
        gerr=0.0;
        for (int i=0; i<n; i++) {
            temp = abs(g[i]) * MAX(abs(xi[i]),1.0) / den;
            if (temp > gerr) { gerr=temp; }
        }
        // update x and function value fp;
        for (int i=0; i<n; i++) { x[i] = xnew[i]; }
        fp=fnew;
#ifdef _ERROR_ANALYSIS_
        res_errs[iter]=gerr;
        step_errs[iter]=xerr;
        fvals[iter]=fp;
#endif
    }

    // clear memory;
    delete[] idx;
    return;
}


/*--------------------------------------------------------/
// 3. Projected Newton Descent
/--------------------------------------------------------*/


#endif /* _PROJDESC_H_ */
