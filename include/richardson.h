#ifndef _RICHARDSON_H_
#define _RICHARDSON_H_

/*------------------------------------------------------/
// we implement a suite of classical Richardson methods 
//   using the Numerical Recipes interface;
/------------------------------------------------------*/

//#include <float.h>  // already included in nr3math.h
//#include <cmath>  // already included in nr3math.h

#include "nr3math.h"

//---------------------------------------------------------
// 2. classical methods
//---------------------------------------------------------

void Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax);

void JOR_Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax);

void GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err, 
    const Int norm, const Doub tol, const Int itmax);

void SOR_GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax);



//---------------------------------------------------------
// 3. function definitions of classical methods
//---------------------------------------------------------

void Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax)
{
    Int n=A.nrows(), m=A.ncols();
    VecDoub r(n), xold(m);

    *err = (Doub) DBL_MAX;
    *iter = 0;
    while ( (abs(*err) > tol) && (++*iter < itmax) ) {
        // stores old x;
        for (int j=0; j<m; j++) { xold[j] = x[j]; }
        // updates new x;
        for (int j=0; j<m; j++) {
            // updates entirely from xold;
            x[j] = b[j] / A[j][j];
            for (int k=0; k<j; k++) {
                x[j] += -(A[j][k] / A[j][j])*xold[k];
            }
            for (int k=j+1; k<m; k++) {
                x[j] += -(A[j][k] / A[j][j])*xold[k];
            }
        }

        // calculate residual r=b-Ax;
        for (int i=0; i<n; i++) {
            r[i] = b[i];
            for (int j=0; j<m; j++) {
                r[i] -= A[i][j]*x[j];
            }
        }
        // get new residual error;
        *err = VECNORM(r, norm);
    }
}

void JOR_Jacobi(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax)
{
    Int n=A.nrows(), m=A.ncols();
    VecDoub r(n), xold(m);

    *err = (Doub) DBL_MAX;
    *iter = 0;
    while ( (abs(*err) > tol) && (++*iter < itmax) ) {
        // stores old x;
        for (int j=0; j<m; j++) { xold[j] = x[j]; }
        // updates new x;
        for (int j=0; j<m; j++) {
            // updates entirely from xold;
            x[j] = x[j] + wparam * b[j] / A[j][j];
            for (int k=0; k<m; k++) {
                x[j] += -wparam * (A[j][k] / A[j][j])*xold[k];
            }
        }

        // calculate residual r=b-Ax;
        for (int i=0; i<n; i++) {
            r[i] = b[i];
            for (int j=0; j<m; j++) {
                r[i] -= A[i][j]*x[j];
            }
        }
        // get new residual error;
        *err = VECNORM(r, norm);
    }
}


void GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax) 
{
    Int n=A.nrows(), m=A.ncols();
    VecDoub r(n);

    *err = (Doub) DBL_MAX;
    *iter = 0;
    while ( (abs(*err) > tol) && (++*iter < itmax) ) {
        // update x;
        for (int j=0; j<m; j++) {
            x[j] = b[j] / A[j][j];
            for (int k=0; k<j; k++) {
                // updates using the newer values of x;
                x[j] += -(A[j][k] / A[j][j]) * x[k];
            }
            for (int k=j+1; k<n; k++) {
                // updates using the older values of x;
                x[j] += -(A[j][k] / A[j][j]) * x[k];
            }
        }

        // calculate residual r=b-Ax;
        for (int i=0; i<n; i++) {
            r[i] = b[i];
            for (int j=0; j<m; j++) {
                r[i] -= A[i][j]*x[j];
            }
        }
        // get new residual error;
        *err = VECNORM(r, norm);
    }
}

void SOR_GaussSeidel(MatDoub_I A, VecDoub_I &b, VecDoub_IO &x, Doub wparam,
    Int *iter, Doub *err, const Int norm, const Doub tol, const Int itmax) 
{
    Int n=A.nrows(), m=A.ncols();
    VecDoub r(n);

    *err = (Doub) DBL_MAX;
    *iter = 0;
    while ( (abs(*err) > tol) && (++*iter < itmax) ) {
        // update x;
        for (int j=0; j<m; j++) {
            x[j] += x[j] - wparam*x[j] + wparam*b[j]/A[j][j];
            for (int k=0; k<j; k++) {
                // updates using the newer values of x;
                x[j] += -wparam * (A[j][k] / A[j][j]) * x[k];
            }
            for (int k=j; k<n; k++) {
                // updates using the older values of x;
                x[j] += -wparam * (A[j][k] / A[j][j]) * x[k];
            }
        }

        // calculate residual r=b-Ax;
        for (int i=0; i<n; i++) {
            r[i] = b[i];
            for (int j=0; j<m; j++) {
                r[i] -= A[i][j]*x[j];
            }
        }
        // get new residual error;
        *err = VECNORM(r, norm);
    }
}
    



#endif /* _RICHARDSON_H_ */