#ifndef _QUASINEWTON_H_
#define _QUASINEWTON_H_

/*------------------------------------------------------/
// we implement a suite of quasi-Newton methods 
//   using the Numerical Recipes interface;
/------------------------------------------------------*/

//#include <float.h>  // already included in nr3math.h
//#include <cmath>  // already included in nr3.h
#include "nr3math.h"
#include "linsearch.h"
#include "linbcg.h"


/*--------------------------------------------------------/
// 1. static (and global) scope parameters
/--------------------------------------------------------*/

// starting parameters;
static int ITMAX=250;  // default maximum iterations
static const double EPS32=(double)(numeric_limits<float>::epsilon());  // EPS32 is float32 precision
static double TOLX=1.0e-13, GTOL=1.0e-13;  // tolerances for descent-step and residual errors
static const double TOLZERO=1.0e-13;  // tolerance for zero value
static double DATOL=1.0e-13;  // minimum tolerance for square-cosine of descent angle
static double STPMX=1.0;  // maximum step length alpha

#ifdef _ERROR_ANALYSIS_
extern double* step_errs;
extern double* fvals;
#endif


/**
 * -------------------------------------------------------------------------/
 * NewtonMethods - General Documentation
 *
 * @param p is an initial seeded value which contains the final iterate upon
 * completion of the algorithm.  @param func is the objection function. 
 * @param maxiter is the maximum descent iterations.  @param iter stores the
 * descent iteration count.  @param fp stores the function value for the
 * current iterate.  @param lsparams specifies the inexact line search 
 * parameters.
 *--------------------------------------------------------------------------/
 */


/*--------------------------------------------------------/
// 2. modified Newton methods
/--------------------------------------------------------*/

/**
 * TruncatedCGDupdate
 * Uses conjugate gradient descent CGD with an approximated Hessian
 * at each iteration.
 * @see NewtonMethods;
 * @exceptsafe No exceptions are thrown.
 */
void TruncatedCGDupdate(
	VecDoub_IO &p, objPhi &func, 
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams, MatDoub* P)
{
	int n=p.size();
	VecDoub g(n), pnew(n), d(n), xi(n);
	MatDoub hessian(n,n);

	// estimate the gradient through finite differences;
	ObjFuncd funcd(func);
	fp=funcd(p);
	funcd.df(p,g);
	// estimate the Hessian through finite differences;
	funcd.Hess(p,hessian);
	// set up conjugate gradient descent object;
	Linbcg TruncCGD(hessian);
	if (P != nullptr) { TruncCGD.assignPreconditioner(P); }

	double fnew;  // stores new function values;
	double temp, den, gerr=DBL_MAX, xerr=DBL_MAX;  // for tolerance tests
	for (iter=0; ( ((xerr > TOLX) && (gerr > GTOL)) && (iter<maxiter) ); iter++) {
		// update direction dk;
		int cgsiter=0; double cgsreserr=DBL_MAX; double cgssteperr=DBL_MAX;
		TruncCGD.CGDsolve(g, d, &cgsiter, &cgsreserr, &cgssteperr, 2, 1.0e-11, 1.0e-11, n+50);
		// normalise and reverse the direction dk so that it solves hessian*dk = -gk;
		for (int i=0; i<n; i++) { d[i]=-d[i]/MAX(VECNORM(d,2),1); }
		// line search along direction dk;
		inexact_linsearch(p, d, func, lsparams);
		if (!lsparams->check) { lsparams->alpha=1.0; }
		for (int i=0; i<n; i++) {
			pnew[i]=p[i]+lsparams->alpha*d[i];
			xi[i]=pnew[i]-p[i];
		}
#ifdef _DEBUG_PRINT_
		fprintf(stderr, "iter: %d :: alpha: %.3f :: check: %s :: line search iterations: %d\n", 
			iter, lsparams->alpha, lsparams->check ? "true" : "false", lsparams->iters);
#endif
		fnew=funcd(pnew);
		if (fnew >= fp) {
			return;
		} else {
			xerr=0.0;
			for (int i=0; i<n; i++) {
				temp = abs(xi[i]) / MAX(abs(p[i]),1.0);
				if (temp > xerr) { xerr=temp; }
			}
			// calculate relative change in function;
			gerr=0.0;
			den = MAX(abs(fp),1.0);
			for (int i=0; i<n; i++) {
				temp = abs(g[i]) * MAX(abs(d[i]),1.0) / den;
				if (temp > gerr) { gerr=temp; }
			}
			// update p to pnew;
			for (int i=0; i<n; i++) { p[i]=pnew[i]; }
			funcd.df(p,g);  // new gradient
			fp=fnew;  // new function value
			funcd.Hess(p,hessian);  // update Hessian
#ifdef _ERROR_ANALYSIS_
			step_errs[iter]=VECNORM(xi,2);
			fvals[iter]=fp;
#endif
		}
	}
	return;
//	throw("too many iterations in CGSupdate");
}



static void GMSCholesky(MatDoub_I &G, MatDoub_IO &L, VecDoub_IO &D, VecDoub_IO &E, double delta)
{
	int n=D.size();

	MatDoub C(n,n,0.0);
	for (int i=0; i<n; i++) { C[i][i]=G[i][i]; }
	for (int i=0; i<n; i++) { L[i][i]=1.0; }

	double gamma=0, xi=0;
	for (int i=0; i<n; i++) { 
		gamma = abs(G[i][i]) > gamma ? abs(G[i][i]) : gamma;
		for (int j=i+1; j<n; j++) {
			xi = abs(G[i][j]) > xi ? abs(G[i][j]) : xi;
		}
	}
	const double EPS=numeric_limits<double>::epsilon();  // EPS is machine precision
	double beta2 = std::max( { gamma, xi/sqrt(n*n-1), EPS } );

	double sum, theta, cmax;
	for (int i=0; i<n; i++) {
		// <-- ... there is "supposed" to be a cmax swap here ... --->
		// <-- ... likely a stability concern, so maybe not necessary ... --->
		for (int j=0; j<i; j++) { L[i][j]=C[i][j]/D[j]; }
		sum=0.0;
		for (int j=i+1; j<n; j++) {
			for (int s=0; s<i; s++) {
				sum += L[i][s]*C[j][s];
			}
			C[j][i]=G[j][i]-sum;
		}
		theta=0.0;
		for (int j=i+1; j<n; j++) {
			// theta remains 0 when i=n;
			theta = abs(C[j][i]) > theta ? abs(C[j][i]) : theta;
		}
		D[i] = std::max( { delta, abs(C[i][i]), theta*theta/(beta2) } );
		E[i] = D[i]-C[i][i];
		for (int j=i+1; j<n; j++) {
			// update to C does not occur when i=n;
			C[j][j] = C[j][j] - C[j][i]*C[j][i]/D[i];
		}
	}
}

/**
 * GMSupdate
 * Performs Gill-Murray stable update.
 * @see NewtonMethods;
 * @exceptsafe No exceptions are thrown.
 */
void GMSupdate(
	VecDoub_IO &p, objPhi &func,
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams) 
{ 
	int n=p.size();
	VecDoub g(n), pnew(n), d(n), xi(n), D(n), E(n);
	MatDoub hessian(n,n), L(n,n);
	double delta=1.0;

	// estimate the gradient through finite differences;
	ObjFuncd funcd(func);
	fp=funcd(p);
	funcd.df(p,g);
	// set the initial line direction to opposite the gradient;
	for (int i=0; i<n; i++) { d[i] = -g[i]; }
	// estimate the Hessian through finite differences;
	funcd.Hess(p,hessian);

	double fnew;  // stores new function value;
	double gnorm=DBL_MAX, sum=0.0;  // gdnorm stores gradient norm;
	double temp, den, gerr=DBL_MAX, xerr=DBL_MAX;  // for tolerance tests
	for (iter=0; iter<maxiter; iter++) {
		// compute the Bill & Murray modified Cholesky factorisation;
		GMSCholesky(hessian, L, D, E, delta);

		gnorm=VECNORM(g,2);
		// update direction dk;
		if (gnorm > TOLZERO) {
			// solve LDL^t*dk = -gk;
			for (int i=0; i<n; i++) {
				sum=-g[i];
				for (int k=i-1; k>=0; k--) { sum -= L[i][k]*d[k]; }
				d[i]=sum/L[i][i];
			}
			for (int i=0; i<n; i++) { d[i] /= D[i]; }
			for (int i=n-1; i>=0; i--) {
				sum=d[i];
				for (int k=i+1; k<n; k++) { sum -= L[k][i]*d[k]; }
				d[i]=sum/L[i][i];
			}
		} else {
			// solve L^tdk = et, where t is index
			//   of minimum dj-ej;
			VecDoub et(n,0.0);
			double psi=DBL_MAX, val; int idx;
			for (int i=0; i<n; i++) {
				val = D[i]-E[i];
				if (psi < val) { 
					psi=val;
					idx=i;
				}
			}
			if (psi < 0) {
				et[idx]=1.0;
				for (int i=n-1; i>=0; i--) {
					sum=et[i];
					for (int k=i+1; k<n; k++) { sum -= L[k][i]*d[k]; }
					d[i]=sum/L[i][i];
				}
			} else {
				return;
			}
		}

		// orient dk in opposite cone of gradient;
		double dtg=0;
		for (int i=0; i<n; i++) { dtg += d[i]*g[i]; }
		if (dtg > 0) {
			for (int i=0; i<n; i++) { d[i]=-d[i]; }
		}

		// line search along direction dk;
		inexact_linsearch(p, d, func, lsparams);
		if (!lsparams->check) { lsparams->alpha=1.0; }
		for (int i=0; i<n; i++) {
			pnew[i]=p[i]+lsparams->alpha*d[i];
			xi[i]=pnew[i]-p[i];
		}
#ifdef _DEBUG_PRINT_
		fprintf(stderr, "iter: %d :: alpha: %.3f :: check: %s :: line search iterations: %d\n", 
			iter, lsparams->alpha, lsparams->check ? "true" : "false", lsparams->iters);
#endif
		fnew=funcd(pnew);
		if (fnew >= fp) {
			return;
		} else {
			xerr=0.0;
			for (int i=0; i<n; i++) {
				temp = abs(xi[i]) / MAX(abs(p[i]),1.0);
				if (temp > xerr) { xerr=temp; }
			}
			// calculate relative change in function;
			gerr=0.0;
			den = MAX(abs(fp),1.0);
			for (int i=0; i<n; i++) {
				temp = abs(g[i]) * MAX(abs(d[i]),1.0) / den;
				if (temp > gerr) { gerr=temp; }
			}
			// update p to pnew;
			for (int i=0; i<n; i++) { p[i]=pnew[i]; }
			funcd.df(p,g);  // new gradient
			fp=fnew;  // new function value
#ifdef _ERROR_ANALYSIS_
			step_errs[iter]=VECNORM(xi,2);
			fvals[iter]=fp;
#endif
			funcd.Hess(p,hessian);  // update Hessian
		}
	}
	return;
//	throw("too many iterations in GMSupdate");
}


/*--------------------------------------------------------/
// 2. rank-2 quasi Newton methods
/--------------------------------------------------------*/

enum QNtype { DFP=0, BFGS=2 };

/**
 * QNupdate - Quasi Newton update
 * Performs either BFGS or DFP update depending on qntype specified.
 * @see NewtonMethods; 
 * @param qntype specifies whether or not the update is DFP or BFGS. 
 * @relates BFGSupdate, DFPupdate;
 * @return void
 *
 * @exceptsafe No exceptions are thrown.
 */
static void QNupdate(
	VecDoub_IO &p, objPhi &func, QNtype qntype,
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams)
{
	int n=p.size();
	VecDoub dg(n), g(n), hdg(n), pnew(n), xi(n);
	MatDoub hessin(n,n);  // inverse Hessian

	// estimate the gradient through finite differences;
	ObjFuncd funcd(func);
	fp=funcd(p);
	funcd.df(p,g);
	
	for (int i=0; i<n; i++) {
		// set initial inverse Hessian to identity;
		for (int j=0; j<n; j++) { hessin[i][j]=0.0; }
		hessin[i][i]=1.0;
		// set the initial line direction to opposite the gradient;
		xi[i] = -g[i];
	}

	double temp, den, gerr=DBL_MAX, xerr=DBL_MAX;  // for tolerance tests
	for (iter=0; ( ((xerr > TOLX) && (gerr > GTOL)) && (iter<maxiter) ); iter++) {
		inexact_linsearch(p, xi, func, lsparams);
		if (!lsparams->check) { lsparams->alpha=1.0; }
		for (int i=0; i<n; i++) {
			pnew[i]=p[i]+lsparams->alpha*xi[i];
			xi[i]=pnew[i]-p[i];
		}
#ifdef _DEBUG_PRINT_
		fprintf(stderr, "iter: %d :: alpha: %.3f :: check: %s :: line search iters: %d\n", 
			iter, lsparams->alpha, lsparams->check ? "true" : "false", lsparams->iters);
#endif
		// calculate relative change in p;
		xerr=0.0;
		for (int i=0; i<n; i++) {
			temp = abs(xi[i]) / MAX(abs(p[i]),1.0);
			if (temp > xerr) { xerr=temp; }
		}
		// calculate relative change in function func;
		gerr=0.0;
		den = MAX(abs(fp),1.0);
		for (int i=0; i<n; i++) {
			temp = abs(g[i]) * MAX(abs(xi[i]),1.0) / den;
			if (temp > gerr) { gerr=temp; }
		}
		fp=funcd(pnew);  // new function value
#ifdef _ERROR_ANALYSIS_
		step_errs[iter]=VECNORM(xi,2);
		fvals[iter]=fp;
#endif
		// update p to pnew;
		for (int i=0; i<n; i++) { p[i]=pnew[i]; }
		// momentarily store the old gradient in dg;
		for (int i=0; i<n; i++) { dg[i]=g[i]; }
		funcd.df(p,g);  // new gradient
		fp=funcd(pnew);  // new function value
		// calculate gradient difference into dg;		
		for (int i=0; i<n; i++) { dg[i]=g[i]-dg[i]; }
		for (int i=0; i<n; i++) {
			hdg[i]=0.0;
			// calculate H_k*y_k into hdg;
			for (int j=0; j<n; j++) { hdg[i] += hessin[i][j]*dg[j]; }
		}

		double fac=0.0, fae=0.0, sumdg=0.0, sumxi=0.0;
		for (int i=0; i<n; i++) {
			// calculate s_k*y_k into fac;
			fac += dg[i]*xi[i];
			// calculate y_k*H_k*y_k into fae;
			fae += dg[i]*hdg[i];
			// sumdg is square-norm of dg,
			//  and sumxi is square-norm of xi;
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}

		if (fac > sqrt(DATOL*sumdg*sumxi)) {
			double faf=1.0/fac, fad=1.0/fae;
			switch(qntype) {
			case BFGS:
				for (int i=0; i<n; i++) {
					for (int j=i; j<n; j++) {
						// hdg and hdg-transpose have the same elements
						//   due to symmetry of H_k,
						//   i.e.  hdg[j] = hdg-transpose[j];
						hessin[i][j] += (faf*faf)*(fac+fae)*xi[i]*xi[j] - faf*(hdg[i]*xi[j]+xi[i]*hdg[j]);
						hessin[j][i]=hessin[i][j];  // inverse Hessian is symmetric
					}
				}
				break;
			case DFP:
				for (int i=0; i<n; i++) {
					for (int j=i; j<n; j++) {
						// hdg and hdg-transpose have the same elements
						//   due to symmetry of H_k,
						//   i.e.  hdg[j] = hdg-transpose[j];
						hessin[i][j] += (faf)*(xi[i]*xi[j]) - (fad)*(hdg[i]*hdg[j]);
						hessin[j][i]=hessin[i][j];  // inverse Hessian is symmetric
					}
				}
				break;
			default:
				break;
			}
		} else {
			// ... Liu and Fukishima cautious update ... Hessian is unchanged ;
		}
		for (int i=0; i<n; i++) {
			xi[i]=0.0;
			for (int j=0; j<n; j++) { 
				xi[i] -= hessin[i][j]*g[j];  // update xi to direction d_k+1
			}
		}
	}
	return;
//	throw("too many iterations in BFGSupdate");
}

/**
 * BFGSupdate
 * Entry point into the Broyden-Fletcher-Goldfarb-Shanno BFGS update.
 * @see This method calls QNupdate;
 * @see NewtonMethods;
*/
void BFGSupdate(
	VecDoub_IO &p, objPhi &func,
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams) 
{
	QNupdate(p, func, BFGS, maxiter, iter, fp, lsparams);
	return;
}

/**
 * DFPupdate
 * Entry point into the Davidon-Fletcher-Powell DFP update.
 * @see NewtonMethods;
*/
void DFPupdate(
	VecDoub_IO &p, objPhi &func,
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams) 
{ 
	QNupdate(p, func, DFP, maxiter, iter, fp, lsparams);
	return;
}

/*--------------------------------------------------------/
// 2. rank-1 quasi Newton methods
/--------------------------------------------------------*/

void PSBupdate(
	VecDoub_IO &p, objPhi &func, 
	int maxiter, int &iter, double &fp, 
	line_search_params* lsparams);


#endif /* _QUASINEWTON_H_ */
