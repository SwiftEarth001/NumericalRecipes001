#ifndef _LINBCG_H_
#define _LINBCG_H_

/*--------------------------------------------------------/
// we expand this to a class which can handle -
//   steepest gradient descent GD, 
//   conjugate gradient descent CGD,
//   preconditioned conjugate gradient descent PCGD;
// although the caveat is no complex numbers;
/--------------------------------------------------------*/

//#include <float.h>  // already included in nr3math.h
#include "nr3math.h"
#include "cholesky.h"

#ifdef _ERROR_ANALYSIS_RESIDUALS_
extern double* res_errs;
#endif


class Linbcg {
public:
	MatDoub A;  // objective function matrix;
	Int nrow, ncol;

	MatDoub* P;  // preconditioner matrix;
	Cholesky* PP;

	Linbcg();
	Linbcg(MatDoub_IO &a);
	Linbcg(MatDoub_IO &a, MatDoub* P);
    Linbcg(int n, int m, double* a);
	Linbcg(int n, int m, double* a, MatDoub* P);
	~Linbcg();

	void assignPreconditioner(MatDoub* P);

	void GDsolve(VecDoub_I &b, VecDoub_IO &x, int *iter, double *reserr, double *steperr,
		const int norm, const double restol, const double steptol, const int itmax);
	void CGDsolve(VecDoub_I &b, VecDoub_IO &x, int *iter, double *reserr, double *steperr,
		const int norm, const double restol, const double steptol, const int itmax);
	void GaussSouthwellsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, 
		Doub *err, const Int norm, const Doub tol, const Int itmax);

/* 	void PCGDsolve(VecDoub_I &b, VecDoub_IO &x, const Int itol, 
		const Doub tol, const Int itmax, Int &iter, Doub &err); */

private:
	virtual void asolve(VecDoub_I &b, VecDoub_O &x, const Int itrnsp);
	virtual void atimes(VecDoub_I &x, VecDoub_O &r, const Int itrnsp);
	virtual void psolve(VecDoub_I &r, VecDoub_O &z);
};


/*------------------------------------------------------/
// 2.1 constructor and destructor routines
/------------------------------------------------------*/

Linbcg::Linbcg() {}

Linbcg::Linbcg(MatDoub_IO &a) : 
	A(a), nrow(A.nrows()), ncol(A.ncols()), 
	P(nullptr), PP(nullptr) {}

Linbcg::Linbcg(MatDoub_IO &a, MatDoub* P) : 
	A(a), nrow(A.nrows()), ncol(A.ncols()), P(P) 
{
	PP = new Cholesky((MatDoub_I) (*P));
}

Linbcg::Linbcg(int n, int m, double* a) : 
	nrow(n), ncol(m), A(n, m, a), 
	P(nullptr), PP(nullptr) {}

Linbcg::Linbcg(int n, int m, double* a, MatDoub* P) : 
	nrow(n), ncol(m), A(n, m, a), P(P) 
{
	PP = new Cholesky((MatDoub_I) (*P));
}

Linbcg::~Linbcg() 
{
	// deallocating the preconditioner is presumed to be
	//   handled by the routine calling Linbcg;
}


/*------------------------------------------------------/
// 2.2 ancillary routines
/------------------------------------------------------*/

void Linbcg::atimes(VecDoub_I &x, VecDoub_O &r, const Int itrnsp)
{
	for (int i=0; i<nrow; i++) {
		r[i] = 0;
		for (int j=0; j<ncol; j++) {
			r[i] += A[i][j]*x[j]; 
		}
	}
}

void Linbcg::asolve(VecDoub_I &b, VecDoub_O &x, const Int itrnsp) 
{

}

void Linbcg::psolve(VecDoub_I &r, VecDoub_O &z) 
{

}

void Linbcg::assignPreconditioner(MatDoub* PC)
{
	P = PC;
	if (PP != nullptr) {
		delete PP;
	}
	PP = new Cholesky((MatDoub_I) (*P));
}


/*------------------------------------------------------/
// 2.3 gradient descent routines
/------------------------------------------------------*/

void Linbcg::GDsolve(VecDoub_I &b, VecDoub_IO &x, int *iter, double *reserr, double *steperr,
	const int norm, const double restol, const double steptol, const int itmax)
{
	Int n=b.size();
	Doub alphk;
	VecDoub r(n), ar(n), dx(n);
	
	const double EPS=1.0e-07;
	*steperr=DBL_MAX;  // for storing step size;
	*reserr=DBL_MAX;  // initial residual error;
	*iter = 0;
	while ( ((*reserr > restol) && (*steperr > steptol)) && (++*iter < itmax) ) {
		// calculate residual r=b-Ax;
		atimes(x, r, 0);  // this is momentarily r=Ax
		for (int k=0; k<n; k++) {
			r[k] = b[k] - r[k];  // this will give r=b-Ax
		}

		// update alpha;
		alphk=0;
		if (P == nullptr) {
			// without preconditioner P;
			for (int k=0; k<n; k++) {
				alphk += r[k]*r[k];
			}
			atimes(r, ar, 0);
			Doub denom = 0;
			for (int k=0; k<n; k++) {
				denom += r[k]*ar[k];
			}
			alphk = alphk / denom;
			// update relative change in x;
			*steperr = alphk*VECNORM(r,norm) / MAX(VECNORM(x,norm),1);
			// update x;
			for (int k=0; k<n; k++) {
				x[k] = x[k] + alphk*r[k];
			}
		} else {
			// with preconditioner P;
			VecDoub z(n), az(n);
			PP->elsolve(r, z);
			for (int k=0; k<n; k++) {
				alphk += z[k]*r[k];
			}
			atimes(z, az, 0);
			Doub denom = 0;
			for (int k=0; k<n; k++) {
				denom += z[k]*az[k];
			}
			alphk = alphk / denom;
			// update relative change in x;
			*steperr = alphk*VECNORM(z,norm) / MAX(VECNORM(x,norm),1);
			// update x;
			for (int k=0; k<n; k++) {
				x[k] = x[k] + alphk*z[k];
			}
		}

		// update error; <-- use relative error;
		*reserr = VECNORM(r,norm) / MAX(VECNORM(b,norm),1);
#ifdef _ERROR_ANALYSIS_RESIDUALS_
		res_errs[(*iter)-1] = *reserr;
#endif
	}
}

void Linbcg::CGDsolve(VecDoub_I &b, VecDoub_IO &x, int *iter, double *reserr, double *steperr,
	const int norm, const double restol, const double steptol, const int itmax)
{
	Int n=b.size();
	Doub alphk=0, betk=0;
	VecDoub r(n), ar(n), p(n), ap(n);

	// initialize r and p;
	atimes(x, r, 0);
	for (int k=0; k<n; k++) { r[k] = b[k] - r[k]; }
	for (int k=0; k<n; k++) { p[k] = r[k]; }
	
	*steperr=DBL_MAX;  // for storing step size;
	*reserr=DBL_MAX;  // initial residual error;
	*iter=0;
	while ( ( (*reserr > restol) && (*steperr > steptol) ) && (++*iter < itmax) ) {
		// update alpha;
		alphk = 0;
		for (int k=0; k<n; k++) {
			alphk += p[k]*r[k];
		}
		atimes(p, ap, 0);
		Doub denom = 0;
		for (int k=0; k<n; k++) {
			denom += p[k]*ap[k];
		}
		alphk = alphk / denom;

		// update relative change in x;
		*steperr = alphk*VECNORM(p,2) / MAX(VECNORM(x,2),1);

		// update x;
		for (int k=0; k<n; k++) {
			x[k] = x[k] + alphk*p[k];
		}
		
		// update residual r=b-Ax;
		atimes(p, ap, 0);
		for (int k=0; k<n; k++) {
			r[k] = r[k] - alphk*ap[k];
		}

		// update beta;
		betk = 0;
		if (P == nullptr) {
			// without preconditioner P;
			for (int k=0; k<n; k++) {
				betk += ap[k]*r[k];
			}
			betk = betk / denom;  // we are not using complex values...
			// update direction p;
			for (int k=0; k<n; k++) {
				p[k] = r[k] - betk*p[k];
			}
		} else {
			// with preconditioner P;
			VecDoub z(n);
			PP->elsolve(r, z);
			for (int k=0; k<n; k++) {
				betk += ap[k]*z[k];
			}
			betk = betk / denom;  // we are not using complex values...
			// update direction p;
			for (int k=0; k<n; k++) {
				p[k] = z[k] - betk*p[k];
			}
		}

		// update error;
		*reserr = VECNORM(r, norm);
#ifdef _ERROR_ANALYSIS_RESIDUALS
		res_errs[(*iter)-1] = *reserr;
#endif
	}
}


/*------------------------------------------------------/
// 2.3 alternative descent routines
/------------------------------------------------------*/

void Linbcg::GaussSouthwellsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, 
	Doub *err, const Int norm, const Doub tol, const Int itmax) 
{
	Int n=b.size();
	Doub alphk;
	VecDoub r(n), d(n), ad(n);

	*err = (Doub) DBL_MAX;
	*iter = 0;
	while ( (abs(*err) > tol) && (++*iter < itmax) ) {
		// calculate residual r=b-Ax;
		atimes(x, r, 0);
		for (int k=0; k<n; k++) {
			r[k] = b[k] - r[k];
		}

		// update direction;
		Doub maxrk = 0.0;
		Int argmaxrk;
		for (int k=0; k<n; k++) {
			if (maxrk <= abs(r[k])) {
				argmaxrk = k;
				maxrk = r[k];
			}
		}
		for (int k=0; k<argmaxrk; k++) {
			d[k] = 0;
		}
		d[argmaxrk] = r[argmaxrk];
		for (int k=argmaxrk+1; k<n; k++) {
			d[k] = 0;
		}

		// update alpha;
		alphk=maxrk*maxrk;
		atimes(d, ad, 0);  // can use a quicker method for this...
		alphk = alphk / (d[argmaxrk]*ad[argmaxrk]);
		
		// update x;
		for (int k=0; k<n; k++) {
			x[k] = x[k] + alphk*d[k];
		}

		// update error;
		*err = VECNORM(r, norm);
	}
}


/* void Linbcg::PCGDsolve(VecDoub_I &b, VecDoub_IO &x, const Int itol, const Doub tol,
	const Int itmax, Int &iter, Doub &err)
{
	Doub ak,akden,bk,bkden=1.0,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	const Doub EPS=1.0e-14;
	Int j,n=b.size();
	VecDoub p(n),pp(n),r(n),rr(n),z(n),zz(n);
	iter=0;
	atimes(x,r,0);
	for (j=0;j<n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	//atimes(r,rr,0);
	if (itol == 1) {
		bnrm=snrm(b,itol);
		asolve(r,z,0);
	}
	else if (itol == 2) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
	}
	else if (itol == 3 || itol == 4) {
		asolve(b,z,0);
		bnrm=snrm(z,itol);
		asolve(r,z,0);
		znrm=snrm(z,itol);
	} else throw("illegal itol in linbcg");
	while (iter < itmax) {
		++iter;
		asolve(rr,zz,1);
		for (bknum=0.0,j=0;j<n;j++) bknum += z[j]*rr[j];
		if (iter == 1) {
			for (j=0;j<n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		} else {
			bk=bknum/bkden;
			for (j=0;j<n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(p,z,0);
		for (akden=0.0,j=0;j<n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(pp,zz,1);
		for (j=0;j<n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(r,z,0);
		if (itol == 1)
			err=snrm(r,itol)/bnrm;
		else if (itol == 2)
			err=snrm(z,itol)/bnrm;
		else if (itol == 3 || itol == 4) {
			zm1nrm=znrm;
			znrm=snrm(z,itol);
			if (abs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=abs(ak)*snrm(p,itol);
				err=znrm/abs(zm1nrm-znrm)*dxnrm;
			} else {
				err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(x,itol);
			if (err <= 0.5*xnrm) err /= xnrm;
			else {
				err=znrm/bnrm;
				continue;
			}
		}
		if (err <= tol) break;
	}
} */


#endif /* _LINBCG_H_ */
