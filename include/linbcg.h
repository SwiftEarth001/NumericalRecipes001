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


class Linbcg {
public:
	MatDoub A;  // objective function matrix;
	Int nrow, ncol;

	MatDoub P;  // preconditioner matrix;

	Linbcg();
	Linbcg(MatDoub_IO &a);
    Linbcg(int n, int m, double* a);
	~Linbcg();

	void GDsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
		const Int norm, const Doub tol, const Int itmax);
	void CGDsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
		const Int norm, const Doub tol, const Int itmax);
	void GaussSouthwellsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, 
		Doub *err, const Int norm, const Doub tol, const Int itmax);
	// void PCGDsolve(VecDoub_I &b, VecDoub_IO &x, const Int itol, 
	//	const Doub tol, const Int itmax, Int &iter, Doub &err);

private:
	// virtual void asolve(VecDoub_I &b, VecDoub_O &x, const Int itrnsp);
	virtual void atimes(VecDoub_I &x, VecDoub_O &r, const Int itrnsp);
};

Linbcg::Linbcg() {}

Linbcg::Linbcg(MatDoub_IO &a) : 
	A(a), nrow(A.nrows()), ncol(A.ncols()) {}

Linbcg::Linbcg(int n, int m, double* a) : 
	nrow(n), ncol(m), A(n, m, a) {}

Linbcg::~Linbcg() {}



void Linbcg::atimes(VecDoub_I &x, VecDoub_O &r, const Int itrnsp)
{
	for (int i=0; i<nrow; i++) {
		r[i] = 0;
		for (int j=0; j<ncol; j++) {
			r[i] += A[i][j]*x[j]; 
		}
	}
}

void Linbcg::GDsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax)
{
	Int n=b.size();
	Doub alphk;
	VecDoub r(n), ar(n);
	
	*err = (Doub) DBL_MAX;
	*iter = 0;
	while ( (abs(*err) > tol) && (++*iter < itmax) ) {
		// calculate residual r=b-Ax;
		atimes(x, r, 0);
		for (int k=0; k<n; k++) {
			r[k] = b[k] - r[k];
		}

		// update alpha;
		alphk=0;
		for (int k=0; k<n; k++) {
			alphk += r[k]*r[k];
		}
		atimes(r, ar, 0);
		Doub denom = 0;
		for (int k=0; k<n; k++) {
			denom += r[k]*ar[k];
		}
		alphk = alphk / denom;
		
		// update x;
		for (int k=0; k<n; k++) {
			x[k] = x[k] + alphk*r[k];
		}

		// update error;
		*err = VECNORM(r, norm);
	}
}

void Linbcg::CGDsolve(VecDoub_I &b, VecDoub_IO &x, Int *iter, Doub *err,
	const Int norm, const Doub tol, const Int itmax)
{
	Int n=b.size();
	Doub alphk=0, betk=0;
	VecDoub r(n), ar(n), p(n), ap(n);

	// initialize r and p;
	atimes(x, r, 0);
	for (int k=0; k<n; k++) { r[k] = b[k] - r[k]; }
	for (int k=0; k<n; k++) { p[k] = r[k]; }
	
	// get initial error;
	*err = VECNORM(r, norm);
	*iter = 0;
	while ( (abs(*err) > tol) && (++*iter < itmax) ) {
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
		for (int k=0; k<n; k++) {
			betk += ap[k]*r[k];
		}
		betk = betk / denom;  // we are not using complex values...

		// update direction p;
		for (int k=0; k<n; k++) {
			p[k] = r[k] - betk*p[k];
		}

		// update error;
		*err = VECNORM(r, norm);
	}
}

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
