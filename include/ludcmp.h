#ifndef _LUDCMP_H_
#define _LUDCMP_H_

/*------------------------------------------------------/
// ...;
/------------------------------------------------------*/

//#include <cmath>  // already included in nr3.h
#include <float.h>
#include "nr3.h"

struct LUdcmp {
	Int n;  // dimension of system
	MatDoub lu;  // matrix A
	VecInt P;  // row permutations
	VecInt Q;  // column permutations
	Doub di;  // sign on row permutation determinant
	Doub dj;  // sign on column permutation determinant

	// ~~~ constructors and destructors ~~~ //
	LUdcmp(int nrow, int ncol, double* a);
	LUdcmp(MatDoub_I &a);
	~LUdcmp();

	// ~~~ factorisation methods ~~~ //
	void Crout_pivot_factorise();  // uses implicit pivoting with diagonal scaling
	void partial_pivot_factorise();  // without diagonal scaling
	void rook_pivot_factorise();
	void complete_pivot_factorise();

	// ~~~ triangular matrix solve ~~~ //
	void solve(VecDoub_I &b, VecDoub_O &x);
	void solve(MatDoub_I &b, MatDoub_O &x);

	// ~~~ inverse and determinant ~~~ //	
	void inverse(MatDoub_O &ainv);
	Doub det();

	// ~~~ multiply P L U Q ~~~ //
	void multLUerror();

	void mprove(VecDoub_I &b, VecDoub_IO &x);
	
	MatDoub_I &aref;
	double* apref;
};


/*------------------------------------------------------/
// 2.1 constructor and destructor routines
/------------------------------------------------------*/

LUdcmp::LUdcmp(int nrow, int ncol, double* a) : 
	lu(nrow, ncol, a), apref(a), aref(lu), 
	n(nrow), P(n), Q(n)
{
	di=1; dj=1;
	for (int j=0; j<n; j++) {
		P[j] = j;
		Q[j] = j;
	}
}

LUdcmp::LUdcmp(MatDoub_I &a) : n(a.nrows()), lu(a), aref(a), P(n) 
{
	di=1; dj=1;
	for (int j=0; j<n; j++) {
		P[j] = j;
		Q[j] = j;
	}
}

LUdcmp::~LUdcmp() {}


/*------------------------------------------------------/
// 2.2 factorisation routines
/------------------------------------------------------*/

void LUdcmp::partial_pivot_factorise() 
{
	Int i, imax, j, k;
	Doub big, temp;
	di=1.0;

	// perform row reduction with row swapping;
	for (k=0; k<n; k++) {
		big=0.0;
		imax=k;
		// find row with (argmax, max) in column k;
		for (i=k; i<n; i++) {
			temp = abs(lu[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		if (big == 0.0) { throw("singular matrix in LUdcmp::partial_pivot_factorise"); }
		// swap row k with row imax;
		if (k != imax) {
			// the entire row is swapped;
			// in commuting the permutation past an L matrix,
			//   the L matrix rows are also swapped;
			for (j=0; j<n; j++) {
				SWAP(lu[imax][j], lu[k][j]);  // SWAP is inline
			}
			SWAP(P[k],P[imax]);  // set row permutation
			di = -di;  // adjust sign on row permutation determinant
		}

		if (lu[k][k] == 0.0) { throw("singular matrix in LUdcmp::partial_pivot_factorise"); }
		// update the values for the lower portion both
		//   left and right of the diagonal;
		for (i=k+1; i<n; i++) {
			temp = (lu[i][k] /= lu[k][k]);  // this also recomputes lu[i][k]
			for (j=k+1; j<n; j++)
				lu[i][j] -= temp*lu[k][j];
		}
	}
}

void LUdcmp::Crout_pivot_factorise() 
{
	const Doub TINY=1.0e-40;
	Int i, imax, j, k;
	Doub big, temp;
	di=1.0;
	// vv acts as a diagonal scaling matrix
	//   acting from the left;
	VecDoub vv(n);

	// find the largest element for each 
	//   row i and store its reciprocal 
	//   in vv[i];
	for (i=0; i<n; i++) {
		big=0.0;
		for (j=0; j<n; j++) {
			if ((temp=abs(lu[i][j])) > big) { big = temp; }
		}
		if (big == 0.0) { throw("singular matrix in LUdcmp::Crout_pivot_factorise"); }
		vv[i] = 1.0/big;
	}

	// perform row reduction with row swapping;
	for (k=0; k<n; k++) {
		big=0.0;
		imax=k;
		// get largest row after multiplying
		//   vv from the left;
		for (i=k; i<n; i++) {
			temp = vv[i]*abs(lu[i][k]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		// swap row k with row imax;
		if (k != imax) {
			for (j=0; j<n; j++) {
				temp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=temp;
			}
			di = -di;  // adjust sign on row permutation determinant
			SWAP(P[k],P[imax]);
			vv[imax]=vv[k];
		}
		if (lu[k][k] == 0.0) { lu[k][k]=TINY; }
		// update the values for the lower portion both
		//   left and right of the diagonal;
		for (i=k+1; i<n; i++) {
			temp = (lu[i][k] /= lu[k][k]);  // this also recomputes lu[i][k]
			for (j=k+1; j<n; j++)
				lu[i][j] -= temp*lu[k][j];
		}
	}
}

void LUdcmp::rook_pivot_factorise() 
{
	Int i, imax, jmax, j, k;
	Doub big, max, temp;
	di=1.0, dj=1.0;

	// perform row reduction with row swapping;
	for (k=0; k<n; k++) {
		// find the largest element for row k;
		big=0.0;
		jmax=k;
		for (j=k; j<n; j++) {
			if ((temp=abs(lu[k][j])) > big) { 
				big=temp; 
				jmax=j;
			}
		}
		if (big == 0.0) { throw("singular matrix in LUdcmp::rook_pivot_factorise"); }
		// perform the column swap implicitly;
		SWAP(Q[k],Q[jmax]);
		dj = -dj;

		big=0.0;
		imax=k;
		// get largest row after the column swap;
		for (i=k; i<n; i++) {
			temp = abs(lu[i][Q[k]]);
			if (temp > big) {
				big=temp;
				imax=i;
			}
		}
		// perform row swap explicitly;
		// this operation is stable and should be 
		//   invariant to the column swap;
		if (k != imax) {
			for (j=0; j<n; j++) {
				temp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=temp;
			}
			SWAP(P[k],P[imax]);
			di = -di;
		}
		if (lu[k][Q[k]] == 0.0) { throw("singular matrix in LUdcmp::rook_pivot_factorise"); }
		// update the values for the lower portion both
		//   left and right of the diagonal;
		for (i=k+1; i<n; i++) {
			temp = (lu[i][Q[k]] /= lu[k][Q[k]]);  // this also recomputes lu[i][Q[k]]
			for (j=k+1; j<n; j++)
				lu[i][Q[j]] -= temp*lu[k][Q[j]];
		}
	}
}

void LUdcmp::complete_pivot_factorise() 
{
	Int i, imax, jmax, j, k;
	Doub big, max, temp;
	di=1.0, dj=1.0;

	// perform row reduction with row swapping;
	//   there is no need for diagonal scaling
	//   in this case;
	for (k=0; k<n; k++) {
		max=0.0;
		imax=k;
		jmax=k;
		// find the largest element in the submatrix
		//   and record its index;
		for (i=k; i<n; i++) {
			for (j=k; j<n; j++) {
				if ((temp=abs(lu[i][j])) > max) { 
					max = temp;
					imax = i;
					jmax = j;
				}
			}
		}
		// perform the column swap implicitly;
		SWAP(Q[k],Q[jmax]);
		dj = -dj;
		// perform the row swap explicitly;
		// this operations is stable and should be 
		//   invariant to the column swap;
		if (k != imax) {
			for (j=0; j<n; j++) {
				temp=lu[imax][j];
				lu[imax][j]=lu[k][j];
				lu[k][j]=temp;
			}
			SWAP(P[k],P[imax]);
			di = -di;
		}
		
		if (lu[k][Q[k]] == 0.0) { throw("singular matrix in LUdcmp::complete_pivot_factorise"); }
		// update the values for the lower portion both
		//   left and right of the diagonal;
		for (i=k+1; i<n; i++) {
			temp = (lu[i][Q[k]] /= lu[k][Q[k]]);  // this also recomputes lu[k][Q[k]]
			for (j=k+1; j<n; j++)
				lu[i][Q[j]] -= temp*lu[k][Q[j]];
		}
	}
}


/*------------------------------------------------------/
// 2.3 triangular matrix solve routines
/------------------------------------------------------*/

void LUdcmp::solve(VecDoub_I &b, VecDoub_O &x)
{
	Int i,ii=0,ip,j;
	Doub sum;
	if (b.size() != n || x.size() != n)
		throw("LUdcmp::solve bad sizes");
	for (i=0; i<n; i++) { x[i] = b[i]; }
	for (i=0; i<n; i++) {
		ip=P[i];
		sum=x[ip];
		x[ip]=x[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
		else if (sum != 0.0)
			ii=i+1;
		x[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=x[i];
		for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
		x[i]=sum/lu[i][i];
	}
}

void LUdcmp::solve(MatDoub_I &b, MatDoub_O &x)
{
	int i,j,m=b.ncols();
	if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
		throw("LUdcmp::solve bad sizes");
	VecDoub xx(n);
	for (j=0;j<m;j++) {
		for (i=0;i<n;i++) xx[i] = b[i][j];
		solve(xx,xx);
		for (i=0;i<n;i++) x[i][j] = xx[i];
	}
}


/*------------------------------------------------------/
// 2.4 matrix inverse and determinant routines
/------------------------------------------------------*/

void LUdcmp::inverse(MatDoub_O &ainv)
{
	Int i,j;
	ainv.resize(n,n);
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) ainv[i][j] = 0.;
		ainv[i][i] = 1.;
	}
	solve(ainv,ainv);
}
Doub LUdcmp::det()
{
	Doub dd = di;
	for (Int i=0;i<n;i++) dd *= lu[i][i];
	return dd;
}


void LUdcmp::mprove(VecDoub_I &b, VecDoub_IO &x)
{
	Int i,j;
	VecDoub r(n);
	for (i=0;i<n;i++) {
		Ldoub sdp = -b[i];
		for (j=0;j<n;j++)
			sdp += (Ldoub)aref[i][j] * (Ldoub)x[j];
		r[i]=sdp;
	}
	solve(r,r);
	for (i=0;i<n;i++) x[i] -= r[i];
}


/*------------------------------------------------------/
// 2.5 calculate relative error
/------------------------------------------------------*/

void LUdcmp::multLUerror() {
	double* A = new double[(n)*(n)];
	delete[] A;
}

#endif /* _LUDCMP_H_ */
