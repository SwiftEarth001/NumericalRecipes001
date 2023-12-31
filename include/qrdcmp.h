#ifndef _QRDCMP_H_
#define _QRDCMP_H_

/*------------------------------------------------------/
// our essential modification to the code is that
//   the factorization scheme can apply to non-square
//   matrices; pivoting should not be necessary as
//   far as we can tell; we do, however, presume that
//   the count of rows should either exceed or equal
//   the count of columns;
/------------------------------------------------------*/

#include "nr3.h"

class QRdcmp 
{
private:
	Int n;
	Int m;
	Bool sing;
public:
	MatDoub qt, r;  // qt is Q-transpose, whereas A=QR
	
	QRdcmp(int nrow, int ncol, double* a);
	QRdcmp(MatDoub_I &a);
	~QRdcmp();
	
	void factorise();

	void solve(VecDoub_I &b, VecDoub_O &x);
	void qtmult(VecDoub_I &b, VecDoub_O &x);
	void qmult(VecDoub_I &b, VecDoub_O &x);
	void rsolve(VecDoub_I &b, VecDoub_O &x);
	void rtsolve(VecDoub_I &b, VecDoub_O &x);
	
	void update(VecDoub_I &u, VecDoub_I &v);
	void rotate(const Int i, const Doub a, const Doub b);
};

QRdcmp::QRdcmp(int nrow, int ncol, double* a) :
	n(nrow), m(ncol), sing(false),
	qt(n,n), r(n,m,(Doub*)(a)) { }

QRdcmp::QRdcmp(MatDoub_I &a) : 
	n(a.nrows()), m(a.ncols()), sing(false), 
	qt(n,n), r(a) { }

QRdcmp::~QRdcmp() { }

void QRdcmp::factorise()
{
	int i,j,k;
	VecDoub c(n), d(n);
	Doub scale,sigma,sum,tau;
	
	for (k=0;k<m;k++) {
		scale=0.0;
		for (i=k;i<n;i++) scale=MAX(scale,abs(r[i][k]));
		if (scale == 0.0) {
			sing=true;
			c[k]=d[k]=0.0;
		} else {
			for (i=k;i<n;i++) r[i][k] /= scale;
			for (sum=0.0,i=k;i<n;i++) sum += SQR(r[i][k]);
			sigma=SIGN(sqrt(sum),r[k][k]);
			c[k] = 2 * (sum + abs(r[k][k])*sqrt(sum));
			d[k] = -scale*sigma;
			r[k][k] += sigma;
			for (j=k+1;j<m;j++) {
				for (sum=0.0,i=k;i<n;i++) sum += r[i][k]*r[i][j];
				tau=sum/c[k];  // <-- beware, floating point division
				for (i=k;i<n;i++) r[i][j] -= 2*tau*r[i][k];
			}
		}
	}
	
	if (d[m-1] == 0.0) sing=true;
	for (i=0;i<n;i++) {
		for (j=0;j<n;j++) qt[i][j]=0.0;
		qt[i][i]=1.0;
	}
	
	for (k=0;k<m;k++) {
		if (c[k] != 0.0) {
			for (j=0;j<n;j++) {
				sum=0.0;
				for (i=k;i<n;i++)
					sum += r[i][k]*qt[i][j];
				sum /= c[k];
				for (i=k;i<n;i++)
					qt[i][j] -= 2*sum*r[i][k];
			}
		}
	}
	for (i=0;i<m;i++) {
		r[i][i]=d[i];
		for (j=0;j<i;j++) r[i][j]=0.0;
	}
	for (i=m;i<n;i++) {
		for (j=0;j<m;j++) r[i][j]=0.0;
	}
}

void QRdcmp::solve(VecDoub_I &b, VecDoub_O &x) {
	qtmult(b,x);
	rsolve(x,x);
}

void QRdcmp::qtmult(VecDoub_I &b, VecDoub_O &x) {
	int i,j;
	Doub sum;
	for (i=0;i<n;i++) {
		sum = 0.0;
		for (j=0;j<n;j++) sum += qt[i][j]*b[j];
		x[i] = sum;
	}
}

void QRdcmp::qmult(VecDoub_I &b, VecDoub_O &x) {
	int i,j;
	Doub sum;
	for (i=0;i<n;i++) {
		sum = 0.0;
		for (j=0;j<n;j++) sum += qt[j][i]*b[j];
		x[i] = sum;
	}
}

void QRdcmp::rsolve(VecDoub_I &b, VecDoub_O &x) {
	int i,j;
	Doub sum;
	if (sing) throw("attempting solve in a singular QR");
	for (i=m-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<m;j++) sum -= r[i][j]*x[j];
		x[i]=sum/r[i][i];
	}
	for (i=m;i<n;i++) { x[i]=0.0; }
}

void QRdcmp::rtsolve(VecDoub_I &b, VecDoub_O &x) {
	int i,j;
	Doub sum;
	if (sing) throw("attempting solve in a singular QR");
	for (i=0;i<m;i++) {
		sum=b[i];
		for (j=i-1;j>=0;j--) sum -= r[j][i]*x[j];
		x[i]=sum/r[i][i];
	}
	for (i=m;i<n;i++) { x[i]=0.0; }
}

void QRdcmp::update(VecDoub_I &u, VecDoub_I &v) {
	Int i,k;
	VecDoub w(u);
	for (k=n-1;k>=0;k--)
		if (w[k] != 0.0) break;
	if (k < 0) k=0;
	for (i=k-1;i>=0;i--) {
		rotate(i,w[i],-w[i+1]);
		if (w[i] == 0.0)
			w[i]=abs(w[i+1]);
		else if (abs(w[i]) > abs(w[i+1]))
			w[i]=abs(w[i])*sqrt(1.0+SQR(w[i+1]/w[i]));
		else w[i]=abs(w[i+1])*sqrt(1.0+SQR(w[i]/w[i+1]));
	}
	for (i=0;i<n;i++) r[0][i] += w[0]*v[i];
	for (i=0;i<k;i++)
		rotate(i,r[i][i],-r[i+1][i]);
	for (i=0;i<n;i++)
		if (r[i][i] == 0.0) sing=true;
}

void QRdcmp::rotate(const Int i, const Doub a, const Doub b)
{
	Int j;
	Doub c,fact,s,w,y;
	if (a == 0.0) {
		c=0.0;
		s=(b >= 0.0 ? 1.0 : -1.0);
	} else if (abs(a) > abs(b)) {
		fact=b/a;
		c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
		s=fact*c;
	} else {
		fact=a/b;
		s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
		c=fact*s;
	}
	for (j=i;j<n;j++) {
		y=r[i][j];
		w=r[i+1][j];
		r[i][j]=c*y-s*w;
		r[i+1][j]=s*y+c*w;
	}
	for (j=0;j<n;j++) {
		y=qt[i][j];
		w=qt[i+1][j];
		qt[i][j]=c*y-s*w;
		qt[i+1][j]=s*y+c*w;
	}
}


#endif /* _QRDCMP_H_ */
