#ifndef _NR3MATH_H_
#define _NR3MATH_H_

//#include <cmath>  // already included in nr3.h
#include <float.h>
#include "nr3.h"


/*--------------------------------------------------------/
// 1. custom structs and classes
/--------------------------------------------------------*/


/*--------------------------------------------------------/
// 2. custom inline routines below;
/--------------------------------------------------------*/

inline auto VECNORM(VecDoub_I &sx, int norm) -> double {
	int i, n=sx.size();
	double ans;
	switch(norm) {
	case 1:
		for (ans=0.0,i=0;i<n;i++) { ans += abs(sx[i]); }
		break;
	case 2:
		for (ans=0.0,i=0;i<n;i++) { ans += SQR(sx[i]); }
		ans = sqrt(ans);
		break;
	case INT_MAX:
		for (ans=0.0,i=0;i<n;i++) { 
			ans = (ans <= abs(sx[i])) ? abs(sx[i]) : ans; 
		}
		break;
	default:
		throw("illegal norm in nr3math");
		break;
	}
	return (ans);
}


/*--------------------------------------------------------/
// 2. function classes
/--------------------------------------------------------*/

template <class T>
struct scalar_function 
{
	double EPS;
	T &func;
	Doub f;
	
    // constructor
    scalar_function(T &funcc) : EPS(1.0e-10), func(funcc) {}

	// assign function
	void assign(T &funcc) 
	{ 
		func=funcc; 
	}

    // function operator (x)
	Doub operator() (VecDoub_I &x)
	{
		return f=func(x);
	}

    // function operator (x,d,alpha)
    Doub operator() (VecDoub_I &x, VecDoub_I &d, double alpha)
    {
        int n = x.size();
        VecDoub v(n);
        for (int i=0; i<n; i++) {
            v[i] = x[i] + alpha*d[i];
        }
        return func(v);
    }

    // gradient Df along Xn
	void df(VecDoub_I &x, VecDoub_O &df)
	{
		int n=df.size();
		VecDoub xh=x;
		f=func(x);
		for (int j=0; j<n; j++) {
			Doub temp=x[j];
			Doub h=EPS*abs(temp);
			if (h == 0.0) { h=EPS; }
			xh[j]=temp+h;
			Doub fh=operator()(xh);
			xh[j]=temp-h;
			Doub fl=operator()(xh);
			xh[j]=temp;
			df[j]=(fh-fl)/(2*h);
		}
	}

    // gradient Df along Xb
	void dfb(VecDoub_I &x, VecDoub_O &df)
	{
		int m=df.size(), n=x.size();
		VecDoub xh=x;
		f=func(x);
		for (int j=n-m; j<n; j++) {
			Doub temp=x[j];
			Doub h=EPS*abs(temp);
			if (h == 0.0) { h=EPS; }
			xh[j]=temp+h;
			Doub fh=operator()(xh);
			xh[j]=temp-h;
			Doub fl=operator()(xh);
			xh[j]=temp;
			df[j]=(fh-fl)/(2*h);
		}
	}

    // directional derivative at alpha
    Doub dfalpha(VecDoub_I &x, VecDoub_I &d, double alpha)
    {
        Doub fa, fnew;
        if (alpha==0) { fa=func(x); }
        else { fa=operator()(x,d,alpha); }
        double h=EPS*VECNORM(x,2);
        if (h == 0.0) { h=EPS; }
        fnew=operator()(x,d,alpha+h);
        return ( (fnew-fa)/h );
    }

    // finite difference Hessian
    void Hess(VecDoub_I &x, MatDoub_O &H)
    {
        // EPS=1.0e-10 is too small a step size when taking its square;
        double h=1.0e-06;  // step size
        double invstep=1/(h*h);
        int n=x.size();
        VecDoub ei(n,0.0), ej(n,0.0), eij(n,0.0);
        for (int i=0; i<n; i++) {
            for (int j=i; j<n; j++) {
                ei[i]=1.0; ej[j]=1.0; eij[i]=1.0; eij[j]=1.0;
                if (i==j) { eij[i]=2.0; }
                H[i][j] = invstep * ( operator()(x,eij,h)-operator()(x,ei,h)-operator()(x,ej,h)+operator()(x) );
                ei[i]=0.0; ej[j]=0.0; eij[i]=0.0; eij[j]=0.0;
                H[j][i]=H[i][j];  // Hessian is symmetric assuming C2 differentiability
            }
        }
    }
};

// definitions for multivariate real objective function;
typedef double (*objPhi)(VecDoub_I &);
typedef struct scalar_function<objPhi> ObjFuncd;

// definitions for real-valued constraint functions;
typedef double (*cnstrPhi)(VecDoub_I &);
typedef struct scalar_function<cnstrPhi> CnstrFuncd;



/*--------------------------------------------------------/
// fin custom inline routines;
/--------------------------------------------------------*/


#endif /* _NR3MATH_H_ */
