#ifndef _INTERP_POLY_H_
#define _INTERP_POLY_H_

/*--------------------------------------------------------------------------------------------------/
// ...
/--------------------------------------------------------------------------------------------------*/

#include "interp_1d.h"

/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/

class LagrangePolyInterp : public Poly_Base_Interp
{
public:
    /* constructor and destructor */
    LagrangePolyInterp(VecDoub_I &xv, VecDoub_I &yv, Int m)
        : Poly_Base_Interp(xv, &yv[0], m) 
    {
        barymx = MatDoub(m+1, m+1, 0.0);
    }
    ~LagrangePolyInterp() { }

    /* get methods */
    void get_inverse_barycweights(VecDoub_IO &b)
    {
        if (b.size() != mm+1) throw("size mismatch in get_barymx");
        for (int i=0; i<=mm; i++) {
            b[i] = barymx[mm][i];
        }
    }

    /* interpolation routines */
    Doub Lagrange_interp(Int jl, Doub x);
    Doub Lagrange_interp_alt(Int jl, Doub x);

private:
    MatDoub barymx;
    Int jsavl=-1;

    inline void update_barymx(const Doub* xa)
    {
        Doub p;
        barymx[0][0] = 1;
        barymx[1][0] = xa[0]-xa[1];
        barymx[1][1] = xa[1]-xa[0];
        for (int k=2; k<=mm; k++) {
            p = 1.0;
            for (int j=0; j<k; j++) {
                barymx[k][j] = barymx[k-1][j] * (xa[j]-xa[k]);
                p = -(xa[j]-xa[k]) * p;
            }
            barymx[k][k] = p;
        }
    }

protected:
    // virtual method rawinterp;
    Doub rawinterp(Int jl, Doub x)
    {
        return Lagrange_interp(jl, x);
    }
};


Doub LagrangePolyInterp::Lagrange_interp(Int jl, Doub x)
{
    const Doub *xa=&xx[jl], *ya=&yy[jl];
    if (jl != jsavl) {
        update_barymx(xa);
        jsavl=jl;
    }

    Doub sumnum=0.0, sumden=0.0, lx;
    for (int i=0; i<=mm; i++) {
        if (x==xa[i]) { return ya[i]; }
        lx = barymx[mm][i]*(x-xa[i]);
        sumnum += ya[i] / lx;
        sumden += 1.0 / lx;
    }

    return (sumnum / sumden);
}


Doub LagrangePolyInterp::Lagrange_interp_alt(Int jl, Doub x)
{
    const Doub *xa=&xx[jl], *ya=&yy[jl];
    if (jl != jsavl) {
        update_barymx(xa);
        jsavl=jl;
    }

    Doub sumnum=0.0, w=1.0, lx;
    for (int i=0; i<=mm; i++) {
        if (x==xa[i]) { return ya[i]; }
        lx = barymx[mm][i]*(x-xa[i]);
        sumnum += ya[i] / lx;
        w *= x-xa[i];
    }

    return (w*sumnum);
}

/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/


class NewtonPolyInterp : public Poly_Base_Interp
{
public:
    /* constructor and destructor */
    NewtonPolyInterp(VecDoub_I &xv, VecDoub_I &yv, Int m)
        : Poly_Base_Interp(xv, &yv[0], m), dy(0.0) 
    {
        divdiff = MatDoub(m+1, m+1, 0.0);
    }
    ~NewtonPolyInterp() { }

    /* get methods */
    void get_divdiff(MatDoub_IO &D)
    {
        if ( (D.nrows() != mm+1) || (D.ncols() != mm+1) ) 
            throw("size mismatch in get_divdiff");
        for (int k=0; k<=mm; k++) {
            for (int j=0; j<=mm; j++) {
                D[k][j] = divdiff[k][j];
            }
        }
    }

    /* interpolation routines */
    Doub Newton_Horner_rule(Int jl, Doub x);
    Doub Neville_scheme(Int jl, Doub x);

private:
    MatDoub divdiff;
    Int jsavn=-1;

    Doub dy;  // used for the Neville scheme

    inline void update_divdiff(const Doub* xa, const Doub* ya) 
    {
        for (int j=0; j<mm; j++) {
            divdiff[0][j]=ya[j];
        }
        for (int k=1; k<=mm; k++) {
            for (int j=0; j<=mm-k; j++) {
                divdiff[k][j] = (divdiff[k-1][j+1]-divdiff[k-1][j]) / (xa[j+k]-xa[j]);
            }
        }
    }

protected:
    // virtual method rawinterp;
    Doub rawinterp(Int jl, Doub x)
    {
        return Neville_scheme(jl, x);
    }
};


Doub NewtonPolyInterp::Newton_Horner_rule(Int jl, Doub x)
{
    const Doub *xa=&xx[jl], *ya=&yy[jl];
    if (jl != jsavn) {
        update_divdiff(xa, ya);
        jsavn=jl;
    }

    Doub p=divdiff[mm][0];
    for (int k=mm-1; k>=0; k--) {
        p = p*(x-xa[k]) + divdiff[k][0];
    }

    return p;
}


Doub NewtonPolyInterp::Neville_scheme(Int jl, Doub x)
{
    Int i, m, ns=0;
    Doub y, den, dif, dift, ho, hp, w;
    const Doub *xa=&xx[jl], *ya=&yy[jl];
    VecDoub c(mm+1), d(mm+1);
    dif=abs(x-xa[0]);
    
    for (i=0; i<=mm; i++) {
        if ((dift=abs(x-xa[i])) < dif) {
            ns=i;
            dif=dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    
    y=ya[ns--];
    for (m=1; m<=mm; m++) {
        for (i=0; i<=mm-m; i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            w=c[i+1]-d[i];
            if ((den=ho-hp) == 0.0) throw("Poly_interp error");
            den = w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
    }
    return y;
}


/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/

struct Hermite_interp : Poly_Base_Interp
{
    /* constructor and destructor */
    Hermite_interp(VecDoub_I &xv, MatDoub_I &yv, VecInt_I &od, int N)
        : Poly_Base_Interp(xv, N-1), NN(N), ydat(yv), order(od)
    {
        maxorder=0; 
        int M=0;
        for (int k=0; k<od.size(); k++) { 
            M += od[k]+1;
            maxorder = (maxorder < od[k]) ? od[k] : maxorder; 
        }
        divdiff=MatDoub(M, M, 0.0);
        deg=M-1;
        xz=VecDoub(deg+1, 0.0);
    }
    Hermite_interp(VecDoub_I &xv, MatDoub_I &yv, int N, int od)
        : Poly_Base_Interp(xv, N-1), NN(N), ydat(yv)
    {
        maxorder=od;
        order = VecInt(N, maxorder);
        deg = N*(od+1)-1;
        divdiff=MatDoub(deg+1, deg+1, 0.0);
        xz=VecDoub(deg+1, 0.0);
    }
    ~Hermite_interp() { }

    /* get methods */
    Int get_degree() { return deg; }
    void get_divdiff(MatDoub_IO &D)
    {
        if ( (D.nrows() != deg+1) || (D.ncols() != deg+1) ) 
            throw("size mismatch in get_divdiff");
        for (int k=0; k<=deg; k++) {
            for (int j=0; j<=deg; j++) {
                D[k][j] = divdiff[k][j];
            }
        }
    }

    /* interpolation routines */
    Doub Newton_Hermite_rule(Int jl, Doub x)
    {
        const Doub *xa=&xx[jl]; Doub *xv=&xz[0];
        if (jl != jsavh) {
            update_divdiff(jl, xa, xv);
            jsavh=jl;
        }

        Doub p=divdiff[deg][0];
        for (int k=deg-1; k>=0; k--) {
            p = p*(x-xv[k]) + divdiff[k][0];
        }

        return p;
    }

private:
    MatDoub ydat;
    VecDoub xz;
    VecInt order;
    MatDoub divdiff;
    Int NN, maxorder, deg;
    Int jsavh=-1;

    inline void update_divdiff(Int jl, const Doub* xa, Doub* xv) 
    {
        VecInt yidx = VecInt(deg+1, 0.0);
        
        for (int idx=0, j=0; j<NN; j++) {
            for (int i=0; i<=order[j]; i++) {
                xv[idx] = xa[j];
                yidx[idx] = j+jl;
                divdiff[0][idx++]=ydat[0][j+jl];
            }
        }

        for (int fac=1, k=1; k<=maxorder; k++) {
            fac *= k;
            for (int j=0; j<=deg-k; j++) {
                if (xv[j]==xv[j+k]) {
                    divdiff[k][j]=ydat[k][yidx[j]]/fac;
                } else {
                    divdiff[k][j] = (divdiff[k-1][j+1]-divdiff[k-1][j]) / (xv[j+k]-xv[j]);
                }
            }
        }

        for (int k=maxorder+1; k<=deg; k++) {
            for (int j=0; j<=deg-k; j++) {
                divdiff[k][j] = (divdiff[k-1][j+1]-divdiff[k-1][j]) / (xv[j+k]-xv[j]);
            }
        }
    }

protected:
    // virtual method rawinterp;
    Doub rawinterp(Int jl, Doub x)
    {
        return Newton_Hermite_rule(jl, x);
    }

};



/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/

struct Spline_interp : Poly_Base_Interp
{
public:
    enum EdgeConstraints { Natural, Periodic, FirsDerivative };

    /* constructors and destructor */
	Spline_interp(VecDoub_I &xv, VecDoub_I &yv, EdgeConstraints ec)
	    : Poly_Base_Interp(xv, &yv[0], 1), y2(xv.size())
	{
        sety2(&xv[0], &yv[0], ec);
        edge=ec;
    }

	Spline_interp(VecDoub_I &xv, const Doub *yv, EdgeConstraints ec)
	    : Poly_Base_Interp(xv, yv, 1), y2(xv.size())
	{
        sety2(&xv[0], yv, ec);
        edge=ec;
    }

	void sety2(const Doub *xv, const Doub *yv, EdgeConstraints ec);

    Doub interpolate(Int jl, Doub x) { return rawinterp(jl, x); }

private:
    VecDoub y2;
    EdgeConstraints edge;

protected:
	Doub rawinterp(Int jl, Doub xv);
};


void Spline_interp::sety2(const Doub *xv, const Doub *yv, EdgeConstraints ec)
{
	Int i, k;
	Doub p, qn, sig, un;
    Doub yp1=1.0, ypn=1.0;
	VecDoub u(n-1);
	
    if (ec == Natural) {
		y2[0]=u[0]=0.0;
	} else {
		y2[0] = -0.5;
		u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
	}

	for (i=1; i<n-1; i++) {
		sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
		u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
	}

	if (ec == Natural) {
		qn=un=0.0;
	} else {
		qn=0.5;
		un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
	}

	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	
    for (k=n-2; k>=0; k--)
	{
        y2[k]=y2[k]*y2[k+1]+u[k];
    }
}


Doub Spline_interp::rawinterp(Int jl, Doub x)
{
	Int klo=jl,khi=jl+1;
	Doub y,h,b,a;
	h=xx[khi]-xx[klo];
	if (h == 0.0) throw("Bad input to routine splint");
	a=(xx[khi]-x)/h;
	b=(x-xx[klo])/h;
	y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]
		+(b*b*b-b)*y2[khi])*(h*h)/6.0;
	return y;
}


#endif /* _INTERP_POLY_H_ */