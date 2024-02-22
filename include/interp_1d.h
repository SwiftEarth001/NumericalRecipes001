#ifndef _INTERP_1D_H_
#define _INTERP_1D_H_

/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/

#include "nr3.h"
#include "vander.h"

static const double PI  =3.141592653589793238463;
static const float  PI_F=3.14159265358979f;

class Poly_Base_Interp
{
public:
	const Doub *xx, *yy;

	/* constructors and destructor */
	Poly_Base_Interp(VecDoub_I &x, const Doub *y, Int m)
		: n(x.size()), mm(m), M(m+1), jsav(0), cor(0), xx(&x[0]), yy(y) 
	{
		dj = MAX(1,(int)pow((Doub)n,0.25));
	}
	Poly_Base_Interp(VecDoub_I &x, Int m)
		: n(x.size()), mm(m), M(m+1), jsav(0), cor(0), xx(&x[0]) 
	{
		dj = MAX(1,(int)pow((Doub)n,0.25));
	}
	~Poly_Base_Interp() { };  // no need for deallocation; 
	
	/* for finding local low index */
	Int findlowindex(Doub x) 
	{
		return (cor ? hunt(x) : locate(x));
	}
	Int locate(const Doub x);
	Int hunt(const Doub x);

	/* run interpolation routine */
	Doub interp(Int jlo, Doub x) 
	{
		return rawinterp(jlo, x);
	}
	Doub linearinterp(Int j, Doub x) {
		if (xx[j]==xx[j+1]) { 
            return yy[j]; 
        } else {
            return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
        }
	}

	/* static interpolation routines */
    static void uniform_mesh(Doub a, Doub b, VecDoub_IO &xv)
    {
        int num = xv.size();
        Doub h = (b-a)/(num-1);
        xv[0]=a; xv[num-1]=b;
        for (int i=1; i<num-1; i++) {
            xv[i]=a+i*h;
        }
    }
    static void Chebyshev_mesh(VecDoub_IO &xv)
    {
        int num = xv.size();
        for (int i=0; i<num; i++) {
            xv[i] = cos( (2*i+1)*PI/(2*(num-1)+2) );
        }
    }

protected:
	Int n, mm, M, jsav, cor, dj;
	Doub virtual rawinterp(Int jlo, Doub x) = 0;
};

Int Poly_Base_Interp::locate(const Doub x)
{
	Int ju, jm, jl;
	if (n < 2 || M < 2 || M > n) throw("locate size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	jl=0;
	ju=n-1;
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-M,jl-((M-2)>>1)));
}


Int Poly_Base_Interp::hunt(const Doub x)
{
	Int jl=jsav, jm, ju, inc=1;
	if (n < 2 || M < 2 || M > n) throw("hunt size error");
	Bool ascnd=(xx[n-1] >= xx[0]);
	if (jl < 0 || jl > n-1) {
		jl=0;
		ju=n-1;
	} else {
		if (x >= xx[jl] == ascnd) {
			for (;;) {
				ju = jl + inc;
				if (ju >= n-1) { ju = n-1; break;}
				else if (x < xx[ju] == ascnd) break;
				else {
					jl = ju;
					inc += inc;
				}
			}
		} else {
			ju = jl;
			for (;;) {
				jl = jl - inc;
				if (jl <= 0) { jl = 0; break;}
				else if (x >= xx[jl] == ascnd) break;
				else {
					ju = jl;
					inc += inc;
				}
			}
		}
	}
	while (ju-jl > 1) {
		jm = (ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	cor = abs(jl-jsav) > dj ? 0 : 1;
	jsav = jl;
	return MAX(0,MIN(n-M,jl-((M-2)>>1)));
}


/*----------------------------------------------------------------------------/
// ...
/----------------------------------------------------------------------------*/


struct Rat_interp : Poly_Base_Interp
{
	Doub dy;
	Rat_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
		: Poly_Base_Interp(xv,&yv[0],m), dy(0.) {}
	Doub rawinterp(Int jl, Doub x);
};


Doub Rat_interp::rawinterp(Int jl, Doub x)
{
	const Doub TINY=1.0e-99;
	Int m,i,ns=0;
	Doub y,w,t,hh,h,dd;
	const Doub *xa = &xx[jl], *ya = &yy[jl];
	VecDoub c(mm),d(mm);
	hh=abs(x-xa[0]);
	for (i=0;i<mm;i++) {
		h=abs(x-xa[i]);
		if (h == 0.0) {
			dy=0.0;
			return ya[i];
		} else if (h < hh) {
			ns=i;
			hh=h;
		}
		c[i]=ya[i];
		d[i]=ya[i]+TINY;
	}
	y=ya[ns--];
	for (m=1;m<mm;m++) {
		for (i=0;i<mm-m;i++) {
			w=c[i+1]-d[i];
			h=xa[i+m]-x;
			t=(xa[i]-x)*d[i]/h;
			dd=t-c[i+1];
			if (dd == 0.0) throw("Error in routine ratint");
			dd=w/dd;
			d[i]=c[i+1]*dd;
			c[i]=t*dd;
		}
		y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
	}
	return y;
}


struct BaryRat_interp : Poly_Base_Interp
{
	VecDoub w;
	Int d;
	BaryRat_interp(VecDoub_I &xv, VecDoub_I &yv, Int dd);
	Doub rawinterp(Int jl, Doub x);
	Doub interp(Doub x);
};


BaryRat_interp::BaryRat_interp(VecDoub_I &xv, VecDoub_I &yv, Int dd)
		: Poly_Base_Interp(xv,&yv[0],xv.size()), w(n), d(dd)
{
	if (n<=d) throw("d too large for number of points in BaryRat_interp");
	for (Int k=0;k<n;k++) {
		Int imin=MAX(k-d,0);
		Int imax = k >= n-d ? n-d-1 : k;
		Doub temp = imin & 1 ? -1.0 : 1.0;
		Doub sum=0.0;
		for (Int i=imin;i<=imax;i++) {
			Int jmax=MIN(i+d,n-1);
			Doub term=1.0;
			for (Int j=i;j<=jmax;j++) {
				if (j==k) continue;
				term *= (xx[k]-xx[j]);
			}
			term=temp/term;
			temp=-temp;
			sum += term;
		}
		w[k]=sum;
	}
}


Doub BaryRat_interp::rawinterp(Int jl, Doub x)
{
	Doub num=0,den=0;
	for (Int i=0;i<n;i++) {
		Doub h=x-xx[i];
		if (h == 0.0) {
			return yy[i];
		} else {
			Doub temp=w[i]/h;
			num += temp*yy[i];
			den += temp;
		}
	}
	return num/den;
}


Doub BaryRat_interp::interp(Doub x) {
	return rawinterp(1,x);
}


#endif /* _INTERP_1D_H_ */
