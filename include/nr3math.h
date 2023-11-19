#ifndef _NR3MATH_H_
#define _NR3MATH_H_

/*--------------------------------------------------------/
// we add in custom inline routines below;
/--------------------------------------------------------*/

//#include <cmath>  // already included in nr3.h
#include <float.h>
#include "nr3.h"

inline auto VECNORM(VecDoub &sx, Int norm) -> Doub {
	Int i, n=sx.size();
	Doub ans;
	switch(norm) {
	case 1:
		for (ans=0.0,i=0;i<n;i++) { ans += abs(sx[i]); }
		break;
	case 2:
		for (ans=0.0,i=0;i<n;i++) { ans += SQR(sx[i]); }
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
// fin custom inline routines;
/--------------------------------------------------------*/


#endif /* _NR3MATH_H_ */
