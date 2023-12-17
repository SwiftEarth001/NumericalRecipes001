#ifndef _LINSEARCH_H_
#define _LINSEARCH_H_

/*--------------------------------------------------------/
// we expand this to a class...;
/--------------------------------------------------------*/

//#include <float.h>  // already included in nr3math.h
#include <random>
#include <time.h>
#include <initializer_list>
#include "nr3math.h"


/*--------------------------------------------------------/
// 1. definitions for objective functions
/--------------------------------------------------------*/


//---------------------------------------------------------
// 2.a exact search headers
//---------------------------------------------------------

double GoldenSection(VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double lbound, double rbound);

double Fibonacci(VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double lbound, double rbound);

double QuadIntrplte(VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double lbound, double rbound);

double CubicIntrplte(VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double lbound, double rbound);


//---------------------------------------------------------
// 2.b inexact search headers
//---------------------------------------------------------

bool GoldsteinRule(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound, 
    float rho, float omega, float backrate, 
    int* iter, const int maxiter);

bool BackTracking(
    VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double* alpha, const float alpha0, float rho, 
    int* iter, int maxiter);

bool WolfePowell(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound,
    float rho, float sigma, float backrate,
    int* iter, int maxiter);

bool StrongWolfe(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound,
    float rho, float sigma, float backrate,
    int* iter, int maxiter);


//---------------------------------------------------------
// 3.a exact line search entry point
//---------------------------------------------------------

enum exact_searches {
    GoldenSectionS=0,
    FibonacciS=2,
    QuadraticInterpolationS=4,
    CubicInterpolationS=8
};

void exact_linsearch(VecDoub_I &x, VecDoub_I &d, VecDoub_O &xnew, objPhi &func, exact_searches ES,
    double* alpha, bool* check, int* iter, double stpmax, int maxiter);


//---------------------------------------------------------
// 3.b inexact line search entry point
//---------------------------------------------------------

enum inexact_searches {
    GoldsteinS=0,
    BackTrackingS=2,
    WolfePowellS=4,
    StrongWolfeS=8,
    NoLNS=256
};

typedef struct linsearch_params {
    double alpha=1.0;
    float stpmax=FLT_MAX;
    float alpha0=1.0;

    inexact_searches IES;
	float rho=0.0;
	float sigma=0.0;
	float omega=0.0;
    
    float backrate=0.9;
    float forate=0.9;

    bool check;
    int maxiter;
    int iters;
} line_search_params;

void inexact_linsearch(VecDoub_I &x, VecDoub_I &d, objPhi &func, line_search_params* params)
{
    srand(static_cast <unsigned> (time(0)));  // time is seeded on millisecond scale
    float randval = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));

    float rho=0.0, omega=0.0, sigma=0.0, alpha0, stpmax, backrate; 
    int maxiter;
    if (params == nullptr) { throw("null pointer provided to inexact line search"); }
    else { 
        rho=params->rho; sigma=params->sigma; omega=params->omega;
        alpha0=params->alpha0; stpmax=params->stpmax; backrate=params->backrate;
        maxiter=params->maxiter;
    }
    if (rho==0.0) { rho = 0.1 + 0.4 * randval; }
    
    switch(params->IES) {
    case 0:
        if (omega==0.0) { omega = 1.0 + 2.0 * randval; }
        params->check = GoldsteinRule(x, d, func, 
            &(params->alpha), alpha0, stpmax, 
            rho, omega, backrate, 
            &(params->iters), maxiter);
        break;
    case 2:
        params->check = BackTracking(x, d, func, 
            &(params->alpha), alpha0, rho,
            &(params->iters), maxiter);
        break;
    case 4:
        params->check = WolfePowell(x, d, func,
            &(params->alpha), alpha0, stpmax,
            rho, sigma, backrate, 
            &(params->iters), maxiter);
        break;
    case 8:
        params->check = StrongWolfe(x, d, func,
            &(params->alpha), alpha0, stpmax,
            rho, sigma, backrate, 
            &(params->iters), maxiter);
        break;
    case 256:
        params->alpha=1.0;
    default:
        break;
    }
}


//---------------------------------------------------------
// 4.a exact search implementation
//---------------------------------------------------------


//---------------------------------------------------------
// 4.b inexact search implementation
//---------------------------------------------------------

bool GoldsteinRule(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound, 
    float rho, float omega, float backrate, 
    int* iter, const int maxiter) 
{
    ObjFuncd funcd(func);
    float lbound=0.0;
    if (alpha0>0) { *alpha=alpha0; }
    else { *alpha=1.0; }

    double phi0 = funcd(p);
    double dphi0 = funcd.dfalpha(p,d,0);
    
    for (int i=0; i<maxiter; i++) {
        double phi_k = funcd(p,d,*alpha);
        double val = phi0 + rho*(*alpha)*dphi0;
        if (phi_k <= val) {
            val = phi0 + (1-rho)*(*alpha)*dphi0;
            if (phi_k >= val) {
                *iter=i+1;
                return true;
            } else {
                lbound = *alpha;  // leave rbound unchanged
                *alpha += omega;
            }
        } else {
            rbound = *alpha;  // leave lbound unchanged
            *alpha = (lbound + rbound)/2;
            omega *= backrate;
        }
    }
    *iter=maxiter;
    return false;
}


bool BackTracking(
    VecDoub_I &p, VecDoub_I &d, objPhi &func,
    double* alpha, const float alpha0, float rho, 
    int* iter, int maxiter)
{
    ObjFuncd funcd(func);
    if (alpha0>0) { *alpha=alpha0; }
    else { *alpha=1.0; }

    Doub phi0 = funcd(p);
    Doub dphi0 = funcd.dfalpha(p,d,0);

    for (int i=0; i<maxiter; i++) {
        Doub phi_k = funcd(p,d,*alpha);
        Doub val = phi0 + rho*(*alpha)*dphi0;
        if (phi_k <= val) {
            *iter=i+1;
            return true;
        } else {
            // new alpha determined from quadratic interpolation;
            *alpha = -dphi0 / (2*(phi_k-phi0-dphi0));
        }
    }
    *iter=maxiter;
    return false;
}


bool WolfePowell(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound,
    float rho, float sigma, float backrate,
    int* iter, int maxiter) 
{
    ObjFuncd funcd(func);
    if (alpha0>0) { *alpha=alpha0; }
    else { *alpha=1.0; }
    float lbound=0.0;
    float rate = 1.0;

    Doub phi0 = funcd(p);
    Doub dphi0 = funcd.dfalpha(p,d,0);
    
    for (int i=0; i<maxiter; i++) {
        Doub phi_k = funcd(p,d,*alpha);
        Doub val = phi0 + rho*(*alpha)*dphi0;
        if (phi_k <= val) {
            Doub dphi_k = funcd.dfalpha(p,d,*alpha);
            if (dphi_k >= sigma*dphi0) {
                *iter=i+1;
                return true;
            } else {
                lbound = *alpha;
                *alpha += rate*(*alpha);
            }
        } else {
            rbound = *alpha;
            *alpha = (*alpha+lbound) / 2;
            rate *= backrate;
        }
    }
    *iter=maxiter;
    return false;
}


bool StrongWolfe(
    VecDoub_I &p, VecDoub_I &d, objPhi &func, 
    double* alpha, const float alpha0, float rbound,
    float rho, float sigma, float backrate,
    int* iter, int maxiter) 
{
    ObjFuncd funcd(func);
    if (alpha0>0) { *alpha=alpha0; }
    else { *alpha=1.0; }
    float lbound=0.0;
    float rate = 1.0;

    Doub phi0 = funcd(p);
    Doub dphi0 = funcd.dfalpha(p,d,0);
    
    for (int i=0; i<maxiter; i++) {
        Doub phi_k = funcd(p,d,*alpha);
        Doub val = phi0 + rho*(*alpha)*dphi0;
        if (phi_k <= val) {
            Doub dphi_k = funcd.dfalpha(p,d,*alpha);
            if (dphi_k >= sigma*dphi0) {
                if (dphi_k <= -sigma*dphi0) {
                    *iter=i+1;
                    return true;
                } else {
                    rbound = *alpha;
                    *alpha = (*alpha+lbound) / 2;
                    rate *= backrate;
                }
            } else {
                lbound = *alpha;
                *alpha += rate*(*alpha);
            }
        } else {
            rbound = *alpha;
            *alpha = (*alpha+lbound) / 2;
            rate *= backrate;
        }
    }
    *iter=maxiter;
    return false;
}


#endif /* _LINSEARCH_H_ */
