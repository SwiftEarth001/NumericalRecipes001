# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 10:59:10 2023

@author: jackd
"""

#### static variable definitions ####

EPS=1.0e-07
MAXITER=500;
IVAL=list([0 for k in range(500)])
XSEQ=list([0 for k in range(500)])
FVAL=list([0 for k in range(500)])


#### numeric methods ####

def bisection(f, a, b, tol, maxiter):
    """    
    Parameters
    ----------
    f : univariate function
    a : left endpoint
    b : right endpoint
    tol : residual tolerance
    maxiter : maximum iterations
    
    Returns
    -------
    root x* of f(x)

    """
    for i,its in enumerate(range(maxiter)):
        x = (a+b)/2;
### << _ERROR_ANALYSIS_
        XSEQ[i] = x;
        FVAL[i] = f(x);
        IVAL[i] = abs(b-a);
### _ERROR_ANALYSIS_ >>
        if (abs(f(x)) < tol):
            return (x,its+1);
        if (f(x)*f(a) < 0):
            b=x;
        if (f(x)*f(b) < 0):
            a=x;
## ~ fin bisection ~ ##


def fixedpoint(f, x0, eps, maxiter):
    """    
    Parameters
    ----------
    f : univariate function
    x0 : initial point
    eps : sequential error tolerance
    maxiter : maximum iterations
    
    Returns
    -------
    fixed point x* of f(x);

    """
    for i,its in enumerate(range(maxiter)):
### << _ERROR_ANALYSIS_
        XSEQ[i] = x0;
### _ERROR_ANALYSIS_ >>
        x = f(x0);
        if (abs(x-x0) < eps):
            return (x,its+1);
        x0 = x;
## ~ fin fixepoint ~ ##


def NewtonMethod(f, x0, tol, maxiter):
    """    
    Parameters
    ----------
    f : univariate function
    x0 : initial point
    tol : residual tolerance
    maxiter : maximum iterations
    
    Returns
    -------
    root x* of f(x);

    """
    h=1.0e-10;
    for i,its in enumerate(range(maxiter)):
### << _ERROR_ANALYSIS_
        XSEQ[i] = x0;
        FVAL[i] = f(x0);
### _ERROR_ANALYSIS_ >>
        df = (f(x0+h)-f(x0-h))/(2*h);
        x = x0-(1/df)*f(x0);
        if (abs(f(x)) < tol):
            return (x,its+1);
        x0 = x;
## ~ fin fixepoint ~ ##


#### errors and residuals ####
def get_seq_errors():
    """
    Calculates and returns the sequential errors
    from the static var XSEQ;
    """
    errors = list([0 for i in range(1,500)]);
    for i in range(499):
        errors[i] = abs(XSEQ[i+1]-XSEQ[i]);
    return errors;

def get_seq_residuals():
    """
    Calculates and returns the sequential residuals
    from the static var FVAL;
    """
    residuals = list([0 for i in range(1,500)]);
    for i in range(499):
        residuals[i] = abs(FVAL[i+1]-FVAL[i]);
    return residuals;

def get_apriori_errors(xf):
    """
    Calculates and returns a priori errors from the
    static var XSEQ;
    """
    errors = list([0 for i in range(500)]);
    for i in range(500):
        errors[i] = abs(XSEQ[i]-xf);
    return errors;

def get_apriori_residuals(fp):
    """
    Calculates and returns a priori residuals from the
    static var FVAL;
    """
    residuals = list([0 for i in range(500)]);
    for i in range(500):
        residuals[i] = abs(FVAL[i]-fp);
    return residuals;
#### ~ fin errors and residuals ~ ####


#### cleanup ####
    """
    Routines for zeroing static vars XSEQ, FVAL, IVAL;
    """
def clear_xseq():
    for i,_ in enumerate(XSEQ):
        XSEQ[i]=0;

def clear_fval():
    for i,_ in enumerate(FVAL):
        FVAL[i]=0;

def clear_ival():
    for i,_ in enumerate(IVAL):
        IVAL[i]=0;

def clear_lists():
    clear_xseq();
    clear_fval();
    clear_ival();
#### ~ cleanup ~ ####


#### main method ####
def main() -> None:
    import numpy as np;
    import matplotlib.pyplot as plt;
    ## < test01 >
    ## < test02 >
    ## < test03 >
## ~ fin main ~ ##    


# __name__="__main__";
if __name__ == "__main__":
    main();
        
    