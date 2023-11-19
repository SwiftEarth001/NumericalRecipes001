# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 03:18:56 2023

@author: jackd
"""

## --- libraries ---
import numpy as np;
import pandas as pd;

## --- global ---
samples = [];


## --- routines ---
def get_LUdcmp(file, size):
    mat = pd.read_csv(file, nrows=size, delim_whitespace=True, header=None)
    perm = pd.read_csv(file, skiprows=size+1, nrows=2, delim_whitespace=True, header=None)
    Pperm = np.array( perm.iloc[0], dtype='int' )
    P = np.zeros((size,size), dtype='float64')
    P[Pperm, list(range(size))] = 1
    Qperm = np.array( perm.iloc[1], dtype='int' )
    Q = np.zeros((size,size), dtype='float64')
    Q[Qperm, list(range(size))] = 1
    mat = np.matmul(mat, Q)
    U = np.triu(mat, k=0)
    L = np.tril(mat, k=-1)
    L = L + np.eye(size)
    return (L, U, P, Q)
## ~ fin get_LUdcmp ~ ##

def get_errors(files, nord):
    relerrs = []; abserrs = [];
    global samples;
    for i,file,size in files:
        L, U, P, Q = get_LUdcmp(file,size)
        A = np.matmul( np.matmul( np.matmul(P, L), U ), Q.transpose() ) - np.array(samples[i])
        abserrs += [np.linalg.norm(A, ord=nord)]
        relerrs += [np.linalg.norm(A, ord=nord) / np.linalg.norm(np.array(samples[i]), ord=nord)]
        errors = np.array( [abserrs, relerrs] )
    return (errors)
## ~ fin get_errors ~ ##


## --- main ---
def main():
    import numpy as np;
    
    filenames = [ 'sample' + str(i) + '_results.txt' for i in range(6) ]
    sizes = [ 10, 25, 50, 100, 250, 500 ]
    files = list( zip( list(range(6)), filenames, sizes ) )
    testerrors_inf = get_errors(files, nord=np.inf)
    testerrors_l1 = get_errors(files, nord=1)
## ~ fin main ~ ##

# __name__ == "__main__"
if (__name__ == "__main__"):
    main()