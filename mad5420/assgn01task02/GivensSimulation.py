# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:44:54 2023

@author: jackd
"""
#--------------------------------
# ...comments...
#--------------------------------

class GivensRotation:
    """
    A Givens rotation, which can act on a matrix or column vector 
    from the left;
    """
    def __init__(self, row1, row2, c=1, s=1):
        self.row1 = row1
        self.row2 = row2
        self.c = c
        self.s = s
        
    def multiply(self, mat, col, start=0, stop=0):
        from math import sqrt
        # stop takes the convention of one more than the final index;
        n, m = mat.shape
        stop = min(stop,m) if (stop > 0) else m
        if (n < max(self.row1, self.row2)):
            raise RuntimeError("incompatible sizes of operands")
        r = sqrt(mat[self.row1][col]**2 + mat[self.row2][col]**2)
        self.c = mat[self.row1][col]/r
        self.s = -mat[self.row2][col]/r
        for j in range(start, stop):
            a, b = mat[self.row1][j], mat[self.row2][j]
            mat[self.row1][j] = self.c*a - self.s*b
            mat[self.row2][j] = self.s*a + self.c*b
    
    def transpose(self):
        self.s = -self.s
### ~~~ Givens Rotation ~~~ ###


def banded_matrix_solve(M, width, b=None, start=0):
    """
    Performs a QR factorization of a banded matrix M
    with band width, and also accordingly transforms
    a column vector b if provided;
    """
    n,m = M.shape
    for j in range(start, m-1):  # do not act on last column
        for w in range(1, width+1):
            g = GivensRotation(j,j+w)
            g.multiply(M, j, start=j, stop=j+2*w+1)
            if (b != None):
                g.multiply(b, 0, start=0, stop=1)
### ~~~ banded_matrix_solve ~~~ ###


def main():
    """
    Simulation routine;
        ?(a discrepancy principle with spectral cutoff is needed)!;
        ?(or a Tichonov regularization, or both)!;
        (we use Tichonov regularization)!;
    """
    from numpy import array, zeros, eye, matmul, append, linspace, float64, sin, cos, inf
    from numpy.linalg import norm
    from scipy.sparse import diags 
    from scipy.linalg import solve_triangular
    
    timegrid = linspace(0, 4, 100)
    xtrue = sin(timegrid) + timegrid*cos(timegrid)**2  # signal xtrue
    # the "precision" of our time grid is on the
    #   scale of 0.0404, and the resulting precision 
    #   of the simulated signal happens to be on the
    #   scale of 0.0904;
    # the least value of xtrue is roughly 0.0807;
    # a spectral cutoff would cut off within this 
    #   resulting precision;
    # noise variances would ideally be within 0.0904; 
    # simulated variances should NOT be disproportionate 
    #   to 0.0904 for the method to make sense;
    # the simulation data is on a different file,
    #   noisySignal.txt;

    Tffparam_vals = [0.0, 0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 1.5, 2.0]
    Gssparam_vals = [0.00000000    , 5.02512563e-03, 5.52763819e-02, 1.05527638e-01, 
                     1.55778894e-01, 2.06030151e-01, 2.56281407e-01, 3.36532663e-01,
                     3.87035176e-01, 4.57286432e-01, 5.07537688e-01, 5.57788945e-01,
                     6.08040201e-01, 6.58291457e-01, 7.08542714e-01, 7.58793970e-01,
                     8.09045226e-01, 8.59296482e-01, 9.09547739e-01, 9.59798995e-01,
                     1.00000000e+00, 1.50000000e+00, 2.10714286e+00, 2.71428571e+00,
                     3.32142857e+00, 3.92857143e+00, 4.53571429e+00, 5.14285714e+00,
                     5.75000000e+00, 6.35714286e+00, 6.96428571e+00, 7.57142857e+00,
                     8.17857143e+00, 8.78571429e+00, 9.39285714e+00, 1.00000000e+01]
    best_params = []
    
    # read simulation data from file;
    with open("noisySignal.txt", 'r') as infile:
        for line in infile.readlines():
            floats = [float(i) for i in line.split(" ") if i.strip()]
            floats = array(floats, dtype=float64)
            normres = zeros((9, 36))
            for j,lambda_param in enumerate(Gssparam_vals):
                for i,alpha_param in enumerate(Tffparam_vals):
                    b = floats[:3]
                    bvec = ( array([b]) ).transpose()
                    M = diags([1,-1],[0,1], shape=(3,3)).toarray()
                    M = lambda_param*matmul(M.transpose(), M) + eye(3)
                    # ...insert Tichonov regularizer or spectral cutoff here...
                    M = M + alpha_param*eye(3)
                    # the delay in Givens rotation is according to
                    #   band width of M, which in this case is 1;
                    g = GivensRotation(0,1)
                    g.multiply(M, 0, 0, 1)
                    g.multiply(bvec, 0, start=0, stop=1)
                    for n in range(4, len(floats)):
                        # append new row and column to M;
                        newrow = zeros((1,n-1))
                        M = append(M, newrow, axis=0)
                        newcol = zeros((n,1))
                        M = append(M, newcol, axis=1)
                        M[n-1][n-2]=-1; 
                        M[n-1][n-1]=3+alpha_param;
                        M[n-2][n-1]=-1;
                        # append new data to bvec;
                        bvec = append(bvec, array([[ floats[n-1] ]]), axis=0)
                        # apply Givens rotation to the third last column;
                        g = GivensRotation(n-3,n-2)
                        g.multiply(M, n-3, n-3, n-1)
                        g.multiply(bvec, 0, start=0, stop=1)                
                    ## end for loop 4 ##
                    # one last Givens rotation is needed;
                    n = len(floats)-1
                    g = GivensRotation(n-2,n-1)
                    g.multiply(M, n-2, n-2, n-1)
                    g.multiply(bvec, 0, start=0, stop=1)
                
                    # solve resulting right triangular system;
                    x = solve_triangular(M, bvec)
                    # add residual of r=x-xtrue to residuals list;
                    normres[i][j] = norm(x-xtrue, ord=2)
                ## end for loop 3 ##
            ## end for loop 2 ##
            
            # get the best fitting parameter alpha given lambda;
            best_alphas = normres.argmin(axis=0)
            best_alphas = [ Tffparam_vals[id] for id in best_alphas ]

            # get the best fitting parameter lambda given alpha;
            best_lambdas = ( normres.argmin(axis=1) ).transpose()
            best_lambdas = [ Gssparam_vals[id] for id in best_lambdas ]
            
            # get the parameters with least residual norm;
            minval = inf; minid = (0,0);
            for i in range(len(Tffparam_vals)):
                for j in range(len(Gssparam_vals)):
                    minval = min(normres[i][j], minval)
                    if (minval == normres[i][j]):
                        minid = (i,j)
            best_params += [(minval, 
                             Tffparam_vals[minid[0]], 
                             Gssparam_vals[minid[1]])]
            
            # print best_alphas to file;
            outfile = open("resultsAlphaParamsl2.txt", "a")
            for alph in best_alphas:
                outfile.write(str(alph) + " ")
            outfile.write("\n")
            outfile.close()
            
            # print best_lambdas to file;
            outfile = open("resultsLambdaParamsl2.txt", "a")
            for lamb in best_lambdas:
                outfile.write(str(lamb) + " ")
            outfile.write("\n")
            outfile.close()
        ## end for loop 1 ##
        
        # print best parameters to file;
        outfile = open("resultsBestParamsl2.txt", 'w')
        for param in best_params:
            outfile.write(str(param[0]) + 
                          " " + str(param[1]) + 
                          " " + str(param[2]) + "\n")
        outfile.close()
    ## end file read infile ##
### ~~~ fin main ~~~ ###

__name__ == "__main__"

if (__name__ == "__main__"):
    main()
                
                
                
                
            
    
    
    
        