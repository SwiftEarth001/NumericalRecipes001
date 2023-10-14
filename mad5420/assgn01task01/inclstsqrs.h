#ifndef _INCLSTSQRS_H_
#define _INCLSTSQRS_H_

/*------------------------------------------------------/
// ...comments...
/------------------------------------------------------*/

#include "../../include/nr3.h"
#include "../../include/qrdcmp.h"
#include "../../include/cholesky.h"


class IncLeastSqrs 
{
private:
    // window parameters;
    // total rows and columns of window;
    int nrows;
    int ncols;
    // moving window states;
    unsigned int start_offset;
    int ndatarows;

    // model and data;
    double** model;
    // pending model and data
    std::vector< std::vector<double> > pending;

    // solution vector and Cholesky factorization;
    double* xvar;
    struct Cholesky* chelsky;

    // private routines;
    void update_data_to_model();
    void update_chelsky_from_model();
    void inc_update_chelsky();
public:
    IncLeastSqrs();
    IncLeastSqrs(int n, int m);
    IncLeastSqrs(int n, int m, double** A, int curr_size);
    void new_data(std::vector<double> &vec);
    void clear_pending();
    void solve();
    double* get_least_sqr();
    ~IncLeastSqrs();
};

//--------------------------------------
// constructor and destructor routines
//--------------------------------------

IncLeastSqrs::IncLeastSqrs() : 
    model(NULL), xvar(NULL), chelsky(NULL) 
{
    nrows = 0;
    ncols = 0;
    start_offset = 0;
    ndatarows = 0;
}

IncLeastSqrs::IncLeastSqrs(int n, int m) : 
    nrows(n), ncols(m), xvar(NULL), chelsky(NULL)
{    
    start_offset = 0;
    ndatarows = 0;
    model = new double*[n];
    for (int k=0; k<n; k++) { model[k] = new double[m]; }
}

IncLeastSqrs::IncLeastSqrs(int n, int m, double** A, int curr_size) :
    nrows(n), ncols(m), model(A), ndatarows(curr_size), xvar(NULL)
{
    start_offset = 0;
    update_chelsky_from_model();
}

IncLeastSqrs::~IncLeastSqrs() 
{
    // delete double array for model
    if (model != NULL) {
        for (int i=0; i<nrows; i++) {
            delete[] model[i];
        }
        delete[] model;
    }
    // delete xvar
    if (xvar != NULL) { delete[] xvar; }
    // delete chelsky
    if (chelsky != NULL) { delete chelsky; }
}

//-------------------------------------------------
// routines for updating model and pending queue
//-------------------------------------------------

void IncLeastSqrs::new_data(std::vector<double> &vec) 
{
    pending.push_back(vec);
}

void IncLeastSqrs::clear_pending() 
{
    while (pending.size() > 0) {
        pending.pop_back();
    }
}

void IncLeastSqrs::update_data_to_model()
{
    for (unsigned int k=0; k<pending.size(); k++) {
        // add new row to model window;
        unsigned int index = (start_offset + k) % (nrows);
        for (int j=0; j<ncols; j++) { model[index][j] = pending[k][j]; }
        // update count of data rows and starting offset;
        ndatarows = std::min(nrows, ndatarows+1);
        if (index <= start_offset) { start_offset = index; }
    }
}

//----------------------------------------------
// routine for updating Cholesky factorization
//----------------------------------------------

void IncLeastSqrs::update_chelsky_from_model() 
{
    // column length of matrix is ncols-1,
    //   since the last element is reserved
    //   for the data vector;
    int m = ncols-1;
   
    // calculate A'*A;
    // the starting offset should not
    //   matter here, since shifting
    //   the rows amounts to a sequence
    //   of row permutations which
    //   cancel out when multiplying
    //   from the left by mat';
    // there is also the peculiarity that
    //   the rows cycle after the data
    //   rows exceed the row size, but
    //   the number of data rows should
    //   cut off to the row size once that
    //   happens;
    double* symmat = new double[(m*m)];
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            symmat[i*m+j] = 0;
            for (int k=0; k<ndatarows; k++) {
                symmat[i*m+j] += model[k][i]*model[k][j];
            }
        }
    }

    // initialize Cholesky matrix from mat'*mat;
    MatDoub model_matrix(m, m, symmat);
    chelsky = new struct Cholesky(model_matrix);

    // delete heap array symmat;
    delete[] symmat;
}

//----------------------------------------
// routine for incremental least squares
//----------------------------------------

void IncLeastSqrs::inc_update_chelsky() 
{
    // write Cholesky lower matrix to array;
    int N = chelsky->el.nrows();
    int s = pending.size();
    int m = ncols-1;
    
    double* mat = new double[(N+s)*(m)];
    for (int i=0; i<N; i++) {
        for (int j=0; j<m; j++) {
            mat[i*m+j] = chelsky->el[j][i]; 
        }
    }
    for (int i=N; i<N+s; i++) {
        int ii = i-N;
        for (int j=0; j<m; j++) {
            mat[i*m+j] = pending[ii][j];
        }
    }

    // perform Householder factorization;
    QRdcmp qrmat(N+s, m, mat);
    
    // write r' to Cholesky matrix;
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            (chelsky->el)[i][j] = (qrmat.r)[j][i];
        }
    }

    // delete mat above;
    delete[] mat;
}

void IncLeastSqrs::solve() 
{
    inc_update_chelsky();   
    update_data_to_model();
    clear_pending();
    
    // write data to vector B=A'*b;
    VecDoub B(ncols-1);
    for (int i=0; i<ncols-1; i++) {
        B[i]=0;
        for (int j=0; j<ndatarows; j++) {
            // the offset should not matter here,
            //   (and nor in calculating A'*A);
            B[i] += model[j][i]*model[j][ncols-1];
        }
    }
    // initialize solution vector X;
    VecDoub_O X(ncols-1);

    // initialize new solution vector xvar;
    if (xvar != NULL) { delete[] xvar; }  // delete any old xvar
    xvar = new double[ncols-1];
    // solve for xvar via forward-backward substitution;
    chelsky->solve((const VecDoub_I)(B), X);
    for (int k=0; k<ncols-1; k++) { xvar[k] = X[k]; }
}

double* IncLeastSqrs::get_least_sqr() 
{
    return (xvar);
}


#endif /* _INCLSTSQRS_H_ */