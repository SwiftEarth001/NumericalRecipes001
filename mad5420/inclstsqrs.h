#ifndef _INCLSTSQRS_H_
#define _INCLSTSQRS_H_

/*------------------------------------------------------/
// ...comments...
/------------------------------------------------------*/

#include "../include/nr3.h"
#include "../include/qrdcmp.h"
#include "../include/cholesky.h"


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
    // pending model data
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

//-----------------------
// constructor routines
//-----------------------

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
    nrows(n), ncols(m), model(A), xvar(NULL)
{
    nrows = n;
    ncols = m;
    start_offset = 0;
    ndatarows = curr_size;
    update_chelsky_from_model();
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
        double** address = (double**) ( ( (unsigned long long int) model + 
            start_offset + k ) % (nrows) );
        for (int j=0; j<ncols; j++) { (*address)[j] = pending[k][j]; }
        // update count of data rows and starting offset;
        ndatarows = std::min(nrows, ndatarows+1);
        if ((unsigned long long int)(address) <= (unsigned long long int)(model) + start_offset) {
            (start_offset)++;
        }
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
    // copy window model to matrix mat;
    double* mat = new double[(ndatarows*m)];
    for (int index=0, i=0; i<ndatarows; i++) {
        for (int j=0; j<m; j++) {
            index = (i+start_offset) % (nrows);
            mat[i*m+j] = model[index][j];
        }
    }

    // calculate mat'*mat;
    double* symmat = new double[(m*m)];
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            symmat[i*m+j] = 0;
            for (int k=0; k<ndatarows; k++) {
                symmat[i*m+j] += mat[k*m+i]*mat[k*m+j];
            }
        }
    }

    // initialize Cholesky matrix from mat'*mat;
    MatDoub* model_matrix = new MatDoub(m, m, symmat);
    chelsky = new struct Cholesky(*model_matrix);
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            printf("%f ", (chelsky->el)[i][j]);
        }
        printf("\n");
    }
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
    MatDoub_I* mmat = new MatDoub_I(N+s, m, mat);
    for (int i=0; i<N+s; i++) {
        for (int j=0; j<m; j++) {
            printf("%f ", (*mmat)[i][j]);
        }
        printf("\n");
    }
    struct QRdcmp* qrfact = new struct QRdcmp(*mmat);
    printf("%d %d\n", qrfact->n, qrfact->m);
    for (int i=0; i<N+s; i++) {
        for (int j=0; j<m; j++) {
            printf("%f ", (qrfact->r)[i][j]);
        }
        printf("\n");
    }
    
    // write r' to Cholesky matrix;
    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            (chelsky->el)[i][j] = (qrfact->r)[j][i];
        }
    }
    for (int i=0; i<N+s; i++) {
        for (int j=0; j<m; j++) {
            printf("%f ", (chelsky->el)[i][j]);
        }
        printf("\n");
    }
}

void IncLeastSqrs::solve() 
{
    inc_update_chelsky();   
    update_data_to_model();
    clear_pending();
    
    // write data to vector;
    double* B = new double[ndatarows];
    for (int index=0, i=0; i<ndatarows; i++) {
        index = (i+start_offset) % (nrows);
        B[i] = model[index][ncols-1];
    }
    VecDoub_I* b = new VecDoub(ndatarows, B);
    
    // initialize solution vector x;
    xvar = new double[ndatarows];
    VecDoub_O* x = new VecDoub(ndatarows, xvar);

    // solve for xvar via forward-backward substitution
    chelsky->solve(*b, *x);
    for (int k=0; k<ndatarows; k++) { xvar[k] = (*x)[k]; }
}

double* IncLeastSqrs::get_least_sqr() 
{
    return (xvar);
}


#endif /* _INCLSTSQRS_H_ */