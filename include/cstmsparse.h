#ifndef _CSTMSPARSE_H_
#define _CSTMSPARSE_H_

#include "nr3.h"

typedef void (*get_t)(int, int);

template <class T>
class CSTMsparseMat {
private:
	int nn;
	int mm;
    T (*get_val)(int, int);
    inline T zero(int i, int j) const;
    inline T identity(int i, int j) const; 
public:
	CSTMsparseMat();    // assumes identity
	CSTMsparseMat(int n, int m, bool id);    // assumes identity or zero
    CSTMsparseMat(int n, int m, get_t func);
    inline T get(int i, int j) const;
	inline int nrows() const;
	inline int ncols() const;
	void resize(int newn, int newm); // resize (contents not preserved)
    typedef T value_type; // make T available externally
	~CSTMsparseMat();
};

template <class T>
CSTMsparseMat<T>::CSTMsparseMat() : nn(0), mm(0), get_val(NULL) {}

template <class T>
CSTMsparseMat<T>::CSTMsparseMat(int n, int m, bool id) 
{
    nn = n;
    mm = m;
    get_val = (id) ? this.identity : this.zero;
}

template <class T>
CSTMsparseMat<T>::CSTMsparseMat(int n, int m, get_t func) 
{
    nn = n;
    mm = m;
    get_val = func;
}

template <class T>
CSTMsparseMat<T>::~CSTMsparseMat() {
    
}

template <class T>
inline T CSTMsparseMat<T>::get(int i, int j) const
{
    return (get_val(i,j));
}

template <class T>
inline int CSTMsparseMat<T>::nrows() const 
{
    return (nn);
}

template <class T>
inline int CSTMsparseMat<T>::ncols() const 
{
    return (mm);
}

template <class T>
void CSTMsparseMat<T>::resize(int n, int m) 
{
    nn = n;
    mm = m;
}


#endif /* _CSTMSPARSE_H_ */