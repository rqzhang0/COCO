// su3.h
// Ben Gamari
// August 2009

#ifndef _SU3_H
#define _SU3_H

#include <string.h>
#include "complex.h"

namespace qcd {

template<typename T>
struct is_scalar {const static bool val = is_numeric<T>::val;};

template<typename T>
struct is_scalar<complex<T> > {const static bool val = is_numeric<T>::val;};


template<typename RT>
struct generic_su3_vector {
  complex<RT> c[3];

  generic_su3_vector() {}

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_su3_vector<RT> >::type 
  operator*(const T& a) const {
    generic_su3_vector<RT> res;
    for (int i=0; i<3; i++)
      res.c[i] = this->c[i] * a;
    return res;
  }

  generic_su3_vector<RT> operator+(const generic_su3_vector<RT>& a) const {
    generic_su3_vector<RT> res;
    for (int i=0; i<3; i++)
      res.c[i] = this->c[i] + a.c[i];
    return res;
  }

  generic_su3_vector<RT> operator-(const generic_su3_vector<RT>& a) const {
    generic_su3_vector<RT> res;
    for (int i=0; i<3; i++)
      res.c[i] = this->c[i] - a.c[i];
    return res;
  }

  generic_su3_vector<RT>& operator+=(const generic_su3_vector<RT>& a) {
    for (int i=0; i<3; i++)
      this->c[i] += a.c[i];
    return *this;
  }

  generic_su3_vector<RT>& operator-=(const generic_su3_vector<RT>& a) {
    for (int i=0; i<3; i++)
      this->c[i] -= a.c[i];
    return *this;
  }

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_su3_vector<RT>&>::type
  operator*=(const T& a) {
    for (int i=0; i<3; i++)
      this->c[i] *= a;
    return *this;
  }
};

typedef generic_su3_vector<double> su3_vector;

template<typename RT>
struct generic_su3_matrix
{
  complex<RT> e[3][3];

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_su3_matrix<RT> >::type
  operator*(const T& a) const 
  {
    generic_su3_matrix<RT> res;
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      res.e[i][j] = this->e[i][j] * a;
    return res;
  }



  generic_su3_matrix<RT> operator+(const generic_su3_matrix<RT>& a) const {
    generic_su3_matrix<RT> res;
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      res.e[i][j] = this->e[i][j] + a.e[i][j];
    return res;
  }

  generic_su3_matrix<RT> operator-(const generic_su3_matrix<RT>& a) const {
    generic_su3_matrix<RT> res;
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      res.e[i][j] = this->e[i][j] - a.e[i][j];
    return res;
  }

  generic_su3_matrix<RT>& operator+=(const generic_su3_matrix<RT>& a) {
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      this->e[i][j] += a.e[i][j];
    return *this;
  }

  generic_su3_matrix<RT>& operator-=(const generic_su3_matrix<RT>& a) {
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      this->e[i][j] -= a.e[i][j];
    return *this;
  }

  template<typename ORT>
  generic_su3_matrix<RT>& operator=(const generic_su3_matrix<ORT>& a) {
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      this->e[i][j] = a.e[i][j];
    return *this;
  }

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_su3_matrix<RT>&>::type
  operator*=(const T& a) 
  {
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      this->e[i][j] *= a;
    return *this;
  }

  // note this is a left multiplication i.e. M *= A gives MA. 
  generic_su3_matrix<RT>& operator*=(const generic_su3_matrix<RT>& a)
  {
    for(int r=0; r<3; ++r)
    {
      complex<RT> row[3];
      memcpy(row, &e[r], sizeof(complex<RT>)*3);
      for(int c=0; c<3; ++c)
      {
	e[r][c] = 0;
	for(int k=0; k<3; ++k) e[r][c] += row[k]*a.e[k][c];
      }
    }
    return *this;
  }

  complex<RT> trace() const 
  { 
    return this->e[0][0] + this->e[1][1] + this->e[2][2]; 
  }
};

typedef generic_su3_matrix<double> su3_matrix;

template<typename RT>
generic_su3_vector<RT> operator*(const generic_su3_matrix<RT>& a, const generic_su3_vector<RT>& b);

template<typename RT>
generic_su3_matrix<RT> adj(const generic_su3_matrix<RT>& a);


template<typename T, typename RT>
typename enable_if<is_scalar<T>::val, generic_su3_matrix<RT> >::type
operator*(const T& a, const generic_su3_matrix<RT>& b) { return b*a; }; 

template<typename RT>
generic_su3_matrix<RT> operator*(const generic_su3_matrix<RT>& a, const generic_su3_matrix<RT>& b);

template<typename RT>
struct generic_wilson_vector {
  generic_su3_vector<RT> d[4];

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_wilson_vector<RT> >::type
  operator*(const T& a) {
    generic_wilson_vector<RT> res;
    for (int i=0; i<4; i++)
      res.d[i] = this->d[i] * a;
    return res;
  }

  generic_wilson_vector<RT> operator+(const generic_wilson_vector<RT>& a) const {
    generic_wilson_vector<RT> res;
    for (int i=0; i<4; i++)
      res.d[i] = this->d[i] + a.d[i];
    return res;
  }

  generic_wilson_vector<RT> operator-(const generic_wilson_vector<RT>& a) const {
    generic_wilson_vector<RT> res;
    for (int i=0; i<4; i++)
      res.d[i] = this->d[i] - a.d[i];
    return res;
  }

  generic_wilson_vector<RT>& operator+=(const generic_wilson_vector<RT>& a) {
    for (int i=0; i<4; i++)
      this->d[i] += a.d[i];
    return *this;
  }

  generic_wilson_vector<RT>& operator-=(const generic_wilson_vector<RT>& a) {
    for (int i=0; i<4; i++)
      this->d[i] -= a.d[i];
    return *this;
  }

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_wilson_vector<RT>&>::type
  operator*=(const T& a) {
    for (int i=0; i<4; i++)
      this->d[i] *= a;
    return *this;
  }

// add by Yi-Bo Yang, for FFT
  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_wilson_vector<RT> >::type
  operator*(const T& a) const{
    generic_wilson_vector<RT> res;
    for (int i=0; i<4; i++)
      res.d[i]= this->d[i] * a;
    return res;
  }
};
// add by Yi-Bo Yang, for FFT
template<typename T, typename RT>
typename enable_if<is_scalar<T>::val, generic_wilson_vector<RT> >::type
operator*(const T& a, const generic_wilson_vector<RT>& b) { return b*a; }; 

typedef generic_wilson_vector<double> wilson_vector;

template<typename RT>
struct generic_half_wilson_vector {
  generic_su3_vector<RT> h[2];

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_half_wilson_vector<RT> >::type
  operator*(const T &a) {
    generic_half_wilson_vector<RT> res;
    for (int i=0; i<2; i++)
      res.h[i] = this->h[i] * a;
    return res;
  }

  template<typename T>
  typename enable_if<is_scalar<T>::val, generic_half_wilson_vector<RT>&>::type
  operator*=(const T &a) {
    for (int i=0; i<2; i++)
      this->h[i] *= a;
    return *this;
  }

  generic_half_wilson_vector<RT> operator+(generic_half_wilson_vector<RT> a) {
    generic_half_wilson_vector<RT> res;
    for (int i=0; i<2; i++)
      res.h[i] = this->h[i] + a.h[i];
    return res;
  }

  generic_half_wilson_vector<RT> operator-(generic_half_wilson_vector<RT> a) {
    generic_half_wilson_vector<RT> res;
    for (int i=0; i<2; i++)
      res.h[i] = this->h[i] - a.h[i];
    return res;
  }
};

typedef generic_half_wilson_vector<double> half_wilson_vector;

/* SU(2) */
typedef struct { double_complex e[2][2]; } su2_matrix;

/*
 * Wilson vectors
 * e.g.
 * wilson_propagator prop;
 * prop.c[ci].d[si].d[sf].c[cf]
 * ----------------------->    complex
 * ----------------->          su3_vector
 * ----------->                wilson_vector
 * ----->                      spin_wilson_vector
 *
 * e.g.
 * wilson_matrix matr;
 * matr.d[si].c[ci].d[sf].c[cf]
 * ----------------------->    complex
 * ----------------->          su3_vector
 * ----------->                wilson_vector
 * ----->                      color_wilson_vector
 */
typedef struct { wilson_vector c[3]; } color_wilson_vector;
typedef struct { wilson_vector d[4]; } spin_wilson_vector;
typedef struct { color_wilson_vector d[4]; } wilson_matrix;
typedef struct { spin_wilson_vector c[3]; } wilson_propagator;

// in exp_map.cpp
template<typename RT>
complex<RT> det(const generic_su3_matrix<RT>& a);

// compute exp(i Q) when Q is hermitian traceless
template<typename RT>
generic_su3_matrix<RT> exp_map(const generic_su3_matrix<RT>& Q);

// compute log(Q) when Q is unitary
template<typename RT>
generic_su3_matrix<RT> log_map(const generic_su3_matrix<RT>& a);

} // qcd namespace

#endif /* _SU3_H */


