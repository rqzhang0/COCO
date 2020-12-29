// complex.h
// Ben Gamari
// August 2009

#ifndef _COMPLEX_H
#define _COMPLEX_H

#include <math.h>
#include "type_support.h"

namespace qcd {

using ::sqrt;
using ::exp;
using ::pow;

template<typename T>
struct complex
{
  T real, imag;
  typedef T EType;
  /*
   * constructors
   */
  complex() : real(0.0), imag(0.0) { }
  complex(T real, T imag=0.0) : real(real), imag(imag) { }

  complex& operator=(const T& real)
  {
    this->real = real;
    this->imag = 0;
    return *this;
  }

  template <typename RT>
  complex(const complex<RT>& b): real(b.real), imag(b.imag) {}

  /*
   * unary operations
   */
  complex operator-() const
  {
    return complex(-this->real, -this->imag);
  }

  /*
   * complex-complex operations
   */
  complex operator+(const complex& b) const
  {
    return complex(real+b.real, imag+b.imag);
  }

  complex operator-(const complex& b) const
  {
    return complex(real-b.real, imag-b.imag);
  }

  complex operator*(const complex& b) const
  {
    return complex( this->real*b.real - this->imag*b.imag,
	this->real*b.imag + this->imag*b.real);
  }

  complex operator/(const complex& b) const
  {
    return *this * conj(b) / norm(b);
  }

  complex& operator+=(const complex& b)
  {
    this->real += b.real;
    this->imag += b.imag;   
    return *this;
  }

  complex& operator-=(const complex& b)
  {
    this->real -= b.real;
    this->imag -= b.imag;   
    return *this;
  }

  complex& operator*=(const complex& b)
  {
    *this = *this * b;
    return *this;
  }

  complex& operator/=(const complex& b)
  {
    *this *= conj(b) / norm(b);
    return *this;
  }

  /*
   * complex-scalar operations
   */
  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  operator+(const U& b) const
  {
    return complex(this->real+(T)b, this->imag);
  }

  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  operator-(const U& b) const
  {
    return complex(this->real-(T)b, this->imag);
  }

  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  operator*(const U& b) const
  {
    return complex(this->real*(T)b, this->imag*(T)b);
  }

  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  operator/(const U& b) const
  {
    return complex(this->real/(T)b, this->imag/(T)b);
  }

  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  & operator*=(const U& b)
  {
    this->real *= (T)b;
    this->imag *= (T)b;
    return *this;
  }

  template<typename U>
  typename enable_if<is_numeric<U>::val, complex<T> >::type
  & operator/=(const U& b)
  {
    this->real /= (T)b;
    this->imag /= (T)b;
    return *this;
  }
};

  /*
   * scalar-complex operations
   */
  template<typename T, typename U>
  complex<typename enable_if<is_numeric<T>::val, U>::type>
  operator+(const T& a, const complex<U>& b)
  {
    return complex<U>((U)a+b.real, b.imag);
  }

  template<typename T, typename U>
  complex<typename enable_if<is_numeric<T>::val, U>::type>
  operator-(const T& a, const complex<U>& b)
  {
    return complex<U>((U)a-b.real, -b.imag);
  }

  template<typename T, typename U>
  complex<typename enable_if<is_numeric<T>::val, U>::type>
  operator*(const T& a, const complex<U>& b)
  {
    return complex<U>((U)a*b.real, (U)a*b.imag);
  }

  template<typename T, typename U>
  complex<typename enable_if<is_numeric<T>::val, U>::type>
  operator/(const T& a, const complex<U>& b)
  {
    return ((U)a)*conj(b)/norm(b);
  }

  template<typename T>
  T norm(const complex<T>& a)
  {
    return a.real*a.real + a.imag*a.imag;
  }

  template<typename T>
  T abs(const complex<T>& a)
  {
    return ::sqrt(norm(a));
  }

  template<typename T>
  complex<T> conj(const complex<T>& a)
  {
    return complex<T>(a.real, -a.imag);
  }

  template<typename T>
  T arg(const complex<T>& a)
  {
    return atan2(a.imag, a.real);
  }

  template<typename T>
  complex<T> sqrt(const complex<T>& a)
  {
    T mag = ::sqrt(abs(a));
    T phase = arg(a) / 2;
    return complex<T>(mag * cos(phase), mag * sin(phase));
  }

  template<typename T>
  complex<T> exp(const complex<T>& a)
  {
    T mag = ::exp(a.real); 
    T phase = a.imag; 
    return complex<T>(mag * cos(phase), mag * sin(phase));
  }

  template<typename T>
  complex<T> pow(const complex<T>& a, const T& p)
  {
    T mag = ::pow(abs(a), p);
    T phase = arg(a) * p;
    return complex<T>(mag * cos(phase), mag * sin(phase));
  }

  typedef complex<float> float_complex;
  typedef complex<double> double_complex;

} // qcd namespace

#endif


