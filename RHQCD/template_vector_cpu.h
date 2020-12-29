// complex_vector_arithmetic_cpu.h
// Andrei Alexandru
// July 2010

#pragma once

#include "comm_low.h"
#include "terminate.h"

// The code here follows the pattern described in the paper
// "Improving Large Vector Operations with C++ Expression Template and ATLAS"
// by L. Plagne and F. Hulsemann.

// This implementation uses the standard algorithms provided by boost library 
// to perform the looping over vector elements.


namespace qcd
{
namespace cpu
{

//--------------------------------------------------------------------------------
//
// Vector class definitions 
//
//--------------------------------------------------------------------------------

// The BaseVector class is used to collect mark all the relevant templates
// that need to be collected by the compiler to for the template expressions.
template <class DERIVED> 
class BaseVectorCpu
{ 
public:
  typedef const DERIVED & CDR; 
  inline CDR getCDR( void ) const { return static_cast<CDR>(*this);}
};


// The Vector class implements a vector of complex numbers
// it relies on the complex_type definition of the VECTOR_TYPE
// and its complex_size.
template <class ELEMENT_TYPE>
class VectorCpu: public BaseVectorCpu< VectorCpu<ELEMENT_TYPE> >
{
public:
  typedef const VectorCpu  StoreType;
  typedef ELEMENT_TYPE ElementType;

  template <class DERIVED> VectorCpu & operator = (const BaseVectorCpu<DERIVED> & right)
  {
    const DERIVED & r=right.getCDR(); 

    // sanity check
    CHECK(veccpu_size == r.size(), "VectorCpu: vectors of unequal size.", ERROR);

#pragma omp parallel for
    for(int i=0; i < (int)veccpu_size; ++i) veccpu_data[i] = r[i];
    return (*this);
  }

  VectorCpu & operator = (const VectorCpu & r) 
  {
    // sanity check
    CHECK(veccpu_size == (size_t)r.size(), "VectorCpu: vectors of unequal size.", ERROR);

#pragma omp parallel for
    for(int i=0; i<(int)veccpu_size; ++i) veccpu_data[i] = r[i];
    return (*this);
  }

  VectorCpu(size_t size, ELEMENT_TYPE* data) : 
  veccpu_size(size), veccpu_data(data) {}


  const ELEMENT_TYPE& operator[](int i) const { return veccpu_data[i]; }
  ELEMENT_TYPE& operator[](int i) { return veccpu_data[i]; }

  inline int size() const { return this->veccpu_size;}
//private:
public:
  size_t        veccpu_size;
  ELEMENT_TYPE* veccpu_data;

private:
  VectorCpu(const VectorCpu&);
};


//--------------------------------------------------------------------------------
//
// Vector operations section
//
//--------------------------------------------------------------------------------

template <class LEFT, class OP, class RIGHT> 
class VectorCpuExpression : public BaseVectorCpu< VectorCpuExpression<LEFT,OP,RIGHT> > 
{
public: 
  typedef const VectorCpuExpression  StoreType; 
  typedef typename RIGHT::ElementType ElementType; 
  
  VectorCpuExpression(const BaseVectorCpu<LEFT> & left, const BaseVectorCpu<RIGHT> & right):
    left_(left.getCDR()), right_(right.getCDR()){}
  
  inline ElementType operator[](size_t i) const { return OP::apply(left_[i], right_[i]); }

  inline size_t size( void ) const { return right_.size();} 

private:
  typename LEFT::StoreType  &left_; 
  typename RIGHT::StoreType &right_;
};


struct Add1
{ 
  template <class T>
  static inline T apply(const T &left, const T &right) { return left+right; } 
};

template <class L, class R > 
VectorCpuExpression< L, Add1, R > operator + (const BaseVectorCpu<L> & left,
    const BaseVectorCpu<R> & right)
{ 
  return VectorCpuExpression<L,Add1,R >(left,right);
}

struct Minus
{ 
  template <class T>
  static inline T apply(const T &left, const T &right) { return left-right; } 
};

template <class L, class R > 
VectorCpuExpression< L, Minus, R > operator - (const BaseVectorCpu<L> & left,
    const BaseVectorCpu<R> & right)
{ 
  return VectorCpuExpression<L,Minus,R >(left,right);
}

struct Times
{ 
  template <class T>
  static inline T apply(const T &left, const T &right) { return left * right; }
};

template <class L, class R > 
VectorCpuExpression< L, Times, R > operator * (const BaseVectorCpu<L> & left,
    const BaseVectorCpu<R> & right)
{ 
  return VectorCpuExpression<L,Times,R >(left,right);
}
//--------------------------------------------------------------------------------
//
// Scalar multiplication section
//
//--------------------------------------------------------------------------------


template <class V, class RT> 
class VectorCpuScalarExpression : public BaseVectorCpu< VectorCpuScalarExpression < V, RT > > 
{ 
public:
  typedef const VectorCpuScalarExpression StoreType; 
  typedef typename V::ElementType ElementType; 

  VectorCpuScalarExpression(const BaseVectorCpu<V> & v, const RT & a):
    a_(a), v_(v.getCDR()) {}
   
  inline ElementType operator[](size_t i) const { return a_*v_[i]; }

  inline size_t size( void ) const { return v_.size();} 

private:
    RT a_; 
    typename V::StoreType &v_;
};

template <class L, class RT> 
VectorCpuScalarExpression<L, RT> operator * (const BaseVectorCpu<L> & v,
    const RT & a)
{ 
  return VectorCpuScalarExpression<L, RT>(v,a);
}


template <class L, class RT> 
VectorCpuScalarExpression<L, RT> operator / (const BaseVectorCpu<L> & v,
    const RT & a)
{ 
  return VectorCpuScalarExpression<L, RT>(v,1.0/a);
}

template <class R, class RT> 
VectorCpuScalarExpression<R, RT> operator * (const RT & a, 
    const BaseVectorCpu<R> & v)
{ 
  return VectorCpuScalarExpression<R, RT>(v,a);
}

template <class L> 
VectorCpuScalarExpression<L, double> operator - (const BaseVectorCpu<L> & v)
{ 
  return VectorCpuScalarExpression<L, double>(v,-1.0);
}

#if 0

//--------------------------------------------------------------------------------
//
// Complex scalar multiplication section
//
//--------------------------------------------------------------------------------


template <class V> 
class VectorCpuComplexScalarExpression : 
  public BaseVectorCpu< VectorCpuComplexScalarExpression < V > > 
{ 
public:
  typedef const VectorCpuComplexScalarExpression StoreType; 
  typedef typename V::ElementType ElementType; 

  VectorCpuComplexScalarExpression(const BaseVectorCpu<V> & v, const complex<ElementType> &a):
    a_(a), v_(v.getCDR()) {}

  inline complex<ElementType>

  inline int size( void ) const { return v_.size();} 

private:
    complex<ElementType> a_; 
    typename V::StoreType & v_;
};

template <class L> 
VectorCpuComplexScalarExpression<L> operator * (const BaseVectorCpu<L> & v,
    const complex<typename L::ElementType> & a)
{ 
  return VectorCpuComplexScalarExpression<L>(v,a);
}

template <class L> 
VectorCpuComplexScalarExpression<L>
operator / (const BaseVectorCpu<L> & v,
    const complex<typename L::ElementType> & a)
{ 
  return VectorCpuComplexScalarExpression<L>(v, 1/a);
}

template <class R> 
VectorCpuComplexScalarExpression<R> operator * (const complex<typename R::ElementType> & a, 
    const BaseVectorCpu<R> & v)
{ 
  return VectorCpuComplexScalarExpression<R>(v,a);
}


//--------------------------------------------------------------------------------
//
// Unary operations section section
//
//--------------------------------------------------------------------------------


template <template<typename> class OP, class RIGHT> 
class VectorCpuUnaryExpression : public BaseVectorCpu< VectorCpuUnaryExpression<OP,RIGHT> > 
{
public: 
  typedef const VectorCpuUnaryExpression StoreType; 
  typedef typename RIGHT::ElementType ElementType; 
  typedef boost::transform_iterator<OP<typename RIGHT::iterator::value_type>, 
	  typename RIGHT::iterator> iterator;

  
  VectorCpuUnaryExpression(const BaseVectorCpu<RIGHT> & right):
    right_(right.getCDR()){}

  iterator begin() const 
  { 
    return iterator(right_.begin(), OP<typename RIGHT::iterator::value_type>()); 
  }
  iterator end() const { return this->begin() + this->size();}
  
  inline int size( void ) const { return right_.size();} 

private:
  typename RIGHT::StoreType & right_;
};

template<typename R>
struct CONJ
{
  typedef R result_type;

  inline result_type operator()(const R &a)
  {
    return boost::make_tuple(boost::get<0>(a),  - boost::get<1>(a));
  }
};

template <class R> 
VectorCpuUnaryExpression<CONJ, R> conj (const BaseVectorCpu<R>& v)
{ 
  return VectorCpuUnaryExpression<CONJ, R>(v);
}

#endif
} // cpu namespace
} //qcd namespace
