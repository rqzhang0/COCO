// vector_util.h
// Andrei Alexandru
// Nov 2003

// Includes declarations for vector operations
// This is what we call big oprations: they take a lot of time and
// are designed to run on multiple nodes.

#ifndef _VECTOR_UTIL_H
#define _VECTOR_UTIL_H

#include "cpu_vector.h"
#include "comm_low.h"
#include <stdio.h>
#include <math.h>

namespace qcd {

//------------------------------------------------------------------------------------
// Packing functions
//------------------------------------------------------------------------------------

void pack_spinor(vector& vec, double* packed_vec);
void unpack_spinor(double* packed_vec, vector& vec);
void pack_links(su3_field& lnk, double* packed_links);

void pack_spinor_eo(vector& vec, double* packed_vec, enum parity parity);
void unpack_spinor_eo(double* packed_vec, vector& vec, enum parity parity);
void pack_links_eo(su3_field& lnk, double* packed_links, enum parity parity);


//------------------------------------------------------------------------------------
// Utility functions
//------------------------------------------------------------------------------------

void gamma5_multiply(const vector& vecSrc, vector& vecDest);
void gamma5_projection(const vector& vecSrc, vector& vecDest, int chirality);
void vector_matrix_multiply(vector* pvecV[], double_complex* pcQ, int iDim, int iCol);
void vector_matrix_multiply_indexed(vector* pvecV[], double_complex* pcQ, int iDim, int iCol, int *piIndex);

double norm2(const vector& vecV);
double_complex cscalar_product2(const vector& vector1, const vector& vector2);
double scalar_product(const vector& vector1, const vector& vector2);

void vec1_plus_vec2(const vector& vec1, const vector& vec2, const vector& result);
void complex_times_vector(const double_complex& scalar, const vector& vec, const vector& result);
void scalar_times_vector(const double scalar, const vector& vec1, const vector& result);
void vec1_plus_complex_times_vec2(const vector& vec1, const double_complex& scalar, const vector& vec2, const vector& result);
void vec1_plus_scalar_times_vec2(const vector& vec1, const double scalar, const vector& vec2, const vector& result);
void cm1_times_vec1_plus_cm2_times_vec2(double_complex& sc1, const vector& vec1, double_complex& sc2, const vector& vec2, vector& result);
void sc1_times_vec1_plus_sc2_times_vec2(double sc1, const vector& vec1, double sc2, const vector& vec2, vector& result);
void copy_vec1_to_vec2(const vector& pvSrc, vector& pvDest);

void clear_vector(vector& vecVector);

void chiral_expand(const chiralvector& src, vector& dest, int chirality);
void gamma5_projection(const vector& src, chiralvector& dest, int chirality);


template<typename EXP>
double norm(const typename qcd::cpu::BaseVectorCpu<EXP>& vecV)
{
  double dRes = 0.0;
  const EXP & vr = vecV.getCDR();
#pragma omp parallel for reduction(+: dRes)
  for(int i=0; i<(int)vr.size(); ++i) dRes += norm(vr[i]); 

  global_sum(dRes);
  return ::sqrt(dRes);
}


template<typename EXP>
double_complex cscalar_product(
    const typename qcd::cpu::BaseVectorCpu<EXP>& vec1,
    const typename qcd::cpu::BaseVectorCpu<EXP>& vec2)
{
  double resre = 0.0, resim = 0.0;
  const EXP & vr1 = vec1.getCDR();
  const EXP & vr2 = vec2.getCDR();
#pragma omp parallel for reduction(+: resre) reduction(+: resim)
  for(int i=0; i<vr1.size(); ++i) 
  {
     double_complex tmp = conj(vr1[i])*vr2[i]; 
     resre += tmp.real;
     resim += tmp.imag;
  }
  double_complex res(resre, resim);
  global_sum(res);
  return res;
}

// Multiply a matrix V which is a collection of k(iDim) vectors on the columns
// with a k x k matrix Q. It only does that for the first iCol vectors.
// Thus V->VQ where only the first iCol are changed, the rest stay the same.
template<typename RT, typename VEC>
void vector_matrix_multiply(VEC* pvecV[], 
    RT* pcQ, int iDim, int iCol)
{
  int len = pvecV[0]->size();
  RT* Q = new RT[iDim*iCol];
  for(int i=0; i<iCol; ++i)
  for(int j=0; j<iDim; ++j) Q[i*iDim+j] = pcQ[j*iDim+i];

#pragma omp parallel shared(Q,iDim,iCol,pvecV,len)
{
  RT *pwvTemp = new RT[iDim];
#pragma omp for
  for(int i=0; i<len; i++)
  {
    for(int j = 0; j < iDim; j++) 
      pwvTemp[j] = (*pvecV[j])[i];

    for(int j = 0; j < iCol; j++)
    {
      RT tmp = 0;
      for(int k = 0; k < iDim; k++)
      {
	//tmp += pwvTemp[k]*pcQ[k*iDim+j];
	tmp += pwvTemp[k]*Q[j*iDim+k];
      }
      (*pvecV[j])[i] = tmp;
    }
  }
  delete [] pwvTemp;
}

  delete [] Q;

}

// Does the same thing as the function above only that it does it for the columns k columns in
// the index vector piIndex.
template<typename RT, typename VEC>
void vector_matrix_multiply_indexed(VEC* pvecV[], 
    RT* pcQ, int iDim, int iCol, int *piIndex)
{
  int len = pvecV[0]->size();
  RT* Q = new RT[iDim*iCol];
  for(int i=0; i<iCol; ++i)
  for(int j=0; j<iDim; ++j) Q[i*iDim+j] = pcQ[j*iDim+piIndex[i]];

#pragma omp parallel shared(pcQ,iDim,iCol,pvecV,piIndex,len)
{
  RT *pwvTemp = new RT[iDim];
#pragma omp for
  for(int i=0; i<len; i++)
  {
    for(int j = 0; j < iDim; j++) 
      pwvTemp[j] = (*pvecV[j])[i];

    for(int j = 0; j < iCol; j++)
    {
      int jindx = piIndex[j];
      RT wvRes = 0.0;
      for(int k = 0; k < iDim; k++)
      {
	//wvRes += pcQ[k*iDim+jindx]*pwvTemp[k];
	wvRes += Q[j*iDim+k]*pwvTemp[k];
      }
      (*pvecV[jindx])[i] = wvRes;
    }
  }
  delete [] pwvTemp;
}
  
  delete [] Q;

}

} // qcd namespace

#endif

