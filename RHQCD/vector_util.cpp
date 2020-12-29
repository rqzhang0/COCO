// vector_util.c
// Andrei Alexandru
// Nov 2003

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc-stub.h>
#include "su3.h"
#include "vector_util.h"
#include "terminate.h"
#include "lattice.h"
#include "comm/comm_low.h"

namespace qcd {

//------------------------------------------------------------------------------------
// Utility functions
//------------------------------------------------------------------------------------


void gamma5_multiply(const vector& src, vector& dest)
{
  int iIndex, iColor;
  wilson_vector* vecSrc = src.data;
  wilson_vector* vecDest = dest.data;

  for(iIndex = 0; iIndex < src.length; iIndex++)
  for(iColor = 0; iColor < 3; iColor++)
  {
    vecDest[iIndex].d[0].c[iColor] =  vecSrc[iIndex].d[0].c[iColor];
    vecDest[iIndex].d[1].c[iColor] =  vecSrc[iIndex].d[1].c[iColor];
    vecDest[iIndex].d[2].c[iColor] = -vecSrc[iIndex].d[2].c[iColor];
    vecDest[iIndex].d[3].c[iColor] = -vecSrc[iIndex].d[3].c[iColor];
  }
}

void gamma5_projection(const vector& src, vector& dest, int chirality)
{
  int iIndex, iColor;
  wilson_vector* vecSrc = src.data;
  wilson_vector* vecDest = dest.data;


  if(chirality == 1)
  for(iIndex = 0; iIndex < src.length; iIndex++)
  for(iColor = 0; iColor < 3; iColor++)
  {
    vecDest[iIndex].d[0].c[iColor] = vecSrc[iIndex].d[0].c[iColor];
    vecDest[iIndex].d[1].c[iColor] = vecSrc[iIndex].d[1].c[iColor];
    vecDest[iIndex].d[2].c[iColor] = 0.0;
    vecDest[iIndex].d[3].c[iColor] = 0.0;
  }
  else if(chirality == -1)
  for(iIndex = 0; iIndex < src.length; iIndex++)
  for(iColor = 0; iColor < 3; iColor++)
  {
    vecDest[iIndex].d[0].c[iColor] = 0.0;
    vecDest[iIndex].d[1].c[iColor] = 0.0;
    vecDest[iIndex].d[2].c[iColor] = vecSrc[iIndex].d[2].c[iColor];
    vecDest[iIndex].d[3].c[iColor] = vecSrc[iIndex].d[3].c[iColor];
  }
  else exit(terminate_error(ERROR_IN_ARGUMENTS));
}


void gamma5_projection(const vector& src, chiralvector& dest, int chirality)
{
  int iIndex;
  wilson_vector* vecSrc = src.data;
  su3_vector* vecDest = dest.data;


  if(chirality == 1)
  {
#pragma omp parallel for
    for(iIndex = 0; iIndex < src.length; iIndex++)
    {
      vecDest[iIndex*2+0] = vecSrc[iIndex].d[0];
      vecDest[iIndex*2+1] = vecSrc[iIndex].d[1];
    }
  }
  else if(chirality == -1)
  {
#pragma omp parallel for
    for(iIndex = 0; iIndex < src.length; iIndex++)
    {
      vecDest[iIndex*2+0] = vecSrc[iIndex].d[2];
      vecDest[iIndex*2+1] = vecSrc[iIndex].d[3];
    }
  }
  else exit(terminate_error(ERROR_IN_ARGUMENTS));
}

void chiral_expand(const chiralvector& src, vector& dest, int chirality)
{
  int idx, oidx;
  if(chirality == 1) { idx = 0; oidx = 2; }
  else if(chirality == -1) { idx = 2; oidx = 0; }
  else exit(terminate_error(ERROR_IN_ARGUMENTS));

  su3_vector zero;

#pragma omp parallel for
  for(int i=0; i< src.length; ++i)
  {
    dest.data[i].d[idx+0] = src.data[2*i+0];
    dest.data[i].d[idx+1] = src.data[2*i+1];
    dest.data[i].d[oidx+0] = zero;
    dest.data[i].d[oidx+1] = zero;
  }
}

#if 0 // use templated definitions which are faster

static void c_scalar_mult_wvec(const wilson_vector& src, const double_complex& phase, wilson_vector& res)
{
  for (int i=0;i<4;i++)
  for (int j=0;j<3;j++)
    res.d[i].c[j] = src.d[i].c[j] * phase;
}

static void clear_wvec(wilson_vector* a)
{
  memset(a, 0, sizeof(wilson_vector));
}

static void add_wvec(const wilson_vector& src1, const wilson_vector& src2, wilson_vector& res)
{
  for (int i=0; i<4; i++)
  for (int j=0; j<3; j++) {
    // These braces heavily influence performance on gcc-4.3.3
    res.d[i].c[j] = src1.d[i].c[j] + src2.d[i].c[j];
  }
}

// Multiply a matrix V which is a collection of k(iDim) vectors on the columns
// with a k x k matrix Q. It only does that for the first iCol vectors.
// Thus V->VQ where only the first iCol are changed, the rest stay the same.
void vector_matrix_multiply(vector* pvecV[], double_complex* pcQ, int iDim, int iCol)
{
  int i, j, k;
  wilson_vector *pwvTemp, wvTemp;
  int len = pvecV[0]->length;

#pragma omp parallel private(i,j,k,pwvTemp,wvTemp) shared(len,iDim,iCol,pcQ,pvecV)
{
  pwvTemp = new wilson_vector[iDim];

#pragma omp for
  for(i=0; i<len; i++)
  {
    for(j = 0; j < iDim; j++) 
      pwvTemp[j] = pvecV[j]->data[i];

    for(j = 0; j < iCol; j++)
    {
      clear_wvec(&pvecV[j]->data[i]);
      for(k = 0; k < iDim; k++)
      {
	c_scalar_mult_wvec(pwvTemp[k], pcQ[k*iDim+j], wvTemp);
        add_wvec(pvecV[j]->data[i], wvTemp, pvecV[j]->data[i]);
      }
    }
  }

  delete [] pwvTemp;
}

}

// Does the same thing as the function above only that it does it for the columns k columns in
// the index vector piIndex.
void vector_matrix_multiply_indexed(vector* pvecV[], double_complex* pcQ, int iDim, int iCol, int *piIndex)
{
  int len = pvecV[0]->length;

#pragma omp parallel shared(len,iDim,iCol,pcQ,pvecV,piIndex)
{
  wilson_vector* pwvTemp = new wilson_vector[iDim];
  wilson_vector wvRes, wvTemp;

#pragma omp for
  for(int i=0; i<len; i++)
  {
    for(int j = 0; j < iDim; j++) 
      pwvTemp[j] = pvecV[j]->data[i];
 
    for(int j = 0; j < iCol; j++)
    {
      int jindx = piIndex[j];
      clear_wvec(&wvRes);
      for(int k = 0; k < iDim; k++)
      {
	c_scalar_mult_wvec(pwvTemp[k], pcQ[k*iDim+jindx], wvTemp);
	add_wvec(wvRes, wvTemp, wvRes);
      }
      pvecV[jindx]->data[i] = wvRes;
    }
  }

  delete [] pwvTemp;
}

}

#else

// Multiply a matrix V which is a collection of k(iDim) vectors on the columns
// with a k x k matrix Q. It only does that for the first iCol vectors.
// Thus V->VQ where only the first iCol are changed, the rest stay the same.
void vector_matrix_multiply(vector* pvecV[], double_complex* pcQ, int iDim, int iCol)
{
  vector_matrix_multiply<double_complex>(pvecV, pcQ, iDim, iCol);
}

// Does the same thing as the function above only that it does it for the columns k columns in
// the index vector piIndex.
void vector_matrix_multiply_indexed(vector* pvecV[], double_complex* pcQ, int iDim, int iCol, int *piIndex)
{
  vector_matrix_multiply_indexed<double_complex>(pvecV, pcQ, iDim, iCol, piIndex);
}

#endif

double norm2(const vector& vecV)
{
  double dRes;

  dRes = scalar_product(vecV, vecV);
  return sqrt(dRes);
}



double_complex cscalar_product2(const vector& vec1, const vector& vec2)
{
  register int i;
  register double real = 0.0;
  register double imag = 0.0;
  register double_complex* v1 = (double_complex*) vec1.data;
  register double_complex* v2 = (double_complex*) vec2.data;
  double_complex cRes;
  int len = vec1.length;


  for(i=len * 12; i>0; --i, ++v1, ++v2) 
  {
    real += (*v1).real * (*v2).real;
    real += (*v1).imag * (*v2).imag;
    imag += (*v1).real * (*v2).imag;
    imag -= (*v1).imag * (*v2).real;
  }
  global_sum(real);
  global_sum(imag);
  return double_complex(real, imag);
}

// This is not a generic scalar product.
// It will return only the real part of it.
// The fact is that the cg inverter only needs 
// scalar products that are expected to be real
double scalar_product(const vector& vector1, const vector& vector2)
{
  register int i;
  register double real = 0.0;
  register double_complex* v1 = (double_complex*) vector1.data;
  register double_complex* v2 = (double_complex*) vector2.data;
  int len = vector1.length;

  for(i=len*12; i>0; --i, v1++, v2++) 
  {
		real += (*v1).real * (*v2).real;
		real += (*v1).imag * (*v2).imag;
  }
  global_sum(real);
  return real;
}


void vec1_plus_vec2(const vector& vec1, const vector& vec2, const vector& result)
{
  register int i = 24*vec1.length;
  register double* v1 = (double *) vec1.data;
  register double* v2 = (double *) vec2.data;
  register double* vr = (double *) result.data;

  for(; i>0; --i, ++v1, ++v2, ++vr)
    (*vr) = (*v1) + (*v2);
}

void complex_times_vector(const double_complex& scalar, const vector& vec, const vector& result)
{
  register int i=12*vec.length;
  register double_complex* v1 = (double_complex *) vec.data;
  register double_complex* vr = (double_complex *) result.data;

  for(; i>0; --i, ++v1, ++vr)
    (*vr)  = scalar * (*v1);
}

void scalar_times_vector(const double scalar, const vector& vec1, const vector& result)
{
  register int i = 24 * vec1.length;
  register double* v1 = (double *) vec1.data;
  register double* vr = (double *) result.data;

  for(; i>0; --i, ++v1, ++vr)
    (*vr) = scalar * (*v1);
}

void vec1_plus_complex_times_vec2(const vector& vec1, const double_complex& scalar, const vector& vec2, const vector& result)
{
  register int i=12 * vec1.length;
  register double_complex* v1 = (double_complex *) vec1.data;
  register double_complex* v2 = (double_complex *) vec2.data;
  register double_complex* vr = (double_complex *) result.data;

  for(; i>0; --i, ++v1, ++v2, ++vr) {
    (*vr)  = scalar * (*v2) + (*v1);
  }
}

void vec1_plus_scalar_times_vec2(const vector& vec1, const double scalar, const vector& vec2, const vector& result)
{
  register int i = 24*vec1.length;
  register double* v1 = (double *) vec1.data;
  register double* v2 = (double *) vec2.data;
  register double* vr = (double *) result.data;
  
  for(; i>0; --i, ++v1, ++v2, ++vr)
    (*vr) = scalar * (*v2) + (*v1); 
}

void cm1_times_vec1_plus_cm2_times_vec2(double_complex& sc1, const vector& vec1, double_complex& sc2, const vector& vec2, vector& result)
{
  register int i = 12*vec1.length;
  register double_complex* v1 = (double_complex *) vec1.data;
  register double_complex* v2 = (double_complex *) vec2.data;
  register double_complex* vr = (double_complex *) result.data;

  for(; i>0; --i, ++v1, ++v2, ++vr)
		(*vr) = sc1 * (*v1) + sc2 * (*v2);
}

void sc1_times_vec1_plus_sc2_times_vec2(double sc1, const vector& vec1, double sc2, const vector& vec2, vector& result)
{
  register int i = 24*vec1.length;
  register double* v1 = (double *) vec1.data;
  register double* v2 = (double *) vec2.data;
  register double* vr = (double *) result.data;

  for(; i>0; --i, ++v1, ++v2, ++vr)
    (*vr) = sc1 * (*v1) + sc2 * (*v2);
}

void copy_vec1_to_vec2(const vector& pvSrc, vector& pvDest)
{
  memcpy(pvDest.data, pvSrc.data, 24*pvSrc.length*sizeof(double));
}

void clear_vector(vector& vecVector)
{
  int i = 24 * vecVector.length;
  double* v = (double *) vecVector.data;

  for(; i>0; --i, ++v) (*v) = 0.0;
}


} // qcd namespace
