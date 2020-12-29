// su3.cpp
// Ben Gamari
// August 2009

#include "su3.h"

namespace qcd {

template<typename RT>
generic_su3_vector<RT> operator*(const generic_su3_matrix<RT>& a, const generic_su3_vector<RT>& b) 
{
  generic_su3_vector<RT> res;
  for (int r=0; r<3; r++)
  for (int c=0; c<3; c++)
    res.c[r] += a.e[r][c] * b.c[c];
  return res;
}

template generic_su3_vector<double> operator*(const generic_su3_matrix<double>& a, 
    const generic_su3_vector<double>& b);
template generic_su3_vector<float> operator*(const generic_su3_matrix<float>& a, 
    const generic_su3_vector<float>& b);

template<typename RT>
generic_su3_matrix<RT> adj(const generic_su3_matrix<RT>& a)
{
  generic_su3_matrix<RT> res;
  for (int r=0; r<3; r++)
  for (int c=0; c<3; c++)
    res.e[r][c] = conj(a.e[c][r]);
  return res;
}

template
generic_su3_matrix<double> adj(const generic_su3_matrix<double>& a);
template
generic_su3_matrix<float> adj(const generic_su3_matrix<float>& a);

template<typename RT>
generic_su3_matrix<RT> operator*(const generic_su3_matrix<RT>& a, const generic_su3_matrix<RT>& b)
{
  generic_su3_matrix<RT> res = a;
  res *= b;
  return res;
}

template generic_su3_matrix<double> operator*(const generic_su3_matrix<double>& a, 
    const generic_su3_matrix<double>& b);
template generic_su3_matrix<float> operator*(const generic_su3_matrix<float>& a, 
    const generic_su3_matrix<float>& b);

} // qcd namespace
