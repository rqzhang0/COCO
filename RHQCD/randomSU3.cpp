// randomSU3.cpp
// Andrei Alexandru
// March 2007

// Generate a random su3matrix distributed uniformly
// on the SU(3) group. The distribution takes into 
// account the invariant measure of the group. We use
// the parametrization and measure of SU(3) described in
// the article "Parametrization of SU(3)" by J.B. Bronzan,
// Phys. Rev. D 38, 1994

#include <math.h>
#include <random.h>
#include <su3.h>
#define PI 3.14159265358979323846

namespace qcd {

void random_su3_matrix(su3_matrix* mat, double_prn* prn_pt)
{
  // We need 8 random numbers from an uniform distribution.
  double z1, z2, z3, cth1, sth1, cth2, sth2, cth3, sth3;
  double f1, f2, f3, f4, f5;

  z1 = prn_pt->rand();
  z2 = prn_pt->rand();
  z3 = prn_pt->rand();
  cth1 = pow(z1, 0.25); // z1 = cos th1 ^ 4
  sth1 = sqrt(1-cth1*cth1);
  cth2 = sqrt(z2);  // z2 = cos th2 ^ 2
  sth2 = sqrt(1-z2);
  cth3 = sqrt(z3);  // z3 = cos th3 ^ 2
  sth3 = sqrt(1-z3);

  f1 = prn_pt->rand()*2*PI;
  f2 = prn_pt->rand()*2*PI;
  f3 = prn_pt->rand()*2*PI;
  f4 = prn_pt->rand()*2*PI;
  f5 = prn_pt->rand()*2*PI;

  // u11 = cos th1 cos th2 exp(i f1) 
  mat->e[0][0].real = cth1*cth2*cos(f1);
  mat->e[0][0].imag = cth1*cth2*sin(f1);
  // u12 = sin th1 exp(i f3)
  mat->e[0][1].real = sth1*cos(f3);
  mat->e[0][1].imag = sth1*sin(f3);
  // u13 = cos th1 sin th2 exp(i f4)
  mat->e[0][2].real = cth1*sth2*cos(f4);
  mat->e[0][2].imag = cth1*sth2*sin(f4);

  // u21 = sin th2 sin th3 exp(-i(f4+f5)) - sin th1 cos th2 cos th3 exp(i(f1+f2-f3))
  mat->e[1][0].real = sth2*sth3*cos(f4+f5)-sth1*cth2*cth3*cos(f1+f2-f3);
  mat->e[1][0].imag = -sth2*sth3*sin(f4+f5)-sth1*cth2*cth3*sin(f1+f2-f3);
  // u22 = cos th1 cos th3 exp(i f2)
  mat->e[1][1].real = cth1*cth3*cos(f2);
  mat->e[1][1].imag = cth1*cth3*sin(f2);
  // u23 = -cos th2 sin th3 exp(-i(f1+f5))-sin th1 sin th2 cos th3 exp(i(f2-f3+f4))
  mat->e[1][2].real = -cth2*sth3*cos(f1+f5)-sth1*sth2*cth3*cos(f2-f3+f4);
  mat->e[1][2].imag = cth2*sth3*sin(f1+f5)-sth1*sth2*cth3*sin(f2-f3+f4);

  // u31 = -sin th1 cos th2 sin th3 exp(i(f1-f3+f5))-sin th2 cos th3 exp(-i(f2+f4))
  mat->e[2][0].real = -sth1*cth2*sth3*cos(f1-f3+f5)-sth2*cth3*cos(f2+f4);
  mat->e[2][0].imag = -sth1*cth2*sth3*sin(f1-f3+f5)+sth2*cth3*sin(f2+f4);
  // u32 = cos th1 sin th3 exp(i f5)
  mat->e[2][1].real = cth1*sth3*cos(f5);
  mat->e[2][1].imag = cth1*sth3*sin(f5);
  // u33 = cos th2 cos th3 exp(-i(f1+f2)) - sin th1 sin th2 sin th3 exp(i(-f3+f4+f5))
  mat->e[2][2].real = cth2*cth3*cos(f1+f2)-sth1*sth2*sth3*cos(-f3+f4+f5);
  mat->e[2][2].imag = -cth2*cth3*sin(f1+f2)-sth1*sth2*sth3*sin(-f3+f4+f5);

}

} // qcd namespace
