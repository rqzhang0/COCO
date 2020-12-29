#ifndef _RANDOMME_H
#define _RANDOMME_H

#include "su3.h"

namespace qcd {

/* random number structures */

struct double_prn {
  /* We assume long is at least 32 bits */
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;

  double_prn(int seed=123, int index=0){ set_seed(seed, index); }

  void set_seed(int seed, int index=0);

  /* Generic random number generator returning a uniformly distributed
     random value on [0,1] */
  float rand();
};

//void random_su3_matrix(su3_matrix* m, double_prn* rnd);

} // qcd namespace

#endif /* _RANDOMME_H */


