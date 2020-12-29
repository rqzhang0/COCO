// lattice.h
// Andrei Alexandru
// Jan 2006

#ifndef _LATTICE_H
#define _LATTICE_H

#include "random.h"		// for double_prn 

#include "vector_util.h"	// for vector
#include "su3.h" // for lattice_desc
#include "cpu_vector.h"
#include "site.h" // for site description

/* Definition of globals */
#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

namespace qcd {

struct global_data
{
	lattice_desc* desc;
	int startflag;
	int saveflag;
	double_prn* node_prn;
	site* lattice;
	int* nn[8];
	qcd::su3_field* links;
};
EXTERN global_data g;


// Overlap structure
struct overlap_global_data
{
  int iNumberSmallEigenvalues;
  double *pdSmallEigenvalues;
  qcd::vector *pvecSmallEigenvectors;
  double dInnerError;
  int iMaxInnerIterations;
  int iLoaded;
}; 

#ifdef CONTROL
overlap_global_data overlap = {0, 0, 0, 0, 0, 0};
#else
extern overlap_global_data overlap;
#endif

} // qcd namespace

#endif /* _LATTICE_H */

