// cpu_vector.h
// Ben Gamari
// August 2009

#ifndef _CPU_VECTOR_H
#define _CPU_VECTOR_H

#include "malloc-stub.h"
#include <stdlib.h>
#include <string.h>

#include "su3.h"
#include "terminate.h"
#include "generic_vector.h"
#include "template_vector_cpu.h"
#include "complex.h"

namespace qcd {

template <typename T, int length_factor=1>
struct cpu_buffer {
  typedef T TElement;
  typedef T* iterator;
  int l_fac;

  int length;
  TElement* data;

  cpu_buffer(int length) : length(length),l_fac(length_factor) {
    this->data = (TElement*) memalign(16, length_factor*length*sizeof(TElement));
    if (this->data == NULL) {
      terminate_error(ERROR_MEMORY_ALLOCATION);
      exit(ERROR_MEMORY_ALLOCATION);
    }
  }

  ~cpu_buffer() 
  { 
    free(this->data); 
  }

  void clear()
  {
    //memset(data, 0, length_factor*length*sizeof(TElement));
#pragma omp parallel for
    for(int i=0; i<length_factor*length; ++i) data[i] = T();
  }
};

template<typename RT>
struct generic_wilson_field :  generic_vector<cpu_buffer<generic_wilson_vector<RT> > >, qcd::cpu::VectorCpu<complex<RT> >
{
  generic_wilson_field(lattice_desc* desc, bool dirty=false) :
    generic_vector<cpu_buffer<generic_wilson_vector<RT> > >(desc, desc->sites_on_node, dirty),
    qcd::cpu::VectorCpu<complex<RT> >(12*desc->sites_on_node, reinterpret_cast<complex<RT>*>(this->data)) 
    {}
      
  using  qcd::cpu::VectorCpu<complex<RT> >::operator = ; 
  // this assignement operator masks the one from generic_vector
  generic_wilson_field<RT>& operator = (const generic_wilson_field<RT>& v) { qcd::cpu::VectorCpu<complex<RT> >::operator = (v); return *this;}
private:
  generic_wilson_field<RT>(const generic_wilson_field<RT>&);
};

typedef generic_wilson_field<double> vector;


struct chiralvector :  generic_vector<cpu_buffer<su3_vector, 2> >, qcd::cpu::VectorCpu<double_complex>
{
  chiralvector(lattice_desc* desc, bool dirty=false) :
    generic_vector<cpu_buffer<su3_vector, 2> >(desc, desc->sites_on_node, dirty),
    qcd::cpu::VectorCpu<double_complex>(6*desc->sites_on_node, reinterpret_cast<double_complex*>(this->data)) 
    {}
      
  using  qcd::cpu::VectorCpu<double_complex>::operator = ; 
  // this assignement operator masks the one from generic_vector
  chiralvector& operator = (const chiralvector& v) { qcd::cpu::VectorCpu<double_complex>::operator = (v); return *this;}
private:
  chiralvector(const chiralvector&);
};

/*
template<typename RT>
struct generic_su3_field : generic_vector<cpu_buffer<generic_su3_matrix<RT>, 4> >{
  generic_su3_field<RT>(lattice_desc* desc, bool dirty=false) :
    generic_vector<cpu_buffer<generic_su3_matrix<RT>, 4> >(desc, desc->sites_on_node, dirty) { }
}
*/
template<typename RT>
struct generic_su3_field : generic_vector<cpu_buffer<generic_su3_matrix<RT>, 4> >{
  generic_su3_field<RT>(lattice_desc* desc, bool dirty=false) :
    generic_vector<cpu_buffer<generic_su3_matrix<RT>, 4> >(desc, desc->sites_on_node, dirty) { }

  generic_su3_field<RT>& operator = (const generic_su3_field<RT> &src)
  {
#pragma omp parallel for
    for(int i=0; i<4*this->desc->sites_on_node; ++i) this->data[i] = src.data[i];
    return *this;
  }
};

typedef generic_su3_field<double> su3_field;
typedef generic_su3_field<float> float_su3_field;

struct su3_matrix_field : generic_vector<cpu_buffer<su3_matrix, 1> > {
  su3_matrix_field(lattice_desc* desc, bool dirty=false) :
    generic_vector<cpu_buffer<su3_matrix, 1> >(desc, desc->sites_on_node, dirty) { }
  su3_matrix_field& operator = (const su3_matrix_field &src)
  {
#pragma omp parallel for
    for(int i=0; i<desc->sites_on_node; ++i) data[i] = src.data[i];
    return *this;
  }
};

} // qcd namespace

#endif

