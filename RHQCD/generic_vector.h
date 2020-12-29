// generic_vector.h
// Ben Gamari
// August 2009

#ifndef _GENERIC_VECTOR_H
#define _GENERIC_VECTOR_H

#include <iostream>
#include "shared_ptr.h"
#include "layout.h"

namespace qcd {


template <typename TBuffer>
struct generic_vector {
public:
  lattice_desc *desc;
  int length;
  int l_fac;
private:
  shared_ptr<TBuffer> buffer;  
public:
  typename TBuffer::TElement* data;


protected:
  generic_vector(lattice_desc* desc, int length, bool dirty) :
      desc(desc), length(length), 
      buffer(shared_ptr<TBuffer>(new TBuffer(length))),
      data(buffer->data)
  {
    l_fac=buffer->l_fac;
    if (!dirty) clear();
  }

public:
  void clear()
  {
    buffer->clear();
  }


private:
  generic_vector(const generic_vector&);
  generic_vector& operator = (const generic_vector&);
};

} // qcd namespace

#endif

