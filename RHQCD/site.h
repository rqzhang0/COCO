// site.h
// Andrei Alexandru
// August 2009

#ifndef _SITE_H
#define _SITE_H

namespace qcd {

struct site {
  // coordinates of this site 
  unsigned short x,y,z,t;
  
  // is it even or odd? 
  char parity;

  // my index in the array 
  unsigned int index;
};

} // qcd namespace

#endif


