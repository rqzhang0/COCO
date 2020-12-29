#include <stdlib.h>

#ifdef MACOSX
void * memalign(size_t a, size_t b)
{
  return malloc(b);
}
#endif
