#pragma once

#ifndef MACOSX
#include <malloc.h>
#else

#include <stdlib.h>

void* memalign(size_t a, size_t b);
#endif
