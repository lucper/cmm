#ifndef H_DEFS
#define H_DEFS

#ifdef _USE_32
#define INT int32_t
#include "libsais.h"
#endif

#ifdef _USE_64
#define INT int64_t
#include "libsais64.h"
#endif

#define SEP '$'

#endif
