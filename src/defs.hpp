#ifndef H_DEFS
#define H_DEFS

#ifdef _USE_32
#define INT int32_t
#define UINT uint32_t
#define WSIZE 32
#include "libsais.h"
#endif

#ifdef _USE_64
#define INT int64_t
#define UINT uint64_t
#define WSIZE 64
#include "libsais64.h"
#endif

#define SEP '$'

#endif
