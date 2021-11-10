//
// myCallback.h
//
#ifndef __MYCALLBACK_H__
#define __MYCALLBACK_H__

#include "stdafx.h"

#ifdef _WIN32

typedef bool (CALLBACK *FPTR)( char *p, int i );
typedef bool (CALLBACK *UPTR)( double dExp, double dInc, double dSum, int i);

#endif


#endif // __MYCALLBACK_H__