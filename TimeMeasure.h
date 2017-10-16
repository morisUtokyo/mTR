#pragma once

#ifdef _WIN32

#include <Windows.h>
typedef struct {
	LARGE_INTEGER start, freq;
} TimeMeasureType;

#else

#include <stdlib.h>
#include <sys/time.h>
typedef struct {
	struct timeval tv;
} TimeMeasureType;

#endif


const TimeMeasureType* TimeMeasureBegin( );
double TimeMeasureEnd( const TimeMeasureType* );