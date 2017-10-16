#include "TimeMeasure.h"

#ifdef _WIN32

const TimeMeasureType* TimeMeasureBegin( ) {

	TimeMeasureType* p = (TimeMeasureType*) malloc( sizeof( TimeMeasureType ) );

	QueryPerformanceFrequency( &( p->freq ) );
	QueryPerformanceCounter( &( p->start ) );

	return p;

}

double TimeMeasureEnd( const TimeMeasureType* p ) {

	LARGE_INTEGER end;
	QueryPerformanceCounter( &end );

	double ret = (double) ( end.QuadPart - p->start.QuadPart ) / (double) p->freq.QuadPart;

	free( p );

	return ret;

}

#else

const TimeMeasureType* TimeMeasureBegin( ) {

	TimeMeasureType* p = (TimeMeasureType*) malloc( sizeof( TimeMeasureType ) );

	gettimeofday( &( p->tv ), NULL );

	return p;

}

double TimeMeasureEnd( const TimeMeasureType* p ) {

	struct timeval e;
	gettimeofday( &e , NULL);

	double ret = ( e.tv_sec - p->tv.tv_sec ) + ( e.tv_usec - p->tv.tv_usec ) * 1.0E-6;

	free( (void*) p );

	return ret;

}

#endif
