#ifndef _DEFINE_H_
#define _DEFINE_H_


#include <math.h>

#define DOUBLEPRECISION

#if defined(DOUBLEPRECISION)
	typedef double comp;
	typedef double real;
	typedef float real2;
#elif defined(FlOATPRECISION)
	typedef float comp;
	typedef float real;
#elif defined(DOUBLECOMPLEX)
	typedef double complex comp;
	typedef double complex real;
	typedef float real2;
#endif

const real PI  = 3.1415926535897932384626433832795;
const real EP0 = 8.85418782E-12;
const real MU0 = 1.25663706E-06;
const real C0  = 2.99792458E+08;
const real Z0  = 376.730313462;
const real PLANKS = 4.13566553853599E-15;
const real E0 = 1.60217662E-19;
const real ME = 9.10938356E-31;


#endif // _CONSTANTS_H_
