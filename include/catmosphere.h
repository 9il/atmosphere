/**
Extern C bindings.
*/


#ifdef __cplusplus


	#include <cstdlib>


	extern "C" {


#else

	#include <stdlib.h>
	#include <stdbool.h>


#endif


void catmosphere_mix
	(
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		const double* p,
		double* mixture
	);


void catmosphere_gradientDescentIteration
	(
		void (*grad)(size_t, const double* arg, double* result),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* mixture,
		double* pi,
		double* xi,
		double* gamma,
		double* c,
		bool (*tolerance)(double, double)
	);


void catmosphere_coordinateDescentIteration
	(
		void (*grad)(size_t, const double* arg, double* result),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* mixture,
		double* pi,
		double* xi,
		double* gamma,
		bool (*tolerance)(double, double)
	);


void catmosphere_coordinateDescentIterationPartial
	(
		double (*partialDerivative)(double),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* mixture,
		double* pi,
		bool (*tolerance)(double, double)
	);


void catmosphere_coordinateDescentIterationPartial_minusSumOfLogs
	(
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* mixture,
		double* pi,
		bool (*tolerance)(double, double)
	);


#ifdef __cplusplus


	}


#endif
