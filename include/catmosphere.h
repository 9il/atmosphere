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


void catmosphere_gradientDescentIteration
	(
		void (*grad)(size_t, const double* arg, double* result),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* gamma,
		double* c,
		bool (*tolerance)(double, double)
	);


void catmosphere_simpleGradientDescentIteration
	(
		double (*simpleGrad)(double),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
		double* c,
		bool (*tolerance)(double, double)
	);


void catmosphere_minusSumOfLogs_gradientDescentIteration
	(
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		double* xi,
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
		double* chi,
		double* pi,
		double* xi,
		double* gamma,
		bool (*tolerance)(double, double)
	);


void catmosphere_simpleCoordinateDescentIteration
	(
		double (*simpleGrad)(double),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		bool (*tolerance)(double, double)
	);


void catmosphere_minusSumOfLogs_coordinateDescentIteration
	(
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* pi,
		bool (*tolerance)(double, double)
	);


void catmosphere_EMLikeDescentIteration
	(
		void (*grad)(size_t, const double* arg, double* result),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* chi,
		double* xi
	);


void catmosphere_simpleEMLikeDescentIteration
	(
		double (*simpleGrad)(double),
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* xi,
		double* c
	);


void catmosphere_minusSumOfLogs_EMLikeDescentIteration
	(
		const double* WTptr,
		size_t k,
		size_t n,
		size_t WTshift,
		double* p,
		double* xi,
		double* c
	);

#ifdef __cplusplus
}
#endif
