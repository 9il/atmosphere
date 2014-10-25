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


double catmosphere_nvmm_em_and_coordinate(
		size_t k, 
		size_t n,
		double* weights,
		const double* grid, 
		const double* sample,
		bool (*tolerance) (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			const double* weightsPrev, 
			const double* weights
		),
		bool (*findRootTolerance)(double a, double b)
	);


double catmosphere_nvmm_em_and_gradient(
		size_t k, 
		size_t n,
		double* weights,
		const double* grid, 
		const double* sample,
		bool (*tolerance) (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			const double* weightsPrev, 
			const double* weights
		),
		bool (*findRootTolerance)(double a, double b)
	);


double catmosphere_nvmm_em(
		size_t k, 
		size_t n,
		double* weights,
		const double* grid, 
		const double* sample,
		bool (*tolerance) (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			const double* weightsPrev, 
			const double* weights
		),
		bool (*findRootTolerance)(double a, double b)
	);


#ifdef __cplusplus
}
#endif
