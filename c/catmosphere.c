/**
Test (compiles or not compiles).

clang -std=c99 catmosphere.c -c
dmd catmosphere.o ../libatmosphere_gm.a -L-lblas ../../simple_matrix/libsimple_matrix.a ../../cblas/libcblas.a
*/
#include "catmosphere.h"

int main(int argc, char const *argv[])
{
	// needs to call rt_init!
	
	double alpha1 = catmosphere_nvmm_em_and_coordinate(0,0,NULL,NULL,NULL,NULL,NULL);
	double alpha2 = catmosphere_nvmm_em_and_gradient(0,0,NULL,NULL,NULL,NULL,NULL);
	double alpha3 = catmosphere_nvmm_em(0,0,NULL,NULL,NULL,NULL,NULL);
	return 0;
}
