module atmosphere.catmosphere;

import std.functional;

import atmosphere.parametrized.nvmm;

extern(C):


double catmosphere_nvmm(NVMM : NormalVarianceMeanMixture!double)(
		size_t k, 
		size_t n,
		double* weights,
		in double* grid, 
		in double* sample,
		in bool function (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			in double* weightsPrev, 
			in double* weights,
		) 
		tolerance,
		in bool function(double a, double b) @nogc nothrow findRootTolerance
	)
{
	auto op = new NVMM(grid[0..k], n);
	op.weights = weights[0..k];
	op.sample = sample[0..n];
	op.optimize(
		(	
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			in double[] weightsPrev, 
			in double[] weights,
		) => tolerance(
			alphaPrev,
			alpha,
			log2LikelihoodValuePrev,
			log2LikelihoodValue,
			weightsPrev.ptr,
			weights.ptr,
		), 
		(a, b) => findRootTolerance(a, b)
		);
	weights[0..k] = op.weights[];
	return op.alpha;
}


double catmosphere_nvmm_em_and_coordinate(
		size_t k, 
		size_t n,
		double* weights,
		const double* grid, 
		const double* sample,
		bool function (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			const double* weightsPrev, 
			const double* weights,
		) 
		tolerance,
		bool function(double a, double b) @nogc nothrow findRootTolerance
	)
{
	return catmosphere_nvmm!(NormalVarianceMeanMixtureEMAndCoordinate!double)
	(
		k,
		n,
		weights,
		grid,
		sample,
		tolerance,
		findRootTolerance
	);
}


double catmosphere_nvmm_em_and_gradient(
		size_t k, 
		size_t n,
		double* weights,
		const double* grid, 
		const double* sample,
		bool function (
			double alphaPrev, 
			double alpha, 
			double log2LikelihoodValuePrev, 
			double log2LikelihoodValue, 
			const double* weightsPrev, 
			const double* weights,
		) 
		tolerance,
		bool function(double a, double b) @nogc nothrow findRootTolerance
	)
{
	return catmosphere_nvmm!(NormalVarianceMeanMixtureEMAndGradient!double)
	(
		k,
		n,
		weights,
		grid,
		sample,
		tolerance,
		findRootTolerance
	);
}
