module atmosphere.mixtureoptimizer;

import std.traits : isFloatingPoint;

import atmosphere.internal;

/**
Base abstract class for all mixture distribution optimization algorithms.
*/
abstract class MixtureOptimizer(T)
	if(isFloatingPoint!T)
{

	/**
	Perform one iteration of optimization.
	*/
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

	inout(Matrix!T) _componentsT() inout @property;
	inout(T)[] _mixture() inout @property;
	inout(T)[] _distribution() inout @property;
	void _distribution(T[] _distribution_) @property;


	void updateMixture()
	{
		mix(cast(Matrix!(const T))_componentsT, _distribution, _mixture);
	}

//final:


	/**
	Returns:
		Const internal components representation.
	*/
	TransposedMatrix!(const(T)) components() const
	{
		return cast(typeof(return))_componentsT.transposed;
	}



	/**
	Returns:
		Const slice to the internal mixture representation.
	Example:
	-------------
	double objectiveFunction(in double[])
	{
	
	}

	//save slice
	auto mixture = optimizer.mixture;

	auto value0 = objectiveFunction(mixture);
	optimizer.eval;
	auto value1 = objectiveFunction(mixture);


	//use $(D .dup) or copy to save current mixture

	//1: .dup
	auto mixtureSave1 = mixture.dup;

	//2: create array
	auto mixtureSave2 = new double[mixture.length];
	//2: copy
	mixtureSave2[] = mixture[];
	-------------
	*/
	const(T)[] mixture() @property const
	{
		return _mixture;
	}

	/**
	Returns:
		Const slice to the internal distribution representation.
	Example:
	-------------
	//save slice
	auto distribution = optimizer.distribution;

	//use $(D .dup) or copy to save current distribution

	//1: .dup
	auto distributionSave1 = distribution.dup;

	//2: create array
	auto distributionSave2 = new double[distribution.length];
	//2: copy
	distributionSave2[] = distribution[];
	-------------
	*/
	const(T)[] distribution() @property const
	{
		return _distribution;
	}

	/**
	Set the mixture distribution and calls $(MREF updateMixture)
	Params:
		_distribution = new mixture distribution
	*/
	void distribution(in T[] _distribution) @property
	{
		this._distribution[] = _distribution[];
		updateMixture;
	}

	/**
	Performs optimization.
	Params:
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of mixture and distribution. 
			The delegate must return true when mixture and distribution are acceptable. 
		findRootTolerance = Tolerance for inner optimization. See $(STDREF numeric, findRoot).
	*/
	void optimize(
			scope bool delegate(in T[] mixturePrev, in T[] mixture, in T[] distributionPrev, in T[] distribution) tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		scope T[] mixturePrev = new T[mixture.length];
		scope T[] distributionPrev = new T[distribution.length];
		do
		{
			mixturePrev[] = mixture[];
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
		}
		while(!tolerance(mixturePrev, mixture, distributionPrev, distribution));
	}

	/**
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of $(D objectiveFunction(mixture)) and distribution. 
			The delegate must return true when mixture and distribution are acceptable. 
		findRootTolerance = Tolerance for inner optimization. See $(STDREF numeric, findRoot).
	*/
	void optimize(
			scope T delegate(in T[] mixture) objectiveFunction, 
			scope bool delegate (T objectiveFunctionValuePrev, T objectiveFunctionValue, in T[] distributionPrev, in T[] distribution) tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T objectiveFunctionValuePrev;
		T objectiveFunctionValue = objectiveFunction(mixture);
		scope T[] distributionPrev = new T[distribution.length];
		do
		{
			objectiveFunctionValuePrev = objectiveFunctionValue;
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
			objectiveFunctionValue = objectiveFunction(mixture);
		}
		while(!tolerance(objectiveFunctionValuePrev, objectiveFunctionValue, distributionPrev, distribution));
	}
}
