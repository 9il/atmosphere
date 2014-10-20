module atmosphere.mixtureoptimizer;

import std.traits : isFloatingPoint;

/**
Base abstract class for all mixture distribution optimization algorithms.
*/
abstract class MixtureOptimizer(T)
	if(isFloatingPoint!T)
{

package:

	void updateMixture();

public:

	/**
	Returns:
		duplicate of the components
	*/
	abstract T[][] components() @property const;

	/**
	Returns:
		mixture of components corresponding the current mixture distribution
	*/
	abstract const(T)[] mixture() @property const;

	/**
	Returns:
		mixture distribution
	*/
	abstract const(T)[] distribution() @property const;

	/**
	Set the mixture distribution.
	Params:
		_distribution = new mixture distribution
	*/
	abstract void distribution(in T[] _distribution) @property;

	/**
	Perfrorm one iteration of optimization.
	*/
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

final:

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
