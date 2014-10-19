module atmosphere;

import core.stdc.tgmath;


import std.numeric : dotProduct;
import std.algorithm : sum;
import std.range : hasLength, isInputRange, ElementType;
import std.traits : isFloatingPoint, Unqual;

import atmosphere.internal;

/**
*/
abstract class MixtureOptimizer(T)
	if(isFloatingPoint!T)
{

private:

	void updateMixture();

public:

	///
	abstract T[][] components() @property const;

	///
	abstract const(T)[] mixture() @property const;

	///
	abstract const(T)[] distribution() @property const;

	///
	abstract void distribution(in T[] _distribution) @property;

	///
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

	///
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

	///
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


/**
*/
abstract class StationaryOptimizer(T) : MixtureOptimizer!T
{

private:

	Matrix!T _componentsT;
	T[] _distribution;
	T[] _mixture;
	scope T[] pi;

	void updateMixture()
	{
		mix(cast(Matrix!(const T))_componentsT, _distribution, _mixture);
	}

public:

	///
	this(size_t k, size_t n)
	{
		_componentsT = Matrix!T(k, n);
		_distribution = new T[k];
		_distribution[] = T(1)/k;
		_mixture = new T[n];
		pi = new T[n];
	}

	///
	void components(in T[][] _components) @property
	{
		foreach(j, w; _components)
			foreach(i, e; w)
				_componentsT[i, j] = e;
	}

	/**
	*/
	void components(Kernels)(Kernels kernels, in T[] grid) @property
		if(isInputRange!Kernels && hasLength!Kernels)
	{
		auto m = _componentsT;
		foreach(kernel; kernels)
		{
			auto r = m.front;
			foreach(i, e; grid)
				r[i] = kernel(e);
			m.popFront;
		}
	}

override:

	T[][] components() @property const
	{
		return _componentsT.transpose.arrays;
	}

	const(T)[] mixture() @property const
	{
		return _mixture;
	}

	const(T)[] distribution() @property const
	{
		return _distribution;
	}

	void distribution(in T[] _distribution) @property
	{
		this._distribution[] = _distribution[];
	}
}




/**
*/
final class GradientDescent(alias Gradient, T) : StationaryOptimizer!T
{

private:

	scope T[] xi;
	scope T[] gamma;
	scope T[] c;

public:

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		xi = new T[n];
		gamma = new T[n];
		c = new T[k];
	}

override:

	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
*/
class NormalVarianceMeanMixtureSeparator(T) : StationaryOptimizer!T
{

private:

	scope T[] _sample;
	scope const T[] _grid;

	T _mean;
	T _alpha;
	T _sumOfLog2s;

	void updateAlpha()
	in
	{
		import std.string;
		assert(_distribution.length == _grid.length, format("%s %s", _distribution, _grid));
	}
	body
	{
		_alpha =  _mean / dotProduct(_distribution, _grid);
	}

	void updateSumOfLog2s()
	{
		_sumOfLog2s = _mixture.sumOfLog2s;
	}

	void updateComponents()
	{
		import std.algorithm : map;

		static struct Kernel
		{
			T alphau;
			T sqrtu;

			this(T alphau, T sqrtu)
			{
				this.alphau = alphau;
				this.sqrtu = sqrtu;
			}

			T opCall(T x) const
			{
				immutable y = (x - alphau) / sqrtu;
				return exp(y * y / -2) / sqrtu;
			}
		}

		auto kernels = _grid.map!(u => Kernel(alpha*u, sqrt(u)));
		components(kernels, _sample);
	}

	void updateAll()
	{
		updateAlpha;
		updateComponents;
		updateMixture;
		updateSumOfLog2s;
	}

public:

	///
	this(in T[] _grid, in T[] _sample)
	in
	{
		assert(_grid.length);
		assert(_sample.length);
	}
	body
	{
		super(_grid.length, _sample.length);
		this._grid = _grid.dup;
		this._sample = new T[_sample.length];
		sample = _sample;
	}

final:

	///
	T alpha() @property const
	{
		return _alpha;
	}

	///
	T mean() @property const
	{
		return _mean;
	}

	///
	T sumOfLog2s() @property const
	{
		return _sumOfLog2s;
	}

	///
	void sample(in T[] _sample) @property
	{
		this._sample[] = _sample[];
		_mean = sample.sum/sample.length;
		updateAll;
	}

	///
	const(T)[] sample() @property const
	{
		return _sample;
	}

	///
	const(T)[] grid() @property const
	{
		return _grid;
	}

	///
	void optimize(
			scope bool delegate (
				T alphaPrev, 
				T alpha, 
				T sumOfLog2sValuePrev, 
				T sumOfLog2sValue, 
				in T[] distributionPrev, 
				in T[] distribution) 
			tolerance
		)
	{
		T sumOfLog2sPrev, alphaPrev;
		scope T[] distributionPrev = new T[_distribution.length];
		do
		{
			sumOfLog2sPrev = _sumOfLog2s;
			alphaPrev = _alpha;
			assert(distribution.length == distributionPrev.length);
			distributionPrev[] = _distribution[];
			eval;
		}
		while(!tolerance(alphaPrev, _alpha, sumOfLog2sPrev, _sumOfLog2s, distributionPrev, _distribution));
	}

override:

	const(T)[] distribution() const @property
	{
		return _distribution;
	}

	void distribution(in T[] _distribution) @property
	{
		super.distribution(_distribution);
		updateAll;
	}
}


/**
*/
final class NormalVarianceMeanMixtureEMSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{

private:
	scope T[] c;

public:
	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
		c = new T[_grid.length];
	}

override:
	
	///
	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		EMIteration!
			((a, b) {foreach(i, ai; a) b[i] = 1/ai;}, T)
			(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, c);
		updateAll;
	}
}

/**
*/
final class NormalVarianceMeanMixtureEMAndGradientSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{

private:
	scope T[] xi;
	scope T[] gamma;
	scope T[] c;

public:
	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
		xi = new T[_sample.length];
		gamma = new T[_sample.length];
		c = new T[_grid.length];
	}

override:
	
	///
	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!
			((a, b) {foreach(i, ai; a) b[i] = -1/ai;}, T)
			(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateAll;
	}
}

/**
*/
final class NormalVarianceMeanMixtureEMAndCoordinateSeparator(T) : NormalVarianceMeanMixtureSeparator!T
{
	///
	this(in T[] _grid, in T[] _sample)
	{
		super(_grid, _sample);
	}

override:
	
	///
	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!
			(a => -1/a, T)
			(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateAll;
	}
}


unittest
{
	alias C = NormalVarianceMeanMixtureEMAndCoordinateSeparator!double;
}

/**
*/
final class CoordinateDescent(alias Gradient, T) : StationaryOptimizer!T
{

private:

	scope T[] xi;
	scope T[] gamma;

public:

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		xi = new T[n];
		gamma = new T[n];
	}

override:

	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, xi, gamma, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**

*/
final class CoordinateDescentPartial(alias PartialDerivative, T) : StationaryOptimizer!T
{
	///
	this(size_t k, size_t n)
	{
		super(k, n);
	}

override:

	void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_componentsT, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
----------
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a);
----------
*/
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a, T);


/**
----------
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;});
----------
*/
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;}, T);

/**
Returns true if Class use coordinate descend.
*/
template isCoordinateOprimization(Class) 
{
	enum bool isCoordinateOprimization = 
	is(Class : CoordinateDescent!(Gradient, T), alias Gradient, T) ||
	is(Class : CoordinateDescentPartial!(PartialDerivative, T), alias PartialDerivative, T);
}

///
unittest
{
	static assert(isCoordinateOprimization!(LikelihoodMaximizationCoordinate!double) == true);
	static assert(isCoordinateOprimization!(LikelihoodMaximizationGradient!double) == false);
}


private:

/**
Computes accurate sum of binary logarithms of input range $(D r).
TODO: Delete this with DMD 2.068.
 */
ElementType!Range sumOfLog2s(Range)(Range r) 
    if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
	import std.math : frexp; 

    long exp = 0;
    Unqual!(typeof(return)) x = 1; 
    foreach (e; r)
    {
        if (e < 0)
            return typeof(return).nan;
        int lexp = void;
        x *= frexp(e, lexp);
        exp += lexp;
        if (x < 0.5) 
        {
            x *= 2;
            exp--;
        }
    }
    return exp + log2(x); 
}