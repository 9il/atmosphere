module atmosphere.mixtureoptimizer;


import atmosphere.internal;
import std.range;
import std.traits;
import std.numeric;
///**
//In most cases coordinate descent is much more faster then gradient descent.
//*/
//module atmosphere.stationary;

import atmosphere.mixtureoptimizer;
import atmosphere.internal;


/**
Params:
	T = floating point type
*/
abstract class MixtureOptimizer(T)
{
	package SlidingWindow!T _componentsT;
	package T[] _distribution;
	package T[] _mixture;


	/**
	Perform one iteration of optimization.
	*/
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

	///update method is called when mixture changed
	void update() {}

	//inout(Matrix!T) _componentsT() inout @property;
	//inout(T)[] _mixture() inout @property;
	//inout(T)[] _distribution() inout @property;
	//void _distribution(T[] _distribution_) @property;


	package void updateMixture()
	{
		mix(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture[0.._componentsT.matrix.width]);
		update();
	}

	package void updateMixtureBack()
	{
		_mixture[_componentsT.matrix.width-1] = dotProduct(_distribution, _componentsT.back);
		update();
	}

	package void updateMixtureBackN(size_t n)
	{
		mix(cast(Matrix!(const T))_componentsT[$-n..$].matrix, _distribution, _mixture[0.._componentsT.matrix.width]);
		update();
	}

//final:


	//Matrix!(const(T)) componentsT() const
	//{
	//	return cast(typeof(return))_components.matrix;
	//}
	

	/**
	Returns:
		Const internal components representation.
	*/
	Matrix!(const(T)) components() const
	{
		return cast(typeof(return))_componentsT.matrix;
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
		return _mixture[0.._componentsT.length];
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
		T[] distributionPrev = new T[distribution.length];
		do
		{
			objectiveFunctionValuePrev = objectiveFunctionValue;
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
			objectiveFunctionValue = objectiveFunction(mixture);
		}
		while(!tolerance(objectiveFunctionValuePrev, objectiveFunctionValue, distributionPrev, distribution));
	}

	///ditto
	void optimize(
			scope T delegate(in T[] mixture) objectiveFunction, 
			scope bool delegate (T objectiveFunctionValuePrev, T objectiveFunctionValue) tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T objectiveFunctionValuePrev;
		T objectiveFunctionValue = objectiveFunction(mixture);
		do
		{
			objectiveFunctionValuePrev = objectiveFunctionValue;
			eval(findRootTolerance);
			objectiveFunctionValue = objectiveFunction(mixture);
		}
		while(!tolerance(objectiveFunctionValuePrev, objectiveFunctionValue));
	}

	///ditto
	void optimize(
			scope bool delegate (in T[] distributionPrev, in T[] distribution) tolerance,
			scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
		)
	{
		T[] distributionPrev = new T[distribution.length];
		do
		{
			distributionPrev[] = distribution[];
			eval(findRootTolerance);
		}
		while(!tolerance(distributionPrev, distribution));
	}

//final:

	/**
	Params:
		k = number of components
		maxN = length of each component. In terms of likelihood maximization maxN is length of a sample.
	*/
	this(size_t k, size_t maxN)
	{
		_componentsT = SlidingWindow!T(k, maxN);
		_distribution = new T[k];
		_distribution[] = T(1)/k;
		_mixture = new T[maxN];
	}

	//alias components = super.components;
	/**
	*/
	void put(Range)(Range r)
	if(isInputRange!Range && hasLength!Range && isNumeric!(ElementType!Range))
	in
	{
		assert(_componentsT.matrix.length == r.length);
		assert(_componentsT.matrix.width < _componentsT.matrix.shift);
	}
	body
	{
		_componentsT.put(r);
		updateMixtureBack;
	}

	void put(RangeOfRanges)(RangeOfRanges ror)
	if(isInputRange!RangeOfRanges && hasLength!RangeOfRanges && 
		isInputRange!(ElementType!RangeOfRanges) && hasLength!(ElementType!RangeOfRanges) && 
		isNumeric!(ElementType!(ElementType!RangeOfRanges)))
	in
	{
		assert(_componentsT.matrix.length == ror.front.length);
		assert(_componentsT.matrix.width + ror.length <= _componentsT.matrix.shift);
	}
	body
	{
		const n = ror.length;
		.put(_componentsT, ror);
		updateMixtureBackN(n);
	}

	void reset()
	{
		_componentsT.reset;
	}

	/*
	Compute internal components using kernels and grid. 
	The pseudocode of this method: $(D components[i, j] = kernel[i](sample[j])).
	Params:
		kernels = Input range of kernels. kernels length equals k. 
		Element type of kernels can be a $(D function) pointer, $(D delegate), 
			$(D struct) or $(D class) with $(D opCall) method.
			Kernel compute feature for each point of grid.
			In terms of likelihood maximization each kernel is a probability density function.
		grid = pointers of components. 
		In terms of likelihood maximization grid is a sample of length n.
	*/
	//void components(Kernels)(Kernels kernels, in T[] grid)
	//	if(isInputRange!Kernels && hasLength!Kernels)
	//in
	//{
	//	assert(kernels.length == components.front.length, format("kernels.length(%s) != this.components.front.length(%s)", kernels.length, components.front.length));
	//	assert(grid.length == components.length, format("kernels.length(%s) != this.components.length(%s)", grid.length, components.length));
	//}
	//body
	//{
	//	import std.algorithm;
	//	auto m = _componentsT;
	//	foreach(kernel; kernels)
	//	{
	//		grid.map!(a => kernel(a)).copy(m.front);
	//		m.popFront;
	//	}
	//	updateMixture;
	//}

//override:

	//inout(Matrix!T) _componentsT() inout @property
	//{
	//	return _componentsT_;
	//}

	//inout(T)[] _mixture() inout @property 
	//{
	//	return _mixture_;
	//}

	//inout(T)[] _distribution() inout @property 
	//{
	//	return _distribution_;
	//}

	//void _distribution(T[] _distribution_) @property
	//{
	//	this._distribution_ = _distribution_;
	//}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
	T = floating point type
*/
final class GradientDescent(alias Gradient, T) : MixtureOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;
	private T[] c;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
		xi = new T[n];
		gamma = new T[n];
		c = new T[k];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
		c.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, xi, gamma, c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
	T = floating point type
*/
final class CoordinateDescent(alias Gradient, T) : MixtureOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
		xi = new T[n];
		gamma = new T[n];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, xi, gamma, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}


/**
Params:
	PartialDerivative = Partial derivative of the objective function. 
		There are optimization reasons for coordinate descent to use partial derivative instead of gradient.
		It is applicable if and only if all partial derivatives are equal.
		$(D b = PartialDerivative(a) should perform b = (d f/ d x_i)(a)).
	T = floating point type
*/
final class CoordinateDescentPartial(alias PartialDerivative, T) : MixtureOptimizer!T
{
	private T[] pi;

	///
	this(size_t k, size_t n)
	{
		super(k, n);
		pi = new T[n];
	}

	~this()
	{
		pi.destroy;
	}

	override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture, pi, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}
}

/**
Likelihood gradient maximization.
----------
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a, T);
----------
*/
alias LikelihoodMaximizationCoordinate(T) = CoordinateDescentPartial!(a => -1/a, T);


/**
Likelihood gradient maximization.
----------
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;});
----------
*/
alias LikelihoodMaximizationGradient(T) = GradientDescent!((a, b) {foreach(i, ai; a) b[i] = -1/ai;}, T);


unittest {
	alias C0 = CoordinateDescent!((a, b){}, double);
	alias C1 = LikelihoodMaximizationGradient!(double);
	alias C2 = LikelihoodMaximizationCoordinate!(double);
}


/////Example:
//unittest {
//	import std.range;

//	//import atmosphere.stationary;
//	auto optimizer = new LikelihoodMaximizationCoordinate!double(10, 100);

//	auto components = optimizer.MixtureOptimizer.components;
//	alias ComponentsType = typeof(components);
//	alias ComponentType = ElementType!ComponentsType;
	
//	static assert(isRandomAccessRange!ComponentType);
//	static assert(isRandomAccessRange!ComponentsType);
//	static assert(hasAssignableElements!ComponentType == false);
//	static assert(hasAssignableElements!ComponentsType == false);
//}