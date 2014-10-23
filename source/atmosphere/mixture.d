module atmosphere.mixture;


import atmosphere.internal;
import atmosphere.utilities : sumOfLog2s;
import std.range;
import std.traits;
import std.numeric;
///**
//In most cases coordinate descent is much more faster then gradient descent.
//*/
//module atmosphere.stationary;

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
	Params:
		k = number of components
		maxLength = maximal length of each component. In terms of likelihood maximization maxLength is length of a sample.
	*/
	this(size_t k, size_t maxLength)
	{
		_componentsT = SlidingWindow!T(k, maxLength);
		_distribution = new T[k];
		_distribution[] = T(1)/k;
		_mixture = new T[maxLength];
	}

	/**
	Perform k (1) iterations of coordinate (gradient or EM) descent optimization algorithm.
	Params:
	findRootTolerance = Defines an early termination condition. 
			Receives the current upper and lower bounds on the root. 
			The delegate must return true when these bounds are acceptable.
	See_Also:
		$(STDREF numeric, findRoot)
	*/
	abstract void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null);

	/**
	update method is called when mixture changes occur.
	*/
	abstract void update();

final:

	/**
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of $(D objectiveFunction(mixture)) and distribution. 
			The delegate must return true when mixture and distribution are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also:
		$(STDREF numeric, findRoot)
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

	/**
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of $(D objectiveFunction(mixture)). 
			The delegate must return true when mixture and distribution are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also:
		$(STDREF numeric, findRoot)
	*/
	void optimize
	(
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

	/**
	Performs optimization.
	Params:
		objectiveFunction = accepts mixture.
		tolerance = Defines an early termination condition. 
			Receives the current and previous versions of distribution. 
			The delegate must return true when mixture and distribution are acceptable. 
		findRootTolerance = Tolerance for inner optimization.
	See_Also:
		$(STDREF numeric, findRoot)
	*/
	void optimize
	(
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

	///
	size_t length() @property const
	{
		return _componentsT.matrix.width;
	}

	///
	size_t maxLength() @property const
	{
		return _componentsT.matrix.shift;
	}


	///
	void reset()
	{
		_componentsT.reset;
	}

	///
	void popFront()
	{
		_componentsT.popFront;
		updateMixtureBack;
	}

	///
	void popFrontN(size_t n)
	{
		_componentsT.popFrontN(n);
		updateMixtureBackN(n);
	}

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

	package void updateMixturePopBack()
	{
		updateMixturePopBackN(1);
	}

	package void updateMixturePopBackN(size_t n)
	{
		_mixture[0.._componentsT.matrix.width-n] = _mixture[n.._componentsT.matrix.width];
		update();
	}
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
	T = floating point type
*/
class GradientDescent(alias Gradient, T) : MixtureOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;
	private T[] c;

	///
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		xi = new T[maxLength];
		gamma = new T[maxLength];
		c = new T[k];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
		c.destroy;
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		gradientDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, mixture, pi[0.._componentsT.matrix.width], xi[0.._componentsT.matrix.width], gamma[0.._componentsT.matrix.width], c, findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/**
Params:
	Gradient = Gradient of the objective function. $(D Gradient(a, b) should perform b = grad_f(a)).
	T = floating point type
*/
class CoordinateDescent(alias Gradient, T) : MixtureOptimizer!T
{
	private T[] pi;
	private T[] xi;
	private T[] gamma;

	///
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
		xi = new T[maxLength];
		gamma = new T[maxLength];
	}

	~this()
	{
		pi.destroy;
		xi.destroy;
		gamma.destroy;
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIteration!(Gradient, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, mixture, pi[0.._componentsT.matrix.width], xi[0.._componentsT.matrix.width], gamma[0.._componentsT.matrix.width], findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/**
Params:
	PartialDerivative = Partial derivative $(D y) of objective convex function $(D u): $(D du/dω_j = y(ω_j), 1 <= j <= n).
	T = floating point type
*/
class CoordinateDescentPartial(alias PartialDerivative, T) : MixtureOptimizer!T
{
	private T[] pi;

	///
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
		pi = new T[maxLength];
	}

	~this()
	{
		pi.destroy;
	}

	final override void eval(scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null)
	{
		coordinateDescentIterationPartial!(PartialDerivative, T)(cast(Matrix!(const T))_componentsT.matrix, _distribution, _mixture[0.._componentsT.matrix.width], pi[0.._componentsT.matrix.width], findRootTolerance is null ? (a, b) => false : findRootTolerance);
		updateMixture;
	}

	override void update(){};
}


/**
*/
interface LikelihoodMaximization(T)
{
	/**
	See_Also:
		 $(STDREF traits, isCallable)
	*/
	void put(PDFRange, SampleRange)(PDFRange pdfs, SampleRange sample);


	///
	void optimize
	(
		scope bool delegate (T sumOfLog2sValuePrev, T sumOfLog2sValue, in T[] distributionPrev, in T[] distribution) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	);

	///
	void optimize
	(
		scope bool delegate (T sumOfLog2sValuePrev, T sumOfLog2sValue) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	);
}

/**
*/
class CoordinateLikelihoodMaximization(T) : CoordinateDescentPartial!(a => -1/a, T), LikelihoodMaximization!T
{
	///
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
	}

	mixin LikelihoodMaximizationTemplate!T;
}

/**
*/
class GradientLikelihoodMaximization(T) : GradientDescent!((a, b) {foreach(i, ai; a) b[i]=-1/ai;}, T), LikelihoodMaximization!T
{
	///
	this(size_t k, size_t maxLength)
	{
		super(k, maxLength);
	}

	mixin LikelihoodMaximizationTemplate!T;
}

private mixin template LikelihoodMaximizationTemplate(T)
{
	void put(PDFRange, SampleRange)(PDFRange pdfs, SampleRange sample)
	if(isInputRange!PDFRange && hasLength!PDFRange && isCallable!(ElementType!PDFRange))
	in
	{
		assert(pdfs.length == _componentsT.matrix.height);
	}
	body
	{
		super.put(sample.map!(x => pdfs.map!(pdf => pdf(x))));
	}

	void optimize
	(
		scope bool delegate (T sumOfLog2sValuePrev, T sumOfLog2sValue, in T[] distributionPrev, in T[] distribution) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		super.optimize((m => m.sumOfLog2s), tolerance, findRootTolerance);
	}

	void optimize
	(
		scope bool delegate (T sumOfLog2sValuePrev, T sumOfLog2sValue) tolerance,
		scope bool delegate(T a, T b) @nogc nothrow findRootTolerance = null,
	)
	{
		super.optimize((m => m.sumOfLog2s), tolerance, findRootTolerance);
	}
}

unittest {
	alias C0 = CoordinateDescent!((a, b){}, double);
	alias C1 = LikelihoodMaximization!(double);
	alias C10 = GradientLikelihoodMaximization!(double);
	alias C11 = CoordinateLikelihoodMaximization!(double);
	alias C2 = GradientDescent!((a, b){}, double);
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