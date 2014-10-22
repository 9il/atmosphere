module atmosphere.sliding;

import atmosphere.mixtureoptimizer;
import atmosphere.internal;



/**
Params:
	T = floating point type
*/
abstract class SlidingOptimizer(T) : MixtureOptimizer!T, InputRange!T, OutputRange!T
{

	private SlidingWindow!T data;
	private T[] _mixture_;
	private T[] _distribution_;


final:

	/**
	Params:
		k = number of components
		maxN = miximal length of each component. In terms of likelihood maximization n is length of sample.
	*/
	this(size_t k, size_t maxN)
	{
		data = SlidingWindowAllocator!T(k+3, maxN, k+1);
		_distribution_ = new T[k];
		_distribution_[] = T(1)/k;
	}

	/**
	This function accept components and store them in internal transposed form.
	Params:
		_components = components
	*/
	void components(in T[][] _components) @property
	in
	{
		assert(_components.length == _componentsT.width);
		foreach(c; _components)
			assert(c.length == _componentsT.height);
	}
	body
	{
		foreach(j, w; _components)
			foreach(i, e; w)
				_componentsT[i, j] = e;
		updateMixture;
	}

	/**
	Compute internal components using kernels and grid. 
	The pseudocode of this method: $(D components[i, j] = kernel[i](sample[j])).
	Params:
		kernels = Input range of kernels. kernels length equals k. 
		Element type of kernels can be a $(D function) pointer, $(D delegate), 
			$(D struct) or $(D class) with $(D opCall) method.
			Kernel compute feature for each point of grid.
			In terms of likelihood maximization each kernel is a probability density function.
		grid = pointers of components. 
		In terms of likelihood maximization grid is just a sample of length n.
	*/
	void components(Kernels)(Kernels kernels, in T[] grid) @property
		if(isInputRange!Kernels && hasLength!Kernels)
	in
	{
		assert(kernels.length == _componentsT.width);
		assert(grid.length == _componentsT.height);
	}
	body
	{
		auto m = _componentsT;
		foreach(kernel; kernels)
		{
			auto r = m.front;
			foreach(i, e; grid)
				r[i] = kernel(e);
			m.popFront;
		}
		updateMixture;
	}

override:

	inout(Matrix!T) _componentsT() inout @property
	{
		return data.matrix[0..$-1];
	}

	inout(T)[] _mixture() inout @property 
	{
		return data.matrix[$-1];
	}

	inout(T)[] _distribution() inout @property 
	{
		return _distribution_;
	}

	void _distribution(T[] _distribution_) @property
	{
		this._distribution_ = _distribution_;
	}
}

unittest {
	alias f = SlidingWindow!double;
}