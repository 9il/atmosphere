/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

import core.stdc.tgmath;

import std.traits;


/**
Returns x such `log(x) - digamma(x) == y`.
*/
T logmdigammaInverse(T)(T y)
	if(isFloatingPoint!T)
{
	import std.mathspecial: logmdigamma;
	import std.numeric: findRoot;
	enum maxY = logmdigamma(T.min_normal);
	static assert(maxY > 0 && maxY < T.infinity);
	if(y > maxY)
		return 0;
	if(y < 0)
		return T.nan;
	if(y == 0)
		return T.infinity;
	if(y > 0)
		return findRoot((T x)=>logmdigamma(x)-y, y > 1.39325f ? T.min_normal : 0.46163f, 1/y);
	return y; //NaN
}

unittest {
	import std.range : iota;
	import std.math : approxEqual;
	import std.mathspecial: logmdigamma;
	foreach(x; iota(1.3, 10.0, 2.0))
	{
		assert(logmdigammaInverse(logmdigamma(x)).approxEqual(x));
		assert(logmdigammaInverse(logmdigamma(1/x)).approxEqual(1/x));
	}
}
