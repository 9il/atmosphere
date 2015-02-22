/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

import core.stdc.tgmath;

import std.traits;

int j;
int j2;

/**
Returns x such `log(x) - digamma(x) == y`.
*/
real logmdigammaInverse(real y)
{
	alias T = real;
	import std.mathspecial: logmdigamma;
	import std.numeric: findRoot;
	enum maxY = logmdigamma(T.min_normal);
	static assert(maxY > 0 && maxY < T.infinity);
	int i, i2;
	T f(T x)
	{
		i++; j++;
		return logmdigamma(x)-y;
	}
	T f2(T x)
	{
		i2++; j2++;
		return logmdigamma(1 / x) - y;
	}
	import std.stdio;
	//scope(exit) writeln("i = ", i, " j = ", j);
	//scope(exit) writeln("i2 = ", i2, " j2 = ", j2);
	if(y > maxY)
		return 0;
	if(y < 0)
		return T.nan;
	if(y == 0)
		return T.infinity;
	if(y > 0)
	{
		//auto x = findRoot(&f, 1 / (2 * y), 1/y);
		auto x2 = 1 / findRoot(&f2, y, 2 * y);
		//writeln("x = ", x, " x2 = ", x2);
		return x2;
	}
	return y; //NaN
}

unittest {
	import std.range : iota;
	import std.math : approxEqual;
	import std.mathspecial: logmdigamma;
	foreach(x; iota(1.3, 10.0, 2.0))
	{
		assert(logmdigammaInverse(logmdigamma(x)).approxEqual(x));
		//assert(logmdigammaInverse(logmdigamma(1/x)).approxEqual(1/x));
	}
}
