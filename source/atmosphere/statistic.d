/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.statistic;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.algorithm : fsum = sum;


/++
Minimal sufficient and complete statistic for the generalized inverse Gaussin disributoin.
+/
struct GeneralizedInverseGaussinStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] * sample[j] / Σ weights[j]`
	T mean;
	///`Σ weights[j] / sample[j] / Σ weights[j]`
	T meani;
	///`Σ weights[j] * log(sample[j]) / Σ weights[j]`
	T meanl;

	///
	this(in T[] sample)
	in {
	}
	body {
		import std.algorithm : map;
		immutable n = sample.length;
		mean = sample.fsum() / n;
		meani = sample.map!"1/a".fsum() / n;
		meanl = T(LN2) * sample.sumOfLog2s() / n;
	}

	///
	this(in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		import std.algorithm : map, zip;
		immutable n = weights.fsum;
		mean = sample.wfsum(weights) / n;
		meani = sample.map!"1/a".wfsum(weights) / n;
		meanl = T(LN2) * sample.map!log2.wfsum(weights) / n;
	}

	///
	this(T mean, T meani, T meanl)
	{
		this.mean = mean;
		this.meani = meani;
		this.meanl = meanl;
	}
}


struct GeneralizedGammaFixedPowerStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] * sample[j] ^^ power / Σ weights[j]`
	T meanp;
	///`Σ weights[j] * log(sample[j]) / Σ weights[j]`
	T meanl;
	///`Σ weights[j] * log(sample[j]) * sample[j] ^^ power / sample[j] / Σ weights[j]`
	T meanlp;

	///
	this(T power, in T[] sample)
	in {
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		mean = sample.fsum() / n;
		meani = sample.map!"1/a".fsum() / n;
		meanl = T(LN2) * sample.sumOfLog2s() / n;
	}

	///
	this(T power, in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		Summator!T n = 0, a = 0, b = 0, c = 0;
		foreach(i, x; sample)
		{
			immutable w = weights[i];
			immutable l = log2(x);
			immutable y = x.pow(power);
			n += w;
			a += w * y;
			b += w * l;
			c += w * (l * y);
		}
		b *= T(LN2);
		c *= T(LN2);
		meanp = a / n;
		meanl = b / n;
		meanlp = c / n;
	}

	///
	this(T mean, T meani, T meanl)
	{
		this.mean = mean;
		this.meani = meani;
		this.meanl = meanl;
	}
}

bool positiveSampleCheck(T)(in T[] sample)
{
	import std.algorithm : all;
	return sample.all!"a > 0 && isNormal(a)";
}

bool positiveSampleCheck(T)(in T[] sample, in T[] weights)
{
	import std.algorithm : all, any;
	return 
	   weights.length == sample.length
	&& sample.all!"a > 0 && isNormal(a)"
	&& weights.all!"a >= 0 && isFinite(a)"
	&& weights.any!"a > 0";
}

private T wfsum(Range)(in Range sample, in T[] weights)
{
	assert(sample.length == weights.length);
	Summator!T s = 0;
	foreach(x; sample)
	{
		s += x * weights.front;
		weights.popFront;
	}
	return s.sum;
}