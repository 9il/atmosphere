/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.statistic;

import core.stdc.tgmath;

import std.traits;
import std.typecons;
import std.math : LN2;
import std.algorithm : map;
//import std.numeric : sumOfLog2s;
import atmosphere.math: sumOfLog2s;

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
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		mean = sample.wfsum() / n;
		meani = sample.map!"1/a".wfsum() / n;
		meanl = T(LN2) * sample.map!(x => cast(Unqual!T)x).sumOfLog2s() / n;
	}

	///
	this(in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		immutable n = weights.wfsum;
		mean = sample.wfsum(weights) / n;
		meani = sample.map!"1/a".wfsum(weights) / n;
		meanl = T(LN2) * sample.map!log2.wfsum(weights) / n;
	}

	///
	C opCast(C : GeneralizedInverseGaussinFixedLambdaStatistic!T)()
	{
		GeneralizedInverseGaussinFixedLambdaStatistic!T stat = void;
		stat.mean = mean;
		stat.meani = meani;
		return stat;
	}
}

unittest {
	alias st = GeneralizedInverseGaussinStatistic!double;
}

/++
Minimal sufficient and complete statistic for the generalized inverse Gaussin disributoin with fixed lambda.
+/
struct GeneralizedInverseGaussinFixedLambdaStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] * sample[j] / Σ weights[j]`
	T mean;
	///`Σ weights[j] / sample[j] / Σ weights[j]`
	T meani;

	///
	this(in T[] sample)
	in {
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		mean = sample.wfsum() / n;
		meani = sample.map!"1/a".wfsum() / n;
	}

	///
	this(in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		immutable n = weights.wfsum;
		mean = sample.wfsum(weights) / n;
		meani = sample.map!"1/a".wfsum(weights) / n;
	}
}

unittest {
	alias st = GeneralizedInverseGaussinFixedLambdaStatistic!double;
}


/++
Minimal sufficient and complete statistic for the generalized gamma disributoin with fixed `power` parameter.
+/
struct GeneralizedGammaFixedPowerStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] * sample[j] ^^ power / Σ weights[j]`
	T meanp;
	///`Σ weights[j] * log(sample[j]) / Σ weights[j]`
	T meanl;

	///
	this(T power, in T[] sample)
	in {
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		meanp = sample.map!(x => x.pow(power)).wfsum / n;
		meanl = T(LN2) * sample.map!log2.wfsum / n;
	}

	///
	this(T power, in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		immutable n = weights.wfsum;
		meanp = sample.map!(x => x.pow(power)).wfsum(weights) / n;
		meanl = T(LN2) * sample.map!log2.wfsum(weights) / n;
	}
}

unittest {
	alias st = GeneralizedGammaFixedPowerStatistic!double;
}


/++
Minimal sufficient and complete statistic for the gamma disributoin.
+/
struct GammaStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] * sample[j] / Σ weights[j]`
	T mean;
	///`Σ weights[j] * log(sample[j]) / Σ weights[j]`
	T meanl;

	///
	this(in T[] sample)
	in {
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		mean = sample.wfsum / n;
		meanl = T(LN2) * sample.map!(x => cast(Unqual!T)x).sumOfLog2s / n;
	}

	///
	this(in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		immutable n = weights.wfsum;
		mean = sample.wfsum(weights) / n;
		meanl = T(LN2) * sample.map!log2.wfsum(weights) / n;
	}
}

unittest {
	alias st = GammaStatistic!double;
}


/++
Minimal sufficient and complete statistic for the inverse-gamma disributoin.
+/
struct InverseGammaStatistic(T)
	if(isFloatingPoint!T)
{
	///`Σ weights[j] / sample[j] / Σ weights[j]`
	T meani;
	///`Σ weights[j] * log(sample[j]) / Σ weights[j]`
	T meanl;

	///
	this(in T[] sample)
	in {
		assert(positiveSampleCheck(sample));
	}
	body {
		immutable n = sample.length;
		meani = sample.map!"1/a".wfsum / n;
		meanl = T(LN2) * sample.map!(x => cast(Unqual!T)x).sumOfLog2s / n;
	}

	///
	this(in T[] sample, in T[] weights)
	in {
		assert(positiveSampleCheck(sample, weights));
	}
	body {
		immutable n = weights.wfsum;
		meani = sample.map!"1/a".wfsum(weights) / n;
		meanl = T(LN2) * sample.map!log2.wfsum(weights) / n;
	}
}

package:

unittest {
	alias st = InverseGammaStatistic!double;
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

auto wfsum(Range)(Range sample)
{
	import atmosphere.summation;
	return sample.fsum!(Summation.KB2);
}

unittest {
	assert(wfsum([1.0, 2]) == 3);
}

T wfsum(Range, T)(Range sample, in T[] weights)
{
	import atmosphere.summation;
	import std.range;
	assert(sample.length == weights.length);
	Summator!(T, Summation.KB2) s = 0;
	foreach(i, w; weights)
	{
		s += sample.front * w;
		sample.popFront;
	}
	return s.sum;
}

unittest {
	assert(wfsum([1.0, 2], [3.0, 4]) == 11);
}
