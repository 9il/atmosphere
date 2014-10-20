module atmosphere;
/**
This module contains hight level implementation of separating mixtures numeric methods.

You can create derived classes from non-abstract classes.

Copyright: Ilya Yaroshenko 2014.

License: MIT.

Authors: Ilya Yaroshenko
*/
import core.stdc.tgmath;


import std.numeric : dotProduct;
import std.algorithm : sum;
import std.range : hasLength, isInputRange, ElementType;
import std.traits : isFloatingPoint, Unqual;

import atmosphere.internal;

public import atmosphere.normalvariancemeanmixture;
public import atmosphere.mixtureoptimizer;
public import atmosphere.stationary;
public import atmosphere.sliding;


///
unittest {
	auto componentsCount = 50;
	auto sampleCount = 10000;
	auto optimizer = new LikelihoodMaximizationCoordinate!double(componentsCount, sampleCount);

	///distribution is uniform
	foreach(e; optimizer.distribution)	
		assert(e == 1.0/componentsCount);


}