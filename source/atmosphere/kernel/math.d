/**
This module is not part of the library.
It contains some helper functions for creating Cyclones.
Function can be highly slow and not precise.
If you have more efficient version you can make pull request to project on github.
*/
module atmosphere.kernel.math;

import core.stdc.tgmath;
import std.math : PI;

alias gamma = tgamma;

///Normal PDF
auto normalPDF(F)(in F x)
{
	return (1 / sqrt(2 * PI)) * exp(x ^^ 2 / -2);
}
  

///ditto
auto normalPDF(F)(in F x, in F mean, in F scale)
{
	return normalPDF((x - mean) / scale) / scale;
}


///gamma PDF
auto gammaPDF(F)(in F x, in F shape)
{
	return pow(x, (shape - 1) / (gamma(shape) * x.exp));
}


///ditto
auto gammaPDF(F)(in F x, in F shape, in F scale)
{
	return gammaPDF(x / scale, shape) / scale;
}

