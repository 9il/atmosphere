/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.math;

import core.stdc.tgmath;

import std.traits;

/** Inverse of the Log Minus Digamma function
 * 
 *   Returns x such $(D log(x) - digamma(x) == y).
 *
 * References:
 *   1. Abramowitz, M., and Stegun, I. A. (1970).
 *      Handbook of mathematical functions. Dover, New York,
 *      pages 258-259, equation 6.3.18.
 * 
 * Authors: Ilya Yaroshenko
 */
real logmdigammaInverse(real y)
{
    import std.mathspecial: logmdigamma;
    import std.numeric: findRoot;
    enum maxY = logmdigamma(real.min_normal);
    static assert(maxY > 0 && maxY <= real.max);
    if (y >= maxY)
        return 1 / y; //lim x->0 (log(x)-digamma(x))*x == 1
    if (y < 0)
        return real.nan;
    if (y < real.min_normal)
        return 0.5 / y; //6.3.18
    if (y > 0)
        // x/2 <= logmdigamma(1 / x) <= x, x > 0
        // calls logmdigamma ~6 times
        return 1 / findRoot((real x) => logmdigamma(1 / x) - y, y,  2*y); // ~6 times
    return y; //NaN
}

unittest {
    import std.range : iota;
    import std.math : approxEqual, nextDown;
    import std.mathspecial: logmdigamma;
    foreach (x; iota(1.3L, 2L, 0.1L))
    {
        assert(logmdigammaInverse(logmdigamma(x)).approxEqual(x));
        assert(logmdigammaInverse(logmdigamma(1/x)).approxEqual(1/x));
    }
    //WolframAlpha, 22.02.2015
    immutable testData = [
        [1.0L, 0.615556766479594378978099158335549201923L],
        [1.0L/8, 4.15937801516894947161054974029150730555L],
        [1.0L/1024, 512.166612384991507850643277924243523243L],
        [0.000500083333325000003968249801594877323784632117L, 1000.0L],
        [1017.644138623741168814449776695062817947092468536L, 1.0L/1024],
    ];
    foreach(test; testData)
        assert(approxEqual(logmdigammaInverse(test[0]), test[1], 1e-15, 0));
    assert(approxEqual(logmdigamma(logmdigammaInverse(1)), 1, 1e-15, 0));
    assert(approxEqual(logmdigamma(logmdigammaInverse(real.min_normal)), real.min_normal, 1e-15, 0));
    assert(approxEqual(logmdigamma(logmdigammaInverse(real.max/2)), real.max/2, 1e-15, 0));
    assert(approxEqual(logmdigammaInverse(logmdigamma(1)), 1, 1e-15, 0));
    assert(approxEqual(logmdigammaInverse(logmdigamma(real.min_normal)), real.min_normal, 1e-15, 0));
    assert(approxEqual(logmdigammaInverse(logmdigamma(real.max/2)), real.max/2, 1e-15, 0));
}
