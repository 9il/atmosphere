/**
This module contains basic summation algorithms.

License: $(LINK2 http://boost.org/LICENSE_1_0.txt, Boost License 1.0).

Authors: $(WEB 9il.github.io, Ilya Yaroshenko)

Source: $(PHOBOSSRC std/numeric/_summation.d)
*/
module atmosphere.summation;

import std.traits;
import std.typecons;
import std.range;
import std.math;

private template SummationType(F)
    if (isFloatingPoint!F || isComplex!F)
{
    version(X86) //workaround for Issue 13474
    {
        static if (!is(Unqual!F == real) && (isComplex!F || !is(Unqual!(typeof(F.init.re)) == real)))
            pragma(msg, "Warning: Summation algorithms on x86 use 80bit representation for single and double floating point numbers.");
        static if (isComplex!F)
        {
            import std.complex : Complex;
            alias SummationType = Complex!real;
        }
        else
            alias SummationType = real;
    }
    else
        alias SummationType = F;
}

/++
Computes accurate sum of binary logarithms of input range $(D r).
+/
ElementType!Range sumOfLog2s(Range)(Range r)
    if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
    long exp;
    Unqual!(typeof(return)) x = 1;
    foreach (e; r)
    {
        if (e < 0)
            return typeof(return).nan;
        int lexp = void;
        x *= frexp(e, lexp);
        exp += lexp;
        if (x < 0.5)
        {
            x *= 2;
            exp--;
        }
    }
    return exp + log2(x);
}

///
unittest
{
    import std.math, std.numeric;
    assert(sumOfLog2s(new double[0]) == 0);
    assert(sumOfLog2s([0.0L]) == -real.infinity);
    assert(sumOfLog2s([-0.0L]) == -real.infinity);
    assert(sumOfLog2s([2.0L]) == 1);
    assert(sumOfLog2s([-2.0L]).isNaN());
    assert(sumOfLog2s([real.nan]).isNaN());
    assert(sumOfLog2s([-real.nan]).isNaN());
    assert(sumOfLog2s([real.infinity]) == real.infinity);
    assert(sumOfLog2s([-real.infinity]).isNaN());
    assert(sumOfLog2s([ 0.25, 0.25, 0.25, 0.125 ]) == -9);
}


/++
Computes geometric mean of input range $(D r).
+/
ElementType!Range geometricMean(Range)(Range r)
    if (isInputRange!Range && isFloatingPoint!(ElementType!Range))
{
    static if (hasLength!Range)
    {
        immutable n = r.length;
        immutable s = r.sumOfLog2s();
    }
    else
    {
        import std.range : tee;
        size_t n;
        immutable s = r.tee!((_){n++;}).sumOfLog2s();
    }
    return exp2(s / n);
}

///
unittest {
    assert(geometricMean([2.0, 8.0]) == 4.0);
}


/++
Summation algorithms.
+/
enum Summation
{
    /++
    Fast summation algorithm.
    +/
    Fast,

    /++
    Naive algorithm (one by one).
    +/
    Naive,

    /++
    $(LUCKY Pairwise summation) algorithm. Range must be a finite sliceable range.
    +/
    Pairwise,

    /++
    $(LUCKY Kahan summation) algorithm.
    +/
    Kahan,

    /++
    $(LUCKY Kahan-Babuška-Neumaier summation algorithm). $(D KBN) gives more accurate results then $(D Kahan).
    +/
    KBN,

    /++
    $(LUCKY Generalized Kahan-Babuška summation algorithm), order 2. $(D KB2) gives more accurate results then $(D Kahan) and $(D KBN).
    +/
    KB2,

    /++
    Precise summation algorithm.
    The value of the sum is rounded to the nearest representable
    floating-point number using the $(LUCKY round-half-to-even rule).
    Result can be differ from the exact value on $(D X86), $(D nextDown(proir) <= result &&  result <= nextUp(proir)).
    The current implementation re-establish special value semantics across iterations (i.e. handling -inf + inf).

    References: $(LINK2 http://www.cs.cmu.edu/afs/cs/project/quake/public/papers/robust-arithmetic.ps,
        "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates", Jonathan Richard Shewchuk),
        $(LINK2 http://bugs.python.org/file10357/msum4.py, Mark Dickinson's post at bugs.python.org).
    +/
    /+
    Precise summation function as msum() by Raymond Hettinger in
    <http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/393090>,
    enhanced with the exact partials sum and roundoff from Mark
    Dickinson's post at <http://bugs.python.org/file10357/msum4.py>.
    See those links for more details, proofs and other references.
    IEEE 754R floating point semantics are assumed.
    +/
    Precise,

    /++
    Performs $(D Pairwise) summation for random access ranges. 
    Otherwise performs $(D KBN) summation of floating point and complex numbers and 
    $(D Kahan) summation of user defined types.
    +/
    Appropriate,
}


/++
Computes sum of range.
+/
template fsum(F, Summation summation = Summation.Precise)
    if (isFloatingPoint!F && isMutable!F)
{
    F fsum(Range)(Range r)
        if (isSummable!(Range, F))
    {
        alias sum = Algo!(Range, summation);
        return sum!(Range, Unqual!F)(r);
    }

    F fsum(Range)(F seed, Range r)
        if (isSummable!(Range, Unqual!F))
    {
        alias sum = Algo!(Range, summation);
        return sum!(Range, Unqual!F)(r, seed);
    }
}

///ditto
template fsum(Summation summation = Summation.Precise)
{
    ForeachType!Range fsum(Range)(Range r)
        if (isSummable!(Range, Unqual!(ForeachType!Range)))
    {
        alias sum = Algo!(Range, summation);
        return sum!(Range, Unqual!(typeof(return)))(r);
    }

    F fsum(F, Range)(F seed, Range r)
        if (isSummable!(Range, Unqual!F))
    {
        alias sum = Algo!(Range, summation);
        return sum!(Range, F)(r, seed);
    }
}

///
unittest {
    import std.algorithm;
    auto ar = [1, 1e100, 1, -1e100].map!(a => a*10000);
    const r = 20000;
    assert(r == ar.fsum!(Summation.KBN));
    assert(r == ar.fsum!(Summation.KB2));
    assert(r == ar.fsum); //Summation.Precise
}

// FIXME
// Fails for 32bit systems.
// See also https://issues.dlang.org/show_bug.cgi?id=13474#c7
// and https://github.com/D-Programming-Language/phobos/pull/2513
version(X86)
{

}
else
{
    unittest {
        import std.algorithm;
        auto ar = [1, 1e100, 1, -1e100].map!(a => a*10000);
        const r = 20000;
        assert(r != ar.fsum!(Summation.Naive));
        assert(r != ar.fsum!(Summation.Kahan));
        assert(r != ar.fsum!(Summation.Pairwise));
    }
}

///
unittest
{
    import std.math, std.algorithm, std.range;
    auto ar = 1000
        .iota
        .map!(n => 1.7L.pow(n+1) - 1.7L.pow(n))
        .array
        ;
    //Summation.Precise is default
    real d = 1.7L.pow(1000);
    import std.stdio;
    assert(fsum(ar.chain([-d])) == -1);
    assert(fsum(-d, ar.retro) == -1);
}

/++
$(D Naive), $(D Pairwise) and $(D Kahan) algorithms can be used for user defined types.
+/
unittest {
    static struct Quaternion(F)
        if (isFloatingPoint!F)
    {
        F[3] array;

        /// + and - operator overloading
        Quaternion opBinary(string op)(auto ref Quaternion rhs) const
            if (op == "+" || op == "-")
        {
            Quaternion ret = void;
            foreach (i, ref e; ret.array)
                mixin("e = array[i] "~op~" rhs.array[i];");
            return ret;
        }

        /// += and -= operator overloading
        Quaternion opOpAssign(string op)(auto ref Quaternion rhs)
            if (op == "+" || op == "-")
        {
            foreach (i, ref e; array)
                mixin("e "~op~"= rhs.array[i];");
            return this;
        }

        ///constructor with single FP argument
        this(F f)
        {
            array[] = f;
        }

        ///assigment with single FP argument
        void opAssign(F f)
        {
            array[] = f;
        }
    }

    Quaternion!double q, p, r;
    q.array = [0, 1, 2];
    p.array = [3, 4, 5];
    r.array = [3, 5, 7];

    assert(r == [p, q].fsum!(Summation.Naive));
    assert(r == [p, q].fsum!(Summation.Pairwise));
    assert(r == [p, q].fsum!(Summation.Kahan));
}

/++
All summation algorithms available for complex numbers.
+/
unittest
{
    import std.complex;
    Complex!double[] ar = [complex(1.0, 2), complex(2, 3), complex(3, 4), complex(4, 5)];
    Complex!double r = complex(10, 14);
    assert(r == ar.fsum!(Summation.Fast));
    assert(r == ar.fsum!(Summation.Naive));
    assert(r == ar.fsum!(Summation.Pairwise));
    assert(r == ar.fsum!(Summation.Kahan));
    assert(r == ar.fsum!(Summation.KBN));
    assert(r == ar.fsum!(Summation.KB2));
    assert(r == ar.fsum); //Summation.Precise
}


/++
Output range for summation.
$(D Precise), $(D KB2), $(D KBN) and $(D Kahan) algorithms are supported.
+/
struct Summator(T, Summation summation = Summation.Precise)
    if (
        isMutable!T && (
            summation == Summation.Precise && isFloatingPoint!T
             || summation == Summation.Kahan && isSummable!T
             || (summation == Summation.KBN || summation == Summation.KB2) 
                   && (isFloatingPoint!T || isComplex!T)
            )
        )
{
    static if (summation == Summation.Kahan)
        alias F = T;
    else
        alias F = SummationType!T;

    static if (summation == Summation.Precise)
    {
        import std.internal.scopebuffer;

    private:
        enum F M = (cast(F)(2)) ^^ (T.max_exp - 1);
        F[32] scopeBufferArray = void;
        ScopeBuffer!F partials;
        //sum for NaN and infinity.
        F s;
        //Overflow Degree. Count of 2^^F.max_exp minus count of -(2^^F.max_exp)
        sizediff_t o;


        /++
        Compute the sum of a list of nonoverlapping floats.
        On input, partials is a list of nonzero, nonspecial,
        nonoverlapping floats, strictly increasing in magnitude, but
        possibly not all having the same sign.
        On output, the sum of partials gives the error in the returned
        result, which is correctly rounded (using the round-half-to-even
        rule).
        Two floating point values x and y are non-overlapping if the least significant nonzero
        bit of x is more significant than the most significant nonzero bit of y, or vice-versa.
        +/
        static F partialsReduce(F s, in F[] partials)
        in
        {
            debug(numeric) assert(!partials.length || s.isFinite);
        }
        body
        {
            bool _break = void;
            foreach_reverse(i, y; partials)
            {
                s = partialsReducePred(s, y, i ? partials[i-1] : 0, _break);
                if (_break)
                    break;
                debug(numeric) assert(s.isFinite);
            }
            return s;
        }

        static F partialsReducePred(F s, F y, F z, out bool _break)
        out(result)
        {
            debug(numeric) assert(result.isFinite);
        }
        body
        {
            F x = s;
            s = x + y;
            F d = s - x;
            F l = y - d;
            debug(numeric)
            {
                assert(x.isFinite);
                assert(y.isFinite);
                assert(s.isFinite);
                assert(fabs(y) < fabs(x));
            }
            if (l)
            {
            //Make half-even rounding work across multiple partials.
            //Needed so that sum([1e-16, 1, 1e16]) will round-up the last
            //digit to two instead of down to zero (the 1e-16 makes the 1
            //slightly closer to two). Can guarantee commutativity.
                if (z && !signbit(l * z))
                {
                    l *= 2;
                    x = s + l;
                    F t = x - s;
                    if (l == t)
                        s = x;
                }
                _break = true;
            }
            return s;
        }

        //Returns corresponding infinity if is overflow and 0 otherwise.
        F overflow()
        {
            if (o == 0)
                return 0;
            if (partials.length && (o == -1 || o == 1)  && signbit(o * partials[$-1]))
            {
                // problem case: decide whether result is representable
                F x = o * M;
                F y = partials[$-1] / 2;
                F h = x + y;
                F d = h - x;
                F l = (y - d) * 2;
                y = h * 2;
                d = h + l;
                F t = d - h;
                version(X86)
                {
                    if (!.isInfinity(cast(T)y) || !.isInfinity(sum()))
                        return 0;
                }
                else
                {
                    if (!.isInfinity(cast(T)y) || partials.length > 1 && !signbit(l * partials[$-2]) && t == l)
                        return 0;
                }
            }
            return F.infinity * o;
        }        
    }
    else 
    static if (summation == Summation.KB2)
    {
        F s;
        F cs;
        F ccs;
    }
    else 
    static if (summation == Summation.KBN)
    {
        F s;
        F c;
    }
    else
    static if (summation == Summation.Kahan)
    {
        F s;
        F c;
        F y; // do not declare in the loop/put (algo can be used for matrixes and etc)
        F t; // ditto
    }


public:

    ///
    this(T n)
    {
        static if (summation == Summation.Precise)
        {
            partials = scopeBuffer(scopeBufferArray);
            s = 0.0;
            o = 0;
            if (n) put(n);
        }
        else 
        static if (summation == Summation.KB2)
        {
            s = n;
            cs = 0.0;
            ccs = 0.0;

        }
        else 
        static if (summation == Summation.KBN)
        {
            s = n;
            c = 0.0;
        }
        else
        static if (summation == Summation.Kahan)
        {
            s = n;
            c = 0.0;
        }
    }

    ///
    //@disable this();

    // free ScopeBuffer
    static if (summation == Summation.Precise)
    ~this()
    {
        partials.free;
    }

    // copy ScopeBuffer if necessary
    static if (summation == Summation.Precise)
    this(this)
    {
        auto a = partials[];
        if (scopeBufferArray.ptr !is a.ptr)
        {
            partials = scopeBuffer(scopeBufferArray);
            partials.put(a);
        }
    }

    ///Adds $(D x) to the internal partial sums.
    void put(T n)
    {
        F x = n;
        static if (summation == Summation.Precise)
        {            
            if (.isFinite(x))
            {
                size_t i;
                foreach (y; partials[])
                {
                    F h = x + y;
                    if (.isInfinity(cast(T)h))
                    {
                        if (fabs(x) < fabs(y))
                        {
                            F t = x; x = y; y = t;
                        }
                        //h == -F.infinity
                        if (signbit(h))
                        {
                            x += M;
                            x += M;
                            o--;
                        }
                        //h == +F.infinity
                        else
                        {
                            x -= M;
                            x -= M;
                            o++;
                        }
                        debug(numeric) assert(x.isFinite);
                        h = x + y;
                    }
                    debug(numeric) assert(h.isFinite);
                    F l;
                    if (fabs(x) < fabs(y))
                    {
                        F t = h - y;
                        l = x - t;
                    }
                    else
                    {
                        F t = h - x;
                        l = y - t;
                    }
                    debug(numeric) assert(l.isFinite);
                    if (l)
                    {
                        partials[i++] = l;
                    }
                    x = h;
                }
                partials.length = i;
                if (x)
                {
                    partials.put(x);
                }
            }
            else
            {
                s += x;
            }
        }
        else 
        static if (summation == Summation.KB2)
        {
            static if (isFloatingPoint!F)
            {
                F t = s + x;
                F c = void;
                if (fabs(s) >= fabs(x))
                {
                    F d = s - t;
                    c = d + x;
                }
                else
                {
                    F d = x - t;
                    c = d + s;
                }
                s = t;
                t = cs + c;
                if (fabs(cs) >= fabs(c))
                {
                    F d = cs - t;
                    d += c;
                    ccs += d;
                }
                else
                {
                    F d = c - t;
                    d += cs;
                    ccs += d;
                }
                cs = t;
            }
            else
            {
                F t = s + x;
                if (fabs(s.re) < fabs(x.re))
                {
                    auto t_re = s.re;
                    s.re = x.re;
                    x.re = t_re;
                }
                if (fabs(s.im) < fabs(x.im))
                {
                    auto t_im = s.im;
                    s.im = x.im;
                    x.im = t_im;
                }
                F c = (s-t)+x;
                s = t;
                if (fabs(cs.re) < fabs(c.re))
                {
                    auto t_re = cs.re;
                    cs.re = c.re;
                    c.re = t_re;
                }
                if (fabs(cs.im) < fabs(c.im))
                {
                    auto t_im = cs.im;
                    cs.im = c.im;
                    c.im = t_im;
                }
                F d = cs - t;
                d += c;
                ccs += d;
                cs = t;
            }
        }
        else 
        static if (summation == Summation.KBN)
        {
            static if (isFloatingPoint!F)
            {
                F t = s + x;
                if (fabs(s) >= fabs(x))
                {
                    F d = s - t;
                    d += x;
                    c += d;
                }
                else
                {
                    F d = x - t;
                    d += s;
                    c += d;
                }
                s = t;
            }
            else
            {
                F t = s + x;
                if (fabs(s.re) < fabs(x.re))
                {
                    auto t_re = s.re;
                    s.re = x.re;
                    x.re = t_re;
                }
                if (fabs(s.im) < fabs(x.im))
                {
                    auto t_im = s.im;
                    s.im = x.im;
                    x.im = t_im;
                }
                F d = s - t;
                d += x;
                c += d;
                s = t;
            }
        }
        else
        static if (summation == Summation.Kahan)
        {
            y = x - c;
            t = s + y;
            c = t - s;
            c -= y;
            s = t;
        }
    }

    /++
    Adds $(D x) to the internal partial sums.
    This operation doesn't re-establish special
    value semantics across iterations (i.e. handling -inf + inf).
    Preconditions: $(D isFinite(x)).
    +/
    static if (summation == Summation.Precise)
    package void unsafePut(F x)
    in {
        assert(.isFinite(x));
    }
    body {
        size_t i;
        foreach (y; partials[])
        {
            F h = x + y;
            debug(numeric) assert(.isFinite(h));
            F l;
            if (fabs(x) < fabs(y))
            {
                F t = h - y;
                l = x - t;
            }
            else
            {
                F t = h - x;
                l = y - t;
            }
            debug(numeric) assert(.isFinite(l));
            if (l)
            {
                partials[i++] = l;
            }
            x = h;
        }
        partials.length = i;
        if (x)
        {
            partials.put(x);
        }
    }

    ///
    unittest {
        import std.math, std.algorithm, std.range;
        auto r = iota(1000).map!(n => 1.7L.pow(n+1) - 1.7L.pow(n));
        Summator!real s = 0.0;
        put(s, r);
        s -= 1.7L.pow(1000);
        assert(s.sum() == -1);
    }

    ///Returns the value of the sum.
    T sum()
    {
        /++
        Returns the value of the sum, rounded to the nearest representable
        floating-point number using the round-half-to-even rule.
        Result can be differ from the exact value on $(D X86), $(D nextDown(proir) <= result &&  result <= nextUp(proir)).
        +/
        static if (summation == Summation.Precise)
        {
            debug(numeric)
            {
                foreach (y; partials[])
                {
                    assert(y);
                    assert(y.isFinite);
                }
                //TODO: Add Non-Overlapping check to std.math
                import std.algorithm : isSorted, map;
                assert(partials[].map!(a => fabs(a)).isSorted);
            }

            if (s)
                return s;
            auto parts = partials[];
            F y = 0.0;
            //pick last
            if (parts.length)
            {
                y = parts[$-1];
                parts = parts[0..$-1];
            }
            if (o)
            {
                immutable F of = o;
                if (y && (o == -1 || o == 1)  && signbit(of * y))
                {
                    // problem case: decide whether result is representable
                    y /= 2;
                    F x = of * M;
                    immutable F h = x + y;
                    F t = h - x;
                    F l = (y - t) * 2;
                    y = h * 2;
                    if (.isInfinity(cast(T)y))
                    {
                        // overflow, except in edge case...
                        x = h + l;
                        t = x - h;
                        y = parts.length && t == l && !signbit(l*parts[$-1]) ?
                            x * 2 :
                            F.infinity * of;
                        parts = null;
                    }
                    else if (l)
                    {
                        bool _break;
                        y = partialsReducePred(y, l, parts.length ? parts[$-1] : 0, _break);
                        if (_break)
                            parts = null;
                    }
                }
                else
                {
                    y = F.infinity * of;
                    parts = null;
                }
            }
            return partialsReduce(y, parts);
        }
        else 
        static if (summation == Summation.KB2)
        {
            return cast(T)(s+(cs+ccs));
        }
        else 
        static if (summation == Summation.KBN)
        {
            return cast(T)(s+c);
        }
        else
        static if (summation == Summation.Kahan)
        {
            return cast(T)s;
        }
    }

    version(none)
    static if (summation == Summation.Precise)
    F partialsSum()
    {
        debug(numeric) partialsDebug;
        auto parts = partials[];
        F y = 0.0;
        //pick last
        if (parts.length)
        {
            y = parts[$-1];
            parts = parts[0..$-1];
        }
        return partialsReduce(y, parts);
    }

    ///Returns $(D Summator) with extended internal partial sums.
    C opCast(C : Summator!P, P)()
        if (
            isMutable!C &&
            P.max_exp >= T.max_exp &&
            P.mant_dig >= T.mant_dig
        )
    {
        static if (is(P == T))
                return this;
        else
        static if (summation == Summation.Precise)
        {
            typeof(return) ret = void;
            ret.s = s;
            ret.o = o;
            ret.partials = scopeBuffer(ret.scopeBufferArray);
            foreach (p; partials[])
            {
                ret.partials.put(p);
            }
            enum exp_diff = P.max_exp / T.max_exp;
            static if (exp_diff)
            {
                if (ret.o)
                {
                    immutable f = ret.o / exp_diff;
                    immutable t = cast(int)(ret.o % exp_diff);
                    ret.o = f;
                    ret.put((P(2) ^^ T.max_exp) * t);
                }
            }
            return ret;
        }
        else 
        static if (summation == Summation.KB2)
        {
            typeof(return) ret = void;
            ret.s = s;
            ret.cs = cs;
            ret.ccs = ccs;
            return ret;
        }
        else 
        static if (summation == Summation.KBN)
        {
            typeof(return) ret = void;
            ret.s = s;
            ret.c = c;
            return ret;
        }
        else
        static if (summation == Summation.Kahan)
        {
            typeof(return) ret = void;
            ret.s = s;
            ret.c = c;
            return ret;
        }
    }

    ///
    unittest
    {
        import std.math;
        float  M = 2.0f ^^ (float.max_exp-1);
        double N = 2.0  ^^ (float.max_exp-1);
        auto s = Summator!float(0); //float summator
        s += M;
        s += M;
        assert(float.infinity == s.sum());
        auto e = cast(Summator!double) s;
        assert(isFinite(e.sum()));
        assert(N+N == e.sum());
    }

    /++
    $(D cast(C)) operator overloading. Returns $(D cast(C)sum()).
    See also: $(D cast)
    +/
    C opCast(C)() if (is(Unqual!C == T))
    {
        return cast(C)sum();
    }

    ///Operator overloading.
    void opAssign(T rhs)
    {
        static if (summation == Summation.Precise)
        {
            partials.free;
            partials = scopeBuffer(scopeBufferArray);
            s = 0.0;
            o = 0;
            if (rhs) put(rhs);
        }
        else 
        static if (summation == Summation.KB2)
        {
            s = rhs;
            cs = 0.0;
            ccs = 0.0;
        }
        else 
        static if (summation == Summation.KBN)
        {
            s = rhs;
            c = 0.0;
        }
        else
        static if (summation == Summation.Kahan)
        {
            s = rhs;
            c = 0.0;
        }
    }

    ///ditto
    void opOpAssign(string op : "+")(T rhs)
    {
        put(rhs);
    }

    ///ditto
    void opOpAssign(string op : "+")(ref Summator rhs)
    {
        static if (summation == Summation.Precise)
        {
            s += rhs.s;
            o += rhs.o;
            foreach (f; rhs.partials[])
                put(f);
        }
        else 
        static if (summation == Summation.KB2)
        {
            put(rhs.ccs);
            put(rhs.cs);
            put(rhs.s);
        }
        else 
        static if (summation == Summation.KBN)
        {
            put(rhs.c);
            put(rhs.s);
        }
        else
        static if (summation == Summation.Kahan)
        {
            put(rhs.s);
        }
    }

    ///ditto
    void opOpAssign(string op : "-")(T rhs)
    {
        static if (summation == Summation.Kahan)
        {
            y = 0.0;
            y -= rhs;
            y -= c;
            t = s + y;
            c = t - s;
            c -= y;
            s = t;
        }
        else
        {
            put(-rhs);
        }
    }

    ///ditto
    void opOpAssign(string op : "-")(ref Summator rhs)
    {
        static if (summation == Summation.Precise)
        {
            s -= rhs.s;
            o -= rhs.o;
            foreach (f; rhs.partials[])
                put(-f);
        }
        else 
        static if (summation == Summation.KB2)
        {
            put(-rhs.ccs);
            put(-rhs.cs);
            put(-rhs.s);
        }
        else 
        static if (summation == Summation.KBN)
        {
            put(-rhs.c);
            put(-rhs.s);
        }
        else
        static if (summation == Summation.Kahan)
        {
            this -= rhs.s;
        }
    }

    ///
    nothrow unittest {
        import std.math, std.algorithm, std.range;
        auto r1 = iota(500).map!(a => 1.7L.pow(a+1) - 1.7L.pow(a));
        auto r2 = iota(500, 1000).map!(a => 1.7L.pow(a+1) - 1.7L.pow(a));
        Summator!real s1 = 0, s2 = 0.0;
        foreach (e; r1) s1 += e;
        foreach (e; r2) s2 -= e;
        s1 -= s2;
        s1 -= 1.7L.pow(1000);
        assert(s1.sum() == -1);
    }

    nothrow unittest {
        import std.typetuple;
        with(Summation)
        foreach (summation; TypeTuple!(Kahan, KBN, KB2, Precise))
        foreach (T; TypeTuple!(float, double, real))
        {
            Summator!(T, summation) sum = 1;
            sum += 3;
            assert(sum.sum == 4);
            sum -= 10;
            assert(sum.sum == -6);
            Summator!(T, summation) sum2 = 3;
            sum -= sum2;
            assert(sum.sum == -9);
            sum2 = 100;
            sum += 100;
            assert(sum.sum == 91);
            auto sum3 = cast(Summator!(real, summation))sum;
            assert(sum3.sum == 91);
            sum = sum2;
        }
    }

    static if (summation == Summation.Precise)
    {
        ///Returns $(D true) if current sum is a NaN.
        bool isNaN()
        {
            return .isNaN(s);
        }

        ///Returns $(D true) if current sum is finite (not infinite or NaN).
        bool isFinite()
        {
            if (s)
                return false;
            return !overflow;
        }

        ///Returns $(D true) if current sum is ±∞.
        bool isInfinity()
        {
            return .isInfinity(s) || overflow();
        }
    }
    else static if(isFloatingPoint!F)
    {
        ///Returns $(D true) if current sum is a NaN.
        bool isNaN()
        {
            return .isNaN(sum());
        }

        ///Returns $(D true) if current sum is finite (not infinite or NaN).
        bool isFinite()
        {
            return cast(bool).isFinite(sum());
        }

        ///Returns $(D true) if current sum is ±∞.
        bool isInfinity()
        {
            return .isInfinity(sum());
        }
    }
}

///
unittest {
    import std.range;
    import std.algorithm: swap;

    ///Moving mean
    class MovingAverage
    {
        Summator!double summator;
        double[] circularBuffer;
        size_t frontIndex;

        double avg() @property
        {
            return summator.sum() / circularBuffer.length;
        }

        this(double[] buffer)
        {
            assert(!buffer.empty);
            circularBuffer = buffer;
            summator = 0;
            .put(summator, buffer);
        }

        ///operation without rounding
        void put(double x)
        {
            summator += x;
            swap(circularBuffer[frontIndex++], x);
            summator -= x;
            frontIndex %= circularBuffer.length;
        }
    }

    /// ma always keeps pricese average of last 1000 elements
    auto ma = new MovingAverage(iota(0.0, 1000.0).array);
    assert(ma.avg == (1000*999/2) / 1000.0);
    /// move by 10 elements
    put(ma, iota(1000.0, 1010.0));
    assert(ma.avg == (1010*1009/2 - 10*9/2) / 1000.0);
}

// check default constructor
unittest
{
    Summator!double d;
    assert(d.isNaN());
    assert(d.sum().isNaN());
    d += 100;
    assert(d.isNaN());
    assert(d.sum().isNaN());
    d = 1;
    d += 1000;
    assert(d.sum == 1001);
}

nothrow unittest
{
    import std.range;
    import std.algorithm;
    import std.math;

    Summator!double summator = 0.0;

    enum double M = (cast(double)2) ^^ (double.max_exp - 1);
    Tuple!(double[], double)[] tests = [
        tuple(new double[0], 0.0),
        tuple([0.0], 0.0),
        tuple([1e100, 1.0, -1e100, 1e-100, 1e50, -1, -1e50], 1e-100),
        tuple([1e308, 1e308, -1e308], 1e308),
        tuple([-1e308, 1e308, 1e308], 1e308),
        tuple([1e308, -1e308, 1e308], 1e308),
        tuple([M, M, -2.0^^1000], 1.7976930277114552e+308),
        tuple([M, M, M, M, -M, -M, -M], 8.9884656743115795e+307),
        tuple([2.0^^53, -0.5, -2.0^^-54], 2.0^^53-1.0),
        tuple([2.0^^53, 1.0, 2.0^^-100], 2.0^^53+2.0),
        tuple([2.0^^53+10.0, 1.0, 2.0^^-100], 2.0^^53+12.0),
        tuple([2.0^^53-4.0, 0.5, 2.0^^-54], 2.0^^53-3.0),
        tuple([M-2.0^^970, -1, M], 1.7976931348623157e+308),
        tuple([double.max, double.max*2.^^-54], double.max),
        tuple([double.max, double.max*2.^^-53], double.infinity),
        tuple(iota(1, 1001).map!(a => 1.0/a).array , 7.4854708605503451),
        tuple(iota(1, 1001).map!(a => (-1.0)^^a/a).array, -0.69264743055982025), //0.693147180559945309417232121458176568075500134360255254120680...
        tuple(iota(1, 1001).map!(a => 1.0/a).retro.array , 7.4854708605503451),
        tuple(iota(1, 1001).map!(a => (-1.0)^^a/a).retro.array, -0.69264743055982025),
        tuple([double.infinity, -double.infinity, double.nan], double.nan),
        tuple([double.nan, double.infinity, -double.infinity], double.nan),
        tuple([double.infinity, double.nan, double.infinity], double.nan),
        tuple([double.infinity, double.infinity], double.infinity),
        tuple([double.infinity, -double.infinity], double.nan),
        tuple([-double.infinity, 1e308, 1e308, -double.infinity], -double.infinity),
        tuple([M-2.0^^970, 0.0, M], double.infinity),
        tuple([M-2.0^^970, 1.0, M], double.infinity),
        tuple([M, M], double.infinity),
        tuple([M, M, -1], double.infinity),
        tuple([M, M, M, M, -M, -M], double.infinity),
        tuple([M, M, M, M, -M, M], double.infinity),
        tuple([-M, -M, -M, -M], -double.infinity),
        tuple([M, M, -2.^^971], double.max),
        tuple([M, M, -2.^^970], double.infinity),
        tuple([-2.^^970, M, M, -2.^^-1074], double.max),
        tuple([M, M, -2.^^970, 2.^^-1074], double.infinity),
        tuple([-M, 2.^^971, -M], -double.max),
        tuple([-M, -M, 2.^^970], -double.infinity),
        tuple([-M, -M, 2.^^970, 2.^^-1074], -double.max),
        tuple([-2.^^-1074, -M, -M, 2.^^970], -double.infinity),
        tuple([2.^^930, -2.^^980, M, M, M, -M], 1.7976931348622137e+308),
        tuple([M, M, -1e307], 1.6976931348623159e+308),
        tuple([1e16, 1., 1e-16], 10000000000000002.0),
    ];
    foreach (i, test; tests)
    {
        foreach (t; test[0]) summator.put(t);
        auto r = test[1];
        auto s = summator.sum;
        version(X86)
        {
            assert(summator.isNaN() == r.isNaN());
            assert(summator.isFinite() == r.isFinite() || r == -double.max && s == -double.infinity || r == double.max && s == double.infinity);
            assert(summator.isInfinity() == r.isInfinity() || r == -double.max && s == -double.infinity || r == double.max && s == double.infinity);
            assert(nextDown(s) <= r && r <= nextUp(s) || s.isNaN && r.isNaN);            
        }
        else
        {
            assert(summator.isNaN() == r.isNaN());
            assert(summator.isFinite() == r.isFinite());
            assert(summator.isInfinity() == r.isInfinity());
            assert(s == r || s.isNaN && r.isNaN);                        
        }
        summator = 0.0;
    }
}


private:

template isComplex(C)
{
    import std.complex : Complex;
    enum bool isComplex = is(C : Complex!F, F);
}

version(X86)
{

}
else
{
    // FIXME (perfomance issue): fabs in std.math available only for for real.
    F fabs(F)(F f) //+-0, +-NaN, +-inf doesn't matter
    {
        if (__ctfe)
        {
            return f < 0 ? -f : f;
        }
        else
        {
            version(LDC)
            {
                import ldc.intrinsics : llvm_fabs;
                return llvm_fabs(f);
            }
            else
            {
                import core.stdc.tgmath : fabs;
                return fabs(f);
            }
        }
    }    
}

template isSummable(F)
{
    enum bool isSummable =
        __traits(compiles,
        {
            F a = 0.1, b, c;
            b = 2.3;
            c = a + b;
            c = a - b;
            a += b;
            a -= b;
        });
}

template isSummable(Range, F)
{
    enum bool isSummable =
        isInputRange!Range &&
        isImplicitlyConvertible!(Unqual!(ForeachType!Range), F) &&
        !isInfinite!Range &&
        isSummable!F;
}

/++
Naive summation algorithm.
+/
F sumNaive(Range, F = Unqual!(ForeachType!Range))(Range r, F s = 0)
{
    foreach (x; r)
    {
        s += x;
    }
    return s;
}

///TODO
alias sumFast = sumNaive;

/++
$(LUCKY Pairwise summation) algorithm. Range must be a finite sliceable range.
+/
F sumPairwise(Range, F = Unqual!(ForeachType!Range))(Range r)
    if (hasLength!Range && hasSlicing!Range)
{
    import std.range : hasLength, hasSlicing;
    static assert(hasLength!Range && hasSlicing!Range);
    switch (r.length)
    {
        case 0: return F(0);
        case 1: return cast(F)r[0];
        case 2: return cast(F)(r[0] + cast(F)r[1]);
        default: auto n = r.length/2; return cast(F)(sumPairwise!(Range, F)(r[0..n]) + sumPairwise!(Range, F)(r[n..$]));
    }
}

F sumPairwise(Range, F = Unqual!(ForeachType!Range))(Range r, F seed)
{
    F s = seed;
    s += sumPairwise!Range(r);
    return s;
}


template sumGenericKahan(Summation summation)
{
    T sumGenericKahan(Range, T = Unqual!(ForeachType!Range))(Range r, T s = 0.0)
        if (summation == Summation.Kahan && isSummable!T
        || (summation == Summation.KBN || summation == Summation.KB2) 
            && (isFloatingPoint!T || isComplex!T))
    {
        auto sum = Summator!(Unqual!T, summation)(s);
        foreach (e; r)
        {
            sum.put(e);
        }
        return sum.sum;
    }
}

/++
$(LUCKY Kahan summation) algorithm.
+/
/++
---------------------
s := x[1]
c := 0
FOR k := 2 TO n DO
y := x[k] - c
t := s + y
c := (t - s) - y
s := t
END DO
---------------------
+/
alias sumKahan = sumGenericKahan!(Summation.Kahan);


/++
$(LUCKY Kahan-Babuška-Neumaier summation algorithm).
$(D KBN) gives more accurate results then $(D Kahan).
+/
/++
---------------------
s := x[1]
c := 0
FOR i := 2 TO n DO
t := s + x[i]
IF ABS(s) >= ABS(x[i]) THEN
    c := c + ((s-t)+x[i])
ELSE
    c := c + ((x[i]-t)+s)
END IF
s := t
END DO
s := s + c
---------------------
+/
alias sumKBN = sumGenericKahan!(Summation.KBN);


/++
$(LUCKY Generalized Kahan-Babuška summation algorithm), order 2.
$(D KB2) gives more accurate results then $(D Kahan) and $(D KBN).
+/
/++
---------------------
s := 0 ; cs := 0 ; ccs := 0
FOR j := 1 TO n DO
    t := s + x[i]
    IF ABS(s) >= ABS(x[i]) THEN
        c := (s-t) + x[i]
    ELSE
        c := (x[i]-t) + s
    END IF
    s := t
    t := cs + c
    IF ABS(cs) >= ABS(c) THEN
        cc := (cs-t) + c
    ELSE
        cc := (c-t) + cs
    END IF
    cs := t
    ccs := ccs + cc
END FOR
RETURN s+cs+ccs
---------------------
+/
alias sumKB2 = sumGenericKahan!(Summation.KB2);

unittest
{
    import std.typetuple;
    with(Summation)
    foreach (F; TypeTuple!(float, double, real))
    {
        F[] ar = [1, 2, 3, 4];
        F r = 10;
        assert(r == ar.fsum!Fast());
        assert(r == ar.fsum!Pairwise());
        assert(r == ar.fsum!Kahan());
        assert(r == ar.fsum!KBN());
        assert(r == ar.fsum!KB2());
        assert(r == ar.fsum!Appropriate());
    }
}

//@@@BUG@@@: DMD 2.066 Segmentation fault (core dumped)
version(none)
unittest
{
    import core.simd;
    static if (__traits(compiles, double2.init + double2.init))
    {
        double2[] ar = [double2([1.0, 2]), double2([2, 3]), double2([3, 4]), double2([4, 6])];
        assert(ar.sumFast().array == double2([10, 14]).array);
        assert(ar.sumPairwise().array == double2([10, 14]).array);
        assert(ar.sumKahan().array == double2([10, 14]).array);
    }
}


/++
Precise summation.
+/
F sumPrecise(Range, F = Unqual!(ForeachType!Range))(Range r, F seed = 0)
    if (isFloatingPoint!F || isComplex!F)
{
    static if (isFloatingPoint!F)
    {
        auto sum = Summator!F(seed);
        foreach (e; r)
        {
            sum.put(e);
        }
        return sum.sum;
    }
    else
    {
        alias T = typeof(F.init.re);
        static if (isForwardRange!Range)
        {
            auto s = r.save;
            auto sum = Summator!T(seed.re);
            foreach (e; r)
            {
                sum.put(e.re);
            }
            T sumRe = sum.sum;
            sum = seed.im;
            foreach (e; s)
            {
                sum.put(e.im);
            }
            return F(sumRe, sum.sum);
        }
        else
        {
            auto sumRe = Summator!T(seed.re);
            auto sumIm = Summator!T(seed.im);
            foreach (e; r)
            {
                sumRe.put(e.re);
                sumIm.put(e.im);
            }
            return F(sumRe.sum, sumIm.sum);
        }
    }
}

template Algo(Range, Summation summation)
{
    static if (summation == Summation.Fast)
        alias Algo = sumFast;
    else
    static if (summation == Summation.Naive)
        alias Algo = sumNaive;
    else
    static if (summation == Summation.Pairwise)
        alias Algo = sumPairwise;
    else
    static if (summation == Summation.Kahan)
        alias Algo = sumKahan;
    else
    static if (summation == Summation.KBN)
        alias Algo = sumKBN;
    else
    static if (summation == Summation.KB2)
        alias Algo = sumKB2;
    else
    static if (summation == Summation.Precise)
        alias Algo = sumPrecise;
    else
    static if (summation == Summation.Appropriate)
    {
        static if (isRandomAccessRange!Range)
            alias Algo = Algo!(Range, Summation.Pairwise);
        else
        static if (isFloatingPoint!(ForeachType!Range) || isComplex!(ForeachType!Range))
            alias Algo = Algo!(Range, Summation.KBN);
        else
            alias Algo = Algo!(Range, Summation.Kahan);
    }
    else
    static assert(0);
}

unittest {
    static struct UDT
    {
        UDT opBinary(string op)(auto ref UDT rhs) const {
            UDT ret;
            return ret;
        }
        UDT opOpAssign(string op)(auto ref UDT rhs) {
            return this;
        }
        this(double f){}
        void opAssign(double f){}
    }
    import std.range;
    import std.complex : Complex;
    static assert(__traits(isSame, Algo!(double[], Summation.Appropriate), sumPairwise));
    static assert(__traits(isSame, Algo!(Complex!double[], Summation.Appropriate), sumPairwise));
    static assert(__traits(isSame, Algo!(UDT[], Summation.Appropriate), sumPairwise));
    static assert(__traits(isSame, Algo!(ForwardRange!double, Summation.Appropriate), sumKBN));
    static assert(__traits(isSame, Algo!(ForwardRange!(Complex!double), Summation.Appropriate), sumKBN));
    static assert(__traits(isSame, Algo!(ForwardRange!UDT, Summation.Appropriate), sumKahan));
}
