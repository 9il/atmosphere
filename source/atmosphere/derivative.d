/++
Derivatives of probability density functions
+/
/**
Authors: [Ilya Yaroshenko](http://9il.github.io)

Copyright: © 2014-2015 [Ilya Yaroshenko](http://9il.github.io)

License: MIT
*/
module atmosphere.derivative;

import atmosphere.utilities;

/++
Derivative of probability density function interface
+/
interface Derivative(T)
{
	/**
	Call operator
	*/
	T opCall(T x);
}

///
alias toDerivative = convertTo!Derivative;
