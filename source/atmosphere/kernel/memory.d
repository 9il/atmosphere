/**
Code memory managment.
*/
module atmosphere.kernel.memory;

import std.traits;
import core.stdc.string : memmove;

import matrix;



/**
Struct that represent flat matrix like in CBLAS API.
*/
struct MatrixAllocator(F)
{

	/**
	*/
	Matrix!F _matrix;

	///
	Matrix!F matrix;


	///
	this(size_t maxHeight, size_t maxWidth, size_t height)
	{
		_matrix = Matrix!F(maxHeight, maxWidth);
		_matrix.width = _matrix.shift;
		matrix.ptr = _matrix.ptr;
		matrix.shift = _matrix.shift;
		matrix.height = height;
	}



	///
	void popFrontN(size_t n)
	in {
		assert(n <= matrix.width, "n > matrix.width");
	}
	body {
		if(n < matrix.width)
		{
			matrix.width -= n;
			matrix.ptr += n;
		}
		else
		{ 
			reset;
		}
	}


	///	
	void popFront()
	{
		popFrontN(1);
	}



	///
	void reset()
	{
		matrix.ptr = _matrix.ptr;
		matrix.width = 0;
	}



	///
	void putBackN(size_t n)
	in{
		assert(matrix.shift >= matrix.width+n);
	}
	body {
		if(n > _matrix.ptrEnd-matrix.ptrEnd)
		{
			bringToFront();
		}
		matrix.width += n;
	}



	///
	void putBack()
	{
		putBackN(1);
	}



	///
	void bringToFront()
	{
		if(matrix.width)
		{
			memmove(_matrix.ptr, matrix.ptr, (matrix.shift*matrix.height)*F.sizeof);					
		}
		matrix.ptr = _matrix.ptr;
	}
}