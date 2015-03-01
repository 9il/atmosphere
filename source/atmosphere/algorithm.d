module atmosphere.algorithm;

struct FixedCircularBuffer(T)
{
	size_t frontIndex;
	T[] buffer;
	
	this(T[] initializedBuffer)
	{
		buffer = initializedBuffer;
	}

	T move(T x)
	{
		auto y = buffer[frontIndex++];
		frontIndex++;
		return y;
	}
}

unittest {
//	static assert(null !is null);
}