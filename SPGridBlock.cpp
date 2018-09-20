#include "SPGridBlock.h"
#include "SPGridUtility.h"

namespace SPGrid {

	SPGridBlock<3>::SPGridBlock(const uint32_t xsize_input, const uint32_t ysize_input, const uint32_t zsize_input,
		const uint32_t block_xbits, const uint32_t block_ybits, const uint32_t block_zbits)
		: xsize(xsize_input), ysize(ysize_input), zsize(zsize_input),
		/// critical icc compilation error!!!
		//depth(nextPowerOfTwo(maximum(xsize, ysize, zsize, 1u << block_zbits))),
		//height(maximum(nextPowerOfTwo(maximum(ysize, xsize)), depth >> 2, 1u << block_ybits)),
		//width(maximum(nextPowerOfTwo(xsize), depth >> 1, 1u << block_xbits)),
		blockWidth(1 << block_xbits), blockHeight(1 << block_ybits), blockDepth(1 << block_zbits)
	{
		depth = (nextPowerOfTwo(maximum(xsize, ysize, zsize, 1u << block_zbits)));
		height = (maximum(nextPowerOfTwo(maximum(ysize, xsize)), depth >> 2, 1u << block_ybits));
		width = (maximum(nextPowerOfTwo(xsize), depth >> 1, 1u << block_xbits));
	}

	void SPGridBlock<3>::checkBounds(const int i, const int j, const int k) const {
		if (i >= width || j >= height || k >= depth)
			throw std::runtime_error(std::string("array indices out of bound"));
	}

}