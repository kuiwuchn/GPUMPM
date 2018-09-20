#ifndef __SPGRID_BLOCK_H_
#define __SPGRID_BLOCK_H_

#include <stdint.h>
#include <array>

namespace SPGrid {

	template<int DIM>
	struct SPGridBlock;

	/**
	 *	\class SPGridBlock
	 *	\brief the 3D domain of the grid
	 */
	template<>
	struct SPGridBlock<3> {
		uint32_t width, height, depth;				///< actual range of each dimension, adjusted for alignment (next power of 2) and uniformity
		uint32_t xsize, ysize, zsize;					///< domain of each dimension (x, y, z)
		uint32_t blockWidth, blockHeight, blockDepth;	///< domain of each block

		/**
		 *	Constructor.
		 *	determine the allocated grid domain and the associated block granularity
		 */
		SPGridBlock(const uint32_t xsize_input, const uint32_t ysize_input, const uint32_t zsize_input,
			const uint32_t block_xbits, const uint32_t block_ybits, const uint32_t block_zbits);

		/**
		 *	\fn uint64_t gridVolume()
		 *	\brief compute the volume of the entire grid
		 */
		uint64_t gridVolume() const {
			return static_cast<uint64_t>(width)*static_cast<uint64_t>(height)*static_cast<uint64_t>(depth);
		}
		/**
		 *	\fn uint32_t blockVolume()
		 *	\brief compute the volume of each block
		 */
		uint32_t blockVolume() const {
			return blockWidth * blockHeight * blockDepth;
		}
		/**
		 *	\fn std::array<uint32_t, 3> blockVolume()
		 *	\brief retrive the domain of each block
		 */
		auto blockSize() const {
			return std::array<uint32_t, 3>{blockWidth, blockHeight, blockDepth};
		}
		/**
		 *	\fn void checkBounds(const int i, const int j, const int k)
		 *	\brief check if the coordinate (i, j, k) is in the domain of the grid
		 */
		void checkBounds(const int i, const int j, const int k) const;
	};

}

#endif