#ifndef __SPGRID_ITERATOR_H_
#define __SPGRID_ITERATOR_H_

#include <stdint.h>

namespace SPGrid {

	template<class MASK>
	class SPGridIterator {
	public:
		using Mask = MASK;
		const uint64_t* const blockOffsets;
		const uint32_t size;
		uint32_t blockIndex, cellIndex;
		enum { 
			cellsPerBlock = Mask::blockSize,
			cellIndexMask = cellsPerBlock - 1
		};

		SPGridIterator(const uint64_t* const blockOffsetsInput, const uint32_t sizeInput)
			: blockOffsets(blockOffsetsInput), size(sizeInput), blockIndex(0), cellIndex(0) {}

		SPGridIterator(const std::pair<const uint64_t*, uint32_t>& blocks)
			: blockOffsets(blocks.first), size(blocks.second), blockIndex(0), cellIndex(0) {}

		inline bool valid() const noexcept { return blockIndex < size; }
		inline void next() noexcept { cellIndex = (cellIndex + 1) & cellIndexMask; if (!cellIndex) blockIndex++; }
		inline void nextBlock() noexcept { cellIndex = 0; blockIndex++; }
		inline uint64_t offset() const {
			return blockOffsets[blockIndex] + static_cast<uint64_t>(cellIndex * Mask::fieldSize);
		}
		auto index() const {	///< this operation is expensive
			return Mask::offset_to_coord(blockOffsets[blockIndex] + static_cast<uint64_t>(cellIndex * Mask::fieldSize));
		}
		template<class FieldMask> inline uint64_t offset() const {
			return blockOffsets[blockIndex] + static_cast<uint64_t>(cellIndex * FieldMask::fieldSize);
		}

		/// acquire current cell, array parameter version
		template<class T_ARRAY> inline
		typename T_ARRAY::Data& data(T_ARRAY& array) const {
			using T = typename T_ARRAY::Data;
			uint64_t blockAddr = reinterpret_cast<uint64_t>(array.getPtr()) + blockOffsets[blockIndex];
			return reinterpret_cast<T*>(blockAddr)[cellIndex];
		}
		/// acquire current cell, pointer parameter version
		template<class T_ARRAY> inline
		typename T_ARRAY::Data& data(void* ptr) const {
			using T = typename T_ARRAY::Data;
			uint64_t blockAddr = reinterpret_cast<uint64_t>(ptr) + blockOffsets[blockIndex];
			return reinterpret_cast<T*>(blockAddr)[cellIndex];
		}
		/// acquire the cell shifted from the current one, array parameter version
		template<class T_ARRAY> inline
		typename T_ARRAY::Data& data(T_ARRAY& array, uint64_t offsetInput) const {
			using T = typename T_ARRAY::Data;
			uint64_t offset = T_ARRAY::Mask::packed_add(blockOffsets[blockIndex] + cellIndex * sizeof(T), offsetInput);
			uint64_t dataAddr = reinterpret_cast<uint64_t>(array.getPtr()) + offset;
			return *reinterpret_cast<T*>(dataAddr);
		}
		/// acquire the cell shifted from the current one, pointer parameter version
		template<class T_ARRAY> inline
		typename T_ARRAY::Data& data(void* dataPtr, uint64_t offsetInput) const {
			using T = typename T_ARRAY::Data;
			uint64_t offset = T_ARRAY::Mask::packed_add(blockOffsets[blockIndex] + cellIndex * sizeof(T), offsetInput);
			uint64_t dataAddr = reinterpret_cast<uint64_t>(dataPtr) + offset;
			return *reinterpret_cast<T*>(dataAddr);
		}
	};

}

#endif