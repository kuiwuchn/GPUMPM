#ifndef __SPGRID_SET_H_
#define __SPGRID_SET_H_

#include <stdint.h>
#include <vector>
#include <set>
#include <algorithm>
#include "SPGrid.h"

namespace SPGrid {

	/**
	 *	\class SPGridSet
	 *	\brief record the status of currently activated sparse blocks
	 */
	template<typename ARRAY>
	class SPGridSet {
	private:
		ARRAY		_flags;					///< flag channel of each cell,indicates the 
#ifdef SPGRID_SET_FEATURE_PAGE_ORDER
		std::vector<uint64_t>	_testBlocks;
		std::set<uint64_t>		_blockRecords;	
#endif
		std::vector<uint64_t>	_pageMasks;	///< array of page mask. 1 page per entry, 64 blocks per page
		std::vector<uint64_t>	_blockOffsets;	///< store all the linear offsets of the activated blocks lexicographically
		uint64_t	_maxLinearOffset;		///< the maximum linear offset of the grid
		bool		_bDirty;				///< indicate if any mask has been touched

	public:
		enum { Dim = ARRAY::Dim };
		using T = typename ARRAY::Data;	///< data type of each channel
		using Mask = typename ARRAY::Mask;	///< field mask
		using Array = decltype(_flags);
		//using Array = ARRAY;

		/**
		 *	Constructor. 
		 */
		SPGridSet(ARRAY arr) : _flags(arr) {
			auto cellCount = _flags.getDomainRef().gridVolume();
			auto blockCount = cellCount >> Mask::blockBitCount;
			_pageMasks.resize((blockCount + 0x3full) >> 6);	///< 64 blocks per entry
			std::fill(_pageMasks.begin(), _pageMasks.end(), 0);

			_maxLinearOffset = static_cast<uint64_t>(cellCount) * Mask::bytesPerElement();	///< the max offset corresponds to the maximum corner cell
			_bDirty = false;
		}

		inline const auto& getBlocks() const noexcept { return _blockOffsets; }
		inline auto& refArray() noexcept { return _flags; }
		inline void* getPtr() noexcept { return _flags.getPtr(); }
		inline const auto& getDomainRef() const noexcept { return _flags.getDomainRef(); }

		/**
		 *	\fn bool isInBound(uint64_t&& linearOffset)
		 *	\brief check if the linear offset is inside the domain
		 */
		inline bool isInBound(uint64_t&& linearOffset) const noexcept { return linearOffset < _maxLinearOffset; }
		inline bool isInBound(const uint64_t& linearOffset) const noexcept { return linearOffset < _maxLinearOffset; }
		/**
		 *	\fn void markPageActive(uint64_t linearOffset)
		 *	\brief check if the linear offset is inside the domain
		 */
		void markPageActive(uint64_t linearOffset) {
			if (isInBound(linearOffset)) {
				uint64_t pageMask = 1ULL << (linearOffset >> 12 & 0x3f);
				_pageMasks[linearOffset >> 18] |= pageMask;
				_bDirty = true;
			}
			else
				throw std::runtime_error("MarkPageActive:\toffset not in bound");
		}
		/**
		 *	\fn void mask(uint64_t linearOffset, const T m, bool commit = false)
		 *	\brief mark the cell flag with mask m, commmit page if commit if true
		 */
		void mask(uint64_t linearOffset, const T m, bool commit = false) {
			if (isInBound(linearOffset)) {
				uint64_t	pageMask = 1ULL << ((linearOffset >> 12) & 0x3f);
				_pageMasks[linearOffset >> 18] |= pageMask;
				if (commit) {
					auto ptr = reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_flags.getPtr()) + linearOffset);
					SPGridAllocator::loadPage(ptr);
				}
				_flags(linearOffset) |= m;
				_bDirty = true;
			}
			else {
				throw std::runtime_error("Mask:\toffset not in bound");
			}
		}
		/**
		 *	\fn void mask(std::array<int, Dim>	coord, const T m, bool commit = false)
		 *	\brief mark the cell flag with mask m, commmit page if commit if true
		 */
		void mask(std::array<int, Dim>	coord, const T m, bool commit = false) {
			uint64_t	linearOffset = Mask::linear_offset(coord);
			mask(linearOffset, m, commit);
		}
#ifdef SPGRID_SET_FEATURE_PAGE_ORDER
		void testAddBlock(std::array<int, 3> coord) {
			uint64_t offset = Mask::linear_offset(coord) & 0xfffffffffffff000;
			if (_blockRecords.find(offset) == _blockRecords.end()) {
				_blockRecords.insert(offset);
				_testBlocks.push_back(offset);
			}
		}
		auto& testGetBlock() {
			return _testBlocks;
		}
#endif

		/// query 
		bool isMaskMatch(uint64_t linearOffset, const T mask) const {
			if (isInBound(linearOffset))
				return mask & _flags.get(linearOffset);
			else
				return false;
		}
		bool isSet(uint64_t linearOffset, const T mask) {
			if (isInBound(linearOffset)) {
				if (_pageMasks[linearOffset >> 18] & (1ULL << (linearOffset >> 12 & 0x3f)))
					return mask & _flags(linearOffset);
			}
			return false;
		}
		bool isSet(std::array<int32_t, Dim> coord, const T mask) {
			uint64_t	linearOffset = Mask::linear_offset(coord);
			return isSet(linearOffset, mask);
		}
		bool isActive(uint64_t linearOffset) const {
			if (isInBound(linearOffset)) {
				if (_pageMasks[linearOffset >> 18] & (1ULL << (linearOffset >> 12 & 0x3f)))
					return true;
			}
			return false;
		}
		bool isActive(std::array<int32_t, Dim> coord) const {
			uint64_t	linearOffset = Mask::linear_offset(coord);
			return isActive(linearOffset);
		}

		/// mask to offset
		auto generateBlockOffsets() {
			std::vector<uint64_t>	blockOffsets;
			for (uint64_t i = 0; i < _pageMasks.size(); i++)
				if (_pageMasks[i])
					for (uint64_t pos = 0; pos < 64; pos++)
						if (_pageMasks[i] & (1ULL << pos))
							blockOffsets.push_back((i << 18ULL) | (pos << 12ULL));
			return blockOffsets;
		}
		void refreshBlockOffsets() {
			if (_bDirty)	_blockOffsets = std::move(generateBlockOffsets());
			_bDirty = false;
		}
		/// offset to mask
		void setFlags(const std::vector<uint64_t>& blockOffsets) {
			_bDirty = true;
			for (auto offset : blockOffsets)
				_pageMasks[offset >> 18] |= (1ULL << (offset >> 12 & 0x3f));
		}

		void clearMasks() {
			std::fill(_pageMasks.begin(), _pageMasks.end(), 0);
			_bDirty = true;
		}
		void clearBlocks() {
			_blockOffsets.clear();
			_bDirty = true;
		}
		void clear() { clearMasks(); clearBlocks(); }

	};

}

#endif