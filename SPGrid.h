#ifndef __SPGRID_H_
#define __SPGRID_H_

#include "SPGridAllocator.h"
#include "SPGridMask.h"
#include "SPGridBlock.h"
#include "SPGridArray.h"
#include <stdint.h>
#include <algorithm>

#ifdef _WIN32
#include <windows.h>
#endif

#include "Setting.h"

namespace SPGrid {

	/**
	*	\class SparseGrid
	*	\brief Interface of SPGrid. Also manages the lifetime of the structure.
	*/
	template<typename Channels, int DIM>
	class SparseGrid {
	public:
		using Mask = SPGridMask<Channels, DIM>;
		enum : unsigned char { Dim = DIM };
		template<typename T_FIELD = Channels> struct Array {
			using Mask = SPGridFieldMask<Channels, T_FIELD, Dim>;
			using ArrayType = SPGridArray<T_FIELD, Mask>;
		};
		/**
		 *	\fn Array<T_FIELD>::ArrayType getArray(T_FIELD T::* field)
		 *	\brief returns an array that give full access to a certain channel of the grid
		 *	\param field the specific channel
		 */
		template<typename T, typename T_FIELD>
		typename std::enable_if<std::is_same<T, Channels>::value, typename Array<T_FIELD>::ArrayType>::type
			getArray(T_FIELD T::* field) {
			using ArrayType = typename Array<T_FIELD>::ArrayType;
			uint64_t offset = OffsetOfMember(field) << Mask::blockBitCount;
			auto ptr = reinterpret_cast<void*>(reinterpret_cast<uint64_t>(_virtualMem) + offset);
			return ArrayType(ptr, _domain);
		}
		/**
		 *	\fn Array<T_FIELD>::ArrayType getConstArray(T_FIELD T::* field) const
		 *	\brief returns an array that give read-only access to a certain channel of the grid
		 *	\param field the specific channel
		 */
		template<typename T, typename T_FIELD>
		typename std::enable_if<std::is_same<T, Channels>::value, typename Array<const T_FIELD>::ArrayType>::type
			getConstArray(T_FIELD T::* field) const {
			using ArrayType = typename Array<const T_FIELD>::ArrayType;
			uint64_t offset = OffsetOfMember(field) << Mask::blockBitCount;	///< modify field bits
			auto ptr = reinterpret_cast<void*>(reinterpret_cast<uint64_t>(_virtualMem) + offset);
			return ArrayType(ptr, _domain);
		}
		/**
		 *	Constructor.
		 *	Check the compliance of the system and reserve a region of virtual memory for the grid
		 */
		SparseGrid(uint32_t x, uint32_t y, uint32_t z);
		/**
		 *	Destructor.
		 *	Release all the memory associated with the grid
		 */
		~SparseGrid();
		auto getPtr() { return _virtualMem; }

	private:
		SPGridBlock<Dim>	_domain;	///< specify the boundary of the grid
		void*	_virtualMem{ nullptr };	///< the virtual memory base pointer for the grid
	};

	template<typename Channels, int Dim>
	SparseGrid<Channels, Dim>::SparseGrid(uint32_t x, uint32_t y, uint32_t z) :
		_domain(x, y, z, Mask::blockXBitCount, Mask::blockYBitCount, Mask::blockZBitCount) {
		SPGridAllocator::init();

#ifdef _WIN32
		/// check compliance
		SYSTEM_INFO sysinfo;
		GetSystemInfo(&sysinfo);
		auto allocGranularity = sysinfo.dwAllocationGranularity;
		auto dwPageSize = sysinfo.dwPageSize;
		printf("page size is: %d alloc granularity is: %d\n", dwPageSize, allocGranularity);
#endif
		if (SPGridAllocator::getPageSize() != 4096) throw std::runtime_error("page size different than 4KB");
		if (sizeof(void*) != 8) throw std::runtime_error("void* is not 64-bit long");

		_virtualMem = SPGridAllocator::allocate(_domain.width, _domain.height, _domain.depth, sizeof(Channels));
	}
	template<typename Channels, int Dim>
	SparseGrid<Channels, Dim>::~SparseGrid() {
#ifdef _WIN32
		SPGridAllocator::deallocate(_virtualMem);
#endif

#ifdef __linux__
		SPGridAllocator::deallocate(_virtualMem, _domain.gridVolume() * sizeof(Channels));
#endif
	}

}

#endif