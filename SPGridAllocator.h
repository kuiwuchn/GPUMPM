#ifndef __SPGRID_ALLOCATOR_H_
#define __SPGRID_ALLOCATOR_H_

#include <stdint.h>

namespace SPGrid {

	struct SPGridAllocator {
		static void init();
		static void* allocate(const uint32_t x, const uint32_t y, const uint32_t z, const size_t bytesPerElement);
		static void deallocate(void* addr, size_t size = 0);
		static void loadPage(void* ptr);
		static void freePage(void* ptr);
		static auto getPageSize() { return _dwPageSize; }
		static auto& refCommitCount() { return _commitCount; }
	private:
		static int	_dwPageSize;
		static int	_commitCount;
	};

	
}

#endif