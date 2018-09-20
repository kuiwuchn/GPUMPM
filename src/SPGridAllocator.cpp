#include <cstdio>
#include "SPGridAllocator.h"
#include "SPGridUtility.h"

#ifdef _WIN32
#include <windows.h>
#include <Psapi.h>
#include <tchar.h>
#endif

#ifdef __linux__
#include <sys/mman.h>
#include <unistd.h>
#endif


namespace SPGrid {

	int	SPGridAllocator::_dwPageSize = 0;
	int	SPGridAllocator::_commitCount = 0;

	void SPGridAllocator::init() {
#ifdef _WIN32
		SYSTEM_INFO sysinfo;
		GetSystemInfo(&sysinfo);
		_dwPageSize = sysinfo.dwPageSize;
#endif

#ifdef __linux__
		_dwPageSize = sysconf(_SC_PAGESIZE);
#endif
		checkTotalMemory();
		checkProcessMemory();
	}
	void* SPGridAllocator::allocate(const uint32_t x, const uint32_t y, const uint32_t z, const size_t bytesPerElement) {
#ifdef _WIN32
		auto ptr = VirtualAlloc(nullptr, x * y * z * bytesPerElement, MEM_RESERVE, PAGE_READWRITE);
#endif

#ifdef __linux__
		auto ptr = mmap(nullptr, x * y * z * bytesPerElement, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
#endif
		printf("address is: %llx\n\n", ptr);
		checkProcessMemory();
		return ptr;
	}

	void SPGridAllocator::deallocate(void* addr, size_t size) {
#ifdef _WIN32
		VirtualFree(addr, 0, MEM_RELEASE);
#endif

#ifdef __linux__
		munmap(addr, size);
#endif
		checkProcessMemory();
	}

	void SPGridAllocator::loadPage(void* ptr) {
		/// no need pre-check residence
		//if (!checkAddressResident(ptr))
		//	_commitCount++;
#ifdef _WIN32
		VirtualAllocEx(GetCurrentProcess(), ptr, _dwPageSize, MEM_COMMIT, PAGE_READWRITE);
#endif
	}

	void SPGridAllocator::freePage(void* ptr) {
#ifdef _WIN32
		VirtualFreeEx(GetCurrentProcess(), ptr, _dwPageSize, MEM_DECOMMIT);
#endif
	}


}