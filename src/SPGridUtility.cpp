#include "SPGridUtility.h"
#include <cstdio>
#include <stdexcept>
#include <immintrin.h>

#ifdef _WIN32
#include <windows.h>
#include <Psapi.h>
#include <tchar.h>
#include <intrin.h>
#endif

#ifdef __linux__
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <string>
#include <cstdlib>
#include <cstring>
#endif

namespace SPGrid {

#if 1

	uint64_t bit_spread(uint32_t bits, uint64_t mask) noexcept {
		return _pdep_u64(static_cast<uint64_t>(bits), mask);
	}
	uint32_t bit_pack(uint64_t bits, uint64_t mask) noexcept {
		return static_cast<uint32_t>(_pext_u64(bits, mask));
	}

	uint64_t index_to_morton_2d_64(uint32_t x, uint32_t y) noexcept {
		return _pdep_u64(x, 0x5555555555555555) | _pdep_u64(y, 0xaaaaaaaaaaaaaaaa);
	}

	void morton_to_index_2d_64(uint64_t m, uint32_t& x, uint32_t& y) noexcept {
		x = _pext_u64(m, 0x5555555555555555);
		y = _pext_u64(m, 0xaaaaaaaaaaaaaaaa);
	}

	uint32_t index_to_morton_2d_32(uint32_t x, uint32_t y) noexcept {
		return _pdep_u32(x, 0x55555555) | _pdep_u32(y, 0xaaaaaaaa);
	}

	void morton_to_index_2d_32(uint32_t m, uint32_t& x, uint32_t& y) noexcept {
		x = _pext_u32(m, 0x55555555);
		y = _pext_u32(m, 0xaaaaaaaa);
	}

#else
	uint64_t bit_spread(uint32_t bits, uint64_t mask) noexcept {
		// gcc: __builtin_ctz
		// msvd: __tzcnt64
		uint64_t rmask = reverse_bits(mask);	//_byteswap_uint64 doesn't apply
		uint64_t result = 0;
		unsigned char lz = 0;
		while (rmask) {
			#ifdef _WIN32
			lz = __lzcnt64(rmask) + 1;
			#endif
			#ifdef __linux__
			lz = __builtin_clzl(rmask) + 1;
			#endif
			result = result << lz | (bits & 1);
			bits >>= 1, rmask <<= lz;
		}
		#ifdef _WIN32
		result = reverse_bits(result) >> __lzcnt64(mask);
		#endif
		#ifdef __linux__
		result = reverse_bits(result) >> __builtin_clzl(mask);
		#endif
		return result;
	}

	uint32_t bit_pack(uint64_t bits, uint64_t mask) noexcept {
		uint64_t ulresult = 0;
		uint64_t uldata = bits;
		int count = 0;
		ulresult = 0;

		uint64_t rmask = reverse_bits(mask);
		unsigned char lz = 0;

		while (rmask) {
			#ifdef _WIN32
			lz = __lzcnt64(rmask);
			#endif
			#ifdef __linux__
			lz = __builtin_clzl(rmask);
			#endif
			uldata >>= lz;
			ulresult = ulresult << 1 | (uldata & 1);
			uldata >>= 1;
			rmask <<= lz + 1;
			count++;
		}
		ulresult <<= 64 - count; // 64 means 64 bits ... maybe not use a constant 64 ...
		ulresult = reverse_bits(ulresult);
		return static_cast<uint32_t>(ulresult);
	}
#endif

	void checkTotalMemory() {
#ifdef _WIN32
		static MEMORYSTATUSEX memInfo;
		memInfo.dwLength = sizeof(MEMORYSTATUSEX);
		GlobalMemoryStatusEx(&memInfo);
		DWORDLONG totalVirtualMem = memInfo.ullTotalPageFile;
		DWORDLONG totalPhysMem = memInfo.ullTotalPhys;
		printf("Total Virtual Mem: %llu (KB)\tTotal Physical Mem: %llu (KB)\n", totalVirtualMem >> 10, totalPhysMem >> 10);
#endif

#ifdef __linux__
		/// total virtual memory
		struct sysinfo memInfo;
		sysinfo(&memInfo);
		int64_t totalVirtualMem = memInfo.totalram;
		//Add other values in next statement to avoid int overflow on right hand side...
		totalVirtualMem += memInfo.totalswap;
		totalVirtualMem *= memInfo.mem_unit;
		/// total physical memory
		int64_t totalPhysMem = memInfo.totalram;
		//Multiply in next statement to avoid int overflow on right hand side...
		totalPhysMem *= memInfo.mem_unit;
		printf("Total Virtual Mem: %llu (KB)\tTotal Physical Mem: %llu (KB)\n", totalVirtualMem >> 10, totalPhysMem >> 10);
#endif
	}

#ifdef __linux__

	int parseLine(char* line){
	    // This assumes that a digit will be found and the line ends in " Kb".
	    int i = strlen(line);
	    const char* p = line;
	    while (*p <'0' || *p > '9') p++;
	    line[i-3] = '\0';
	    i = atoi(p);
	    return i;
	}

	int getProcessInfo(std::string&& key) { //Note: this value is in KB!
		FILE* file = fopen("/proc/self/status", "r");
		int result = -1;
		char line[128];

		while (fgets(line, 128, file) != NULL) {
			if (strncmp(line, key.data(), key.size()) == 0) {
				result = parseLine(line);
				break;
			}
		}
		fclose(file);
		return result;
	}
#endif

	void checkProcessMemory() {
#ifdef _WIN32
		static PROCESS_MEMORY_COUNTERS_EX pmc;
		GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
		SIZE_T virtualMemUsed = pmc.PrivateUsage;
		SIZE_T physMemUsed = pmc.WorkingSetSize;
		printf("Used Virtual Mem: %llu (KB)\tUsed Physical Mem: %llu (KB)\n", virtualMemUsed >> 10, physMemUsed >> 10);
#endif

#ifdef __linux__
		printf("Used Virtual Mem: %llu (KB)\tUsed Physical Mem: %llu (KB)\n", getProcessInfo("VmSize:"), getProcessInfo("VmRSS:"));
#endif
	}
	bool checkAddressResident(void * const addr) {
#ifdef _WIN32
		void* page_addr = addr;
		static MEMORY_BASIC_INFORMATION meminfo;
		VirtualQuery(page_addr, &meminfo, sizeof(meminfo));
		/*
		printf("address [%llx] state %llx\n", page_addr, meminfo.State);
		if (meminfo.State != MEM_COMMIT)
			printf("page %llx hasn't been committed\n", page_addr);
		printf("base address of %llx: %llx\n", addr, meminfo.BaseAddress);
		*/
		return meminfo.State == MEM_COMMIT;
#endif

#ifdef __linux__
		void* page_addr = reinterpret_cast<void*>(reinterpret_cast<uint64_t>(addr) & 0xfffffffffffff000ULL);
		unsigned char status;
		if (mincore(page_addr, 4096, &status) == -1)
			switch (errno) {
			case ENOMEM: throw std::runtime_error("In Check_Address_Resident() : Input address hasn\'t been mapped\n");
			default: throw std::runtime_error(" In Check_Address_Resident() : mincore() failed\n");
			}
		if (!status) throw std::runtime_error("In Check_Address_Resident() : Input address not resident in memory\n");
		return true;
#endif
	}

}