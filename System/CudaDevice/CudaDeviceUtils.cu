#include "CudaDeviceUtils.cuh"

namespace mn {

    __device__ bool atomicMinf(float* address, float val) {
		int* address_as_i = (int*)address;
		int old = *address_as_i, assumed;
		if (*address <= val) return false;
		do {
			assumed = old;
			old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fminf(val, __int_as_float(assumed))));
		} while (assumed != old);
		return true;
	}

	__device__ bool atomicMaxf(float* address, float val) {
		int* address_as_i = (int*)address;
		int old = *address_as_i, assumed;
		if (*address >= val) return false;
		do {
			assumed = old;
			old = ::atomicCAS(address_as_i, assumed,
				__float_as_int(::fmaxf(val, __int_as_float(assumed))));
		} while (assumed != old);
		return true;
	}

	__device__ bool atomicMinD(double* address, double val) {
#ifdef _WIN32
		uint64_t* address_as_ull = (uint64_t*)address;
		uint64_t old = *address_as_ull, assumed;
#else
		unsigned long long * address_as_ull = (unsigned long long*)address;
		unsigned long long old = *address_as_ull, assumed;
#endif // _WIN32
		if (*address <= val) return false;
		do {
			assumed = old;
			old = ::atomicCAS(address_as_ull, assumed,
				__double_as_longlong(::fmin(val, __longlong_as_double(assumed))));
		} while (assumed != old);
		return true;
	}

	__device__ bool atomicMaxD(double* address, double val) {
#ifdef _WIN32
		uint64_t* address_as_ull = (uint64_t*)address;
		uint64_t old = *address_as_ull, assumed;
#else
		unsigned long long * address_as_ull = (unsigned long long*)address;
		unsigned long long old = *address_as_ull, assumed;
#endif // _WIN32
		if (*address >= val) return false;
		do {
			assumed = old;
			old = ::atomicCAS(address_as_ull, assumed,
				__double_as_longlong(::fmax(val, __longlong_as_double(assumed))));
		} while (assumed != old);
		return true;
	}

    __device__ uint64_t Packed_Add(const uint64_t* masks,const uint64_t i,const uint64_t j) {
        uint64_t x_result=( (i | ~masks[0]) + (j & masks[0]) ) & masks[0];
        uint64_t y_result=( (i | ~masks[1]) + (j & masks[1]) ) & masks[1];
        uint64_t z_result=( (i | ~masks[2]) + (j & masks[2]) ) & masks[2];
        uint64_t w_result=( (i | masks[0]|masks[1]|masks[2]) + (j & ~(masks[0]|masks[1]|masks[2])) ) & ~(masks[0]|masks[1]|masks[2]);
        uint64_t result=x_result | y_result | z_result | w_result;
        return result;
    }

    __device__ uint64_t Bit_Spread_Mine(const uint64_t mask, const int data) {
        uint64_t rmask = __brevll(mask);
        int dat = data;
        uint64_t result = 0;
        unsigned char lz, offset = __clzll(rmask);
        while (rmask) {
            lz = __clzll(rmask) + 1;
            result = result << lz | (dat & 1);
            dat >>= 1, rmask <<= lz;
        }
        result = __brevll(result) >> __clzll(mask);
        return result;
    }

    __device__ int Bit_Pack_Mine(const uint64_t mask, const uint64_t data) {
        union{ uint64_t slresult; uint64_t ulresult; };
		uint64_t uldata=data; int count=0; ulresult=0;

		uint64_t rmask = __brevll(mask);
        unsigned char lz;

        while(rmask) {
            lz = __clzll(rmask) ;
            uldata>>=lz;
            ulresult<<=1;
            count++;
            ulresult |= (uldata & 1);
            uldata>>=1;
            rmask <<= lz+1;
        }
        ulresult<<=64-count; // 64 means 64 bits ... maybe not use a constant 64 ...
        ulresult=__brevll(ulresult);
        return (int)slresult;
    }

}
