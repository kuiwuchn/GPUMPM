#ifndef __CUDA_DEVICE_UTILS_CUH_
#define __CUDA_DEVICE_UTILS_CUH_

#include <stdint.h>
#include <cuda_runtime.h>

namespace mn {

    __device__ bool atomicMinf(float* address, float val);
	__device__ bool atomicMaxf(float* address, float val);
	__device__ bool atomicMinD(double* address, double val);
	__device__ bool atomicMaxD(double* address, double val);

    __device__ uint64_t Packed_Add(const uint64_t* masks,const uint64_t i,const uint64_t j);
    __device__ uint64_t Bit_Spread_Mine(const uint64_t mask, const int data);
    __device__ int Bit_Pack_Mine(const uint64_t mask, const uint64_t data);

}

#endif
