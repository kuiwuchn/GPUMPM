#ifndef __CUDA_HOST_UTILS_HPP_
#define __CUDA_HOST_UTILS_HPP_

#include <string>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

namespace mn {

#define checkThrustErrors(func) \
	try {func;}					\
	catch (thrust::system_error &e) { std::cout << std::string(__FILE__) << ":" << __LINE__ << " " << e.what() << std::endl; }

	template<class T>
	__inline__ __host__ T* getRawPtr(thrust::device_vector<T> &V) {
		return thrust::raw_pointer_cast(V.data());
	}
	template<class T>
	__inline__ __host__ thrust::device_ptr<T> getDevicePtr(thrust::device_vector<T> &V) {
		return thrust::device_ptr<T>(thrust::raw_pointer_cast(V.data()));
	}
	template<class T>
	__inline__ __host__ thrust::device_ptr<T> getDevicePtr(T* V) {
		return thrust::device_ptr<T>(V);
	}

	inline void  reportMemory(std::string msg) {
		size_t free_byte;
		size_t total_byte;
		cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
	
		if ( cudaSuccess != cuda_status ) {
			printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
			exit(1);
		}
	
		double free_db = (double)free_byte; 
		double total_db = (double)total_byte;
		double used_db = total_db - free_db;
		printf("GPU memory usage (%s): used = %f, free = %f MB, total = %f MB\n",
			msg.data(), used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
	}

	void checkCurrentCudaError(std::string msg);
}

#endif