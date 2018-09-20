#include "CudaDevice.h"
#include <cstdio>

#include <cuda_runtime.h>
#include <cuda_occupancy.h>		///<	for optimal kernel launching
#include <cuda_profiler_api.h>	///<	for evaluating kernel performance
#include <helper_cuda.h>

namespace mn {

	cudaDeviceProp*	CudaDevice::_akDeviceProps;
	int				CudaDevice::_iDevID;
	std::unordered_map<std::string, KernelConfig>	
					CudaDevice::_kFuncTable;

	KernelConfig::KernelConfig(KernelFunc f, cudaFuncCache cacheConfig, bool isWave) :
		func(f), cachePreference(cacheConfig), waveFashion(isWave) {
		cudaFuncGetAttributes(&attribs, f);
		maxOccBlockSize = CudaDevice::evalOptimalBlockSize(attribs, cachePreference);
		if (cacheConfig != cudaFuncCachePreferNone)	///< should be different from device cache preference
			checkCudaErrors(cudaFuncSetCacheConfig(f, cacheConfig));
	}

	CudaDevice::CudaDevice() {
		// acquire devices
		int deviceCount = 0;
		cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
		if (error_id != cudaSuccess) {
			printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
			printf("Result = FAIL\n");
			exit(EXIT_FAILURE);
		}
		if (deviceCount == 0)	printf("There are no available device(s) that support CUDA\n");
		else					printf("Detected %d CUDA Capable device(s)\n", deviceCount);
		_akDeviceProps = new cudaDeviceProp[deviceCount];
		for (int i = 0; i < deviceCount; i++) {
			cudaSetDevice(i);
			cudaGetDeviceProperties(&_akDeviceProps[i], i);
		}
		_iDevID = findCudaDevice(0, nullptr);
		printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n",
			_akDeviceProps[_iDevID].multiProcessorCount, _akDeviceProps[_iDevID].major, 
			_akDeviceProps[_iDevID].minor);

		printf("  Finished \'CudaDevice\' initialization\n");
	}

	CudaDevice::~CudaDevice() {
		delete[] _akDeviceProps;
		printf("  Finished \'CudaDevice\' termination\n");
	}

	int CudaDevice::generalGridSize(int& threadNum, int& blockSize) const { return (threadNum + blockSize - 1) / blockSize; }
	int CudaDevice::waveGridSize(int& threadNum, int& blockSize) const {
		return ((threadNum / blockSize / _akDeviceProps[_iDevID].multiProcessorCount) * _akDeviceProps[_iDevID].multiProcessorCount ?
			(threadNum / blockSize / _akDeviceProps[_iDevID].multiProcessorCount) * _akDeviceProps[_iDevID].multiProcessorCount : 1);
	}
	
	/// static methods
	int CudaDevice::evalOptimalBlockSize(cudaFuncAttributes attribs, cudaFuncCache cachePreference, size_t smemBytes) {
		cudaOccDeviceProp prop = _akDeviceProps[_iDevID];	///< cache preference
		cudaOccFuncAttributes occAttribs = attribs;
		cudaOccDeviceState occCache;
		switch(cachePreference) {
		case cudaFuncCachePreferNone:
			occCache.cacheConfig = CACHE_PREFER_NONE;
			break;
		case cudaFuncCachePreferShared:
			occCache.cacheConfig = CACHE_PREFER_SHARED;
			break;
		case cudaFuncCachePreferL1:
			occCache.cacheConfig = CACHE_PREFER_L1;
			break;
		case cudaFuncCachePreferEqual:
			occCache.cacheConfig = CACHE_PREFER_EQUAL;
			break;
		default:
			;	///< should throw error
		}
		int minGridSize, blockSize;
		cudaOccMaxPotentialOccupancyBlockSize(&minGridSize, &blockSize, &prop, &occAttribs, &occCache, nullptr, smemBytes);
		return blockSize;
	}

	ExecutionPolicy CudaDevice::launchConfig(std::string kernelName, int threadNum, bool sync, size_t smemSize, cudaStream_t sid) const {
		if (_kFuncTable.find(kernelName) == _kFuncTable.end()) {
			int bs = 256;
			printf("Warning: Kernel function %s not registered! Use 256 setting!\n");
			return{ generalGridSize(threadNum, bs), bs, smemSize, sync };
		}
		auto&	config = _kFuncTable[kernelName.data()];
		int		bs = config.maxOccBlockSize;
		if (smemSize > 0)
			bs = evalOptimalBlockSize(config.attribs, config.cachePreference, smemSize);
		//printf("configurating for kernel[%s] blocksize: %d\n", kernelName.c_str(), bs);
		if (config.waveFashion)
			return{ waveGridSize(threadNum, bs), bs, smemSize, sync };
		return{ generalGridSize(threadNum, bs), bs, smemSize, sync };
	}

	void CudaDevice::registerKernel(std::string tag, KernelFunc f, cudaFuncCache cacheConfig, bool waveFashion) {
		_kFuncTable.emplace(tag, KernelConfig(f, cacheConfig, waveFashion));
		printf("Kernel[%s](%s) block size configuration: %d\n", tag.data(), waveFashion?"wave":"general", _kFuncTable[tag.data()].maxOccBlockSize);
	}
	const KernelConfig& CudaDevice::findKernel(std::string tag) {
		return _kFuncTable[tag.data()];
	}

	
}
