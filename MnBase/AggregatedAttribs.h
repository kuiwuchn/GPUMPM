#ifndef __AGGREGATED_ATTRIBS_H_
#define __AGGREGATED_ATTRIBS_H_

#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <array>

namespace mn {

	template<int numAttribs>
	class AttribPort {	///< indices to certain attributes
	public:
		__host__ __device__ AttribPort() {}
		__host__ __device__ ~AttribPort() {}
		// selective link copy
		__inline__ __host__ void link(void **_hdAttrs, int *stencil) {
			for (int i = 0; i < numAttribs; i++)
				checkCudaErrors(cudaMemcpy(_ddAttrs + i, _hdAttrs + stencil[i], sizeof(void*), cudaMemcpyHostToHost));
		}
		// continuous link copy
		__inline__ __host__ void link(void **_hdAttrs, int stpos) {
			checkCudaErrors(cudaMemcpy(_ddAttrs, _hdAttrs + stpos, sizeof(void*)*numAttribs, cudaMemcpyHostToHost));
		}
	protected:
		void*	_ddAttrs[numAttribs];
	};

	template<int numAttribs, int numPorts = 1>
	class AttribConnector {
	public:
		AttribConnector() = default;
		~AttribConnector() {
			for (auto& attrib : _attribs)
				cudaFree(attrib);
			for (auto& port: _ports)
				delete port;
		}

	protected:
		void* port(unsigned int i) const {
			//static_assert(i < numPorts && _ports[i] != nullptr);
			return _ports[i];
		}
		/// manage data
		std::array<void*, numAttribs>	_attribs;		///< cuda allocated arrays
		/// distribute ports
		std::array<void*, numPorts>		_ports{ nullptr };
	};

}

#endif
