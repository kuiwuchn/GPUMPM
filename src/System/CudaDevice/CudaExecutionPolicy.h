#ifndef __CUDA_EXECUTION_POLICY_H_
#define __CUDA_EXECUTION_POLICY_H_
#include <string>

namespace mn {

	struct LaunchInput {	///< could contain more information on operation (error checking/ time recording/ etc...)
		LaunchInput() = delete;
		LaunchInput(std::string kernel, int taskNum, size_t sharedMemBytes = 0) :
			kernelName(kernel), numThreads(taskNum), sharedMemBytes(sharedMemBytes) {}
		const std::string&	name() { return kernelName; }
		const int&			threads() { return numThreads; }
		const size_t&		memBytes() { return sharedMemBytes; }
	private:
		const std::string	kernelName;
		const int			numThreads;
		const size_t		sharedMemBytes;
	};

	/// kernel launching configuration
	struct ExecutionPolicy {
		ExecutionPolicy() {}
		ExecutionPolicy(int gs, int bs, size_t memsize, bool s) :
			gridSize(gs), blockSize(bs), sharedMemBytes(memsize), sync(s) {}
		int			getGridSize() const { return gridSize; }
		int			getBlockSize() const { return blockSize; }
		size_t		getSharedMemBytes() const { return sharedMemBytes; }
		bool		needSync() const { return sync; }
	private:
		int			gridSize{ 0 };
		int			blockSize{ 0 };
		size_t		sharedMemBytes{ 0 };
		bool		sync{ false };
	};

}

#endif