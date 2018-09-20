#include <helper_cuda.h>
#include <cuda_runtime.h>
#include <Setting.h>

#include <System/CudaDevice/CudaDevice.h>
#include <System/Log/Logger.hpp>

#include <Simulation/MPM/Simulator.h>
#include <Simulation/MPM/SimulatorKernels.cuh>
#include <Simulation/TimeIntegrator/GridUpdateKernels.cuh>
#include <Simulation/TimeIntegrator/P2GKernels.cuh>
#include <Simulation/TimeIntegrator/G2PKernels.cuh>
#include <Simulation/TimeIntegrator/MPMComputationKernels.cuh>
#include <MnBase/Math/Matrix/MatrixKernels.cuh>
#include <Benchmarks.cuh>

using namespace mn;

void mn::checkCurrentCudaError(std::string msg) {
    cudaDeviceSynchronize(); // for print
    auto error = cudaGetLastError();
    if (error != cudaSuccess) printf("%s %s\n", msg.c_str(), cudaGetErrorString(error));
}

int main() {
    mn::Logger::startup();
    mn::CudaDevice::startup();

    mn::CudaDevice::registerKernel("AosToSoa",							(mn::KernelFunc)aosToSoa,							cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("SoaToAos",							(mn::KernelFunc)soaToAos,							cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("CalcOffset",						(mn::KernelFunc)calcOffset,							cudaFuncCachePreferL1, false);   
    mn::CudaDevice::registerKernel("RegisterPage",                      (mn::KernelFunc)registerPage,                       cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("FindPage",                          (mn::KernelFunc)findPage,                           cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("ReorderKey",                        (mn::KernelFunc)reorderKey,                         cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("UpdateIndices",						(mn::KernelFunc)updateIndices,						cudaFuncCachePreferL1, false);
	mn::CudaDevice::registerKernel("Gather3D",					        (mn::KernelFunc)gather3D,						    cudaFuncCachePreferL1, false);
	mn::CudaDevice::registerKernel("Gather3DShared",					(mn::KernelFunc)gather3DShared,						cudaFuncCachePreferL1, false);
	
    mn::CudaDevice::registerKernel("MarkPageBoundary",					(mn::KernelFunc)markPageBoundary,					cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("MarkCellBoundary",					(mn::KernelFunc)markCellBoundary,					cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("MarkBlockOffset",					(mn::KernelFunc)markBlockOffset,					cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("MarkPageSize",						(mn::KernelFunc)markPageSize,						cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("MarkVirtualPageOffset",				(mn::KernelFunc)markVirtualPageOffset,				cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("CalcMaxVel",						(mn::KernelFunc)calcMaxVel,							cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("CalcIndex",							(mn::KernelFunc)calcIndex,							cudaFuncCachePreferL1, false);

    mn::CudaDevice::registerKernel("InitMatrix",						(mn::KernelFunc)initMatrix,							cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("G2P_Implicit_Compute_dP_dF",		(mn::KernelFunc)G2P_Implicit_Compute_dP_dF,			cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("ComputeContributionFixedCorotated", (mn::KernelFunc)computeContributionFixedCorotated,	cudaFuncCachePreferL1, false);

    mn::CudaDevice::registerKernel("ApplyGravity",						(mn::KernelFunc)applyGravity,						cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("PostP2G",							(mn::KernelFunc)postP2G,							cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("PreG2P",							(mn::KernelFunc)preG2P,								cudaFuncCachePreferL1, false);
    mn::CudaDevice::registerKernel("UpdateVelocity",					(mn::KernelFunc)updateVelocity,						cudaFuncCachePreferL1, false);

    {
        printf("Begin Program!\n");
        mn::MPMSimulator simulator;
        mn::Benchmarks benchmarks(simulator);
        benchmarks.run();
    }

    mn::CudaDevice::shutdown();
    mn::Logger::shutdown();

}

