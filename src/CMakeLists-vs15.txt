#	Required version of CMake
cmake_minimum_required(VERSION 3.0)

#	Project title
project(MPM)

add_subdirectory(VS15CppBuild)
include_directories(./ ${PROJECT_SOURCE_DIR}/)

#	C++

#	CUDA

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("-- _host_name variable is set to: " ${_host_name})

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("-- _host_name variable is set to: " ${_host_name})
if(${_host_name} STREQUAL "mine-win10")

find_package(CUDA REQUIRED)
include_directories(
	"C:/SDKs/cuda/v9.2/Toolkit/include"
	"C:/SDKs/cuda/v9.2/Samples/common/inc")
link_directories("C:/SDKs/cuda/v9.2/Toolkit/lib/x64")

include_directories(C:/SDKs/Partio)   
link_directories(C:/SDKs/Partio)
set(CUDA_SEPARABLE_COMPILATION ON)

else()

find_package(CUDA REQUIRED)
include_directories(D:/Libs/CUDAToolkit/v8.0/include)
link_directories(D:/Libs/CUDAToolkit/v8.0/lib/x64)

include_directories(D:/Libs/partio)   
link_directories(D:/Libs/partio) 
set(CUDA_SEPARABLE_COMPILATION ON)

endif()

set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
message("arch=61")
    
cuda_add_executable(p2g_cuda
    main.cpp

	Benchmarks.cuh
    Benchmarks.cu
  
	helper_cuda.h
	helper_string.h
	setting.h
	SPGrid_Allocator.h
	SPGrid_Array.h
	SPGrid_Geometry.h
	SPGrid_Geometry.cpp
	SPGrid_Mask.h
	SPGrid_Set.h
	SPGrid_Utilities.h
	SPGrid_Utilities.cpp	
	std_array.h
	
	# MnBase
	MnBase/Singleton.h
	MnBase/AggregatedAttribs.h
	MnBase/Math/Matrix/MatrixKernels.cu
	MnBase/Math/Matrix/MatrixKernels.cuh
	MnBase/Math/Matrix/svd3.h
	MnBase/Math/Probability/Probability.h
	MnBase/Math/Probability/Probability.cpp
	MnBase/Meta/AllocMeta.h
	MnBase/Meta/Meta.h
	
	# Simulation
	Simulation/ConstitutiveModel/ConstitutiveModelKernels.cuh
	
	Simulation/DomainTransform/DomainTransformer.cpp
	Simulation/DomainTransform/DomainTransformer.cu
	Simulation/DomainTransform/DomainTransformer.cuh
    Simulation/DomainTransform/DomainTransformKernels.cu
	Simulation/DomainTransform/DomainTransformKernels.cuh
	
	Simulation/Geometry/Geometry.h
	
	Simulation/Grid/GridDomain.cpp
	Simulation/Grid/GridDomain.cu
	Simulation/Grid/GridDomain.cuh
    Simulation/Grid/GridKernels.cu
	Simulation/Grid/GridKernels.cuh
	    
	Simulation/Material/Material.cpp
	Simulation/Material/Material.h
	
	Simulation/Model/Model.cuh
		
	Simulation/MPM/ModelInitSchemes.hpp
	Simulation/MPM/Simulator.cpp
	Simulation/MPM/Simulator.cu
	Simulation/MPM/Simulator.h
    Simulation/MPM/SimulatorKernels.cu
	Simulation/MPM/SimulatorKernels.cuh
	
	Simulation/Particle/ParticleDomain.cpp
    Simulation/Particle/ParticleDomain.cu
	Simulation/Particle/ParticleDomain.cuh
    Simulation/Particle/ParticleKernels.cu
	Simulation/Particle/ParticleKernels.cuh
    
	Simulation/TimeIntegrator/ExplicitTimeIntegrator.cu
	Simulation/TimeIntegrator/G2PKernels.cu
	Simulation/TimeIntegrator/G2PKernels.cuh
	Simulation/TimeIntegrator/GridUpdateKernels.cu
	Simulation/TimeIntegrator/GridUpdateKernels.cuh
	Simulation/TimeIntegrator/ImplicitTimeIntegrator.cu
	Simulation/TimeIntegrator/MPMComputationKernels.cu
	Simulation/TimeIntegrator/MPMComputationKernels.cuh
    Simulation/TimeIntegrator/MPMTimeIntegrator.cu
    Simulation/TimeIntegrator/P2GKernels.cu
    Simulation/TimeIntegrator/P2GKernels.cuh
    Simulation/TimeIntegrator/TimeIntegrator.cuh 
	
	# System
	System/CudaDevice/CudaDevice.h
	System/CudaDevice/CudaDevice.cpp
	System/CudaDevice/CudaDeviceUtils.cuh
	System/CudaDevice/CudaDeviceUtils.cu
	System/CudaDevice/CudaExecutionPolicy.h  
    System/CudaDevice/CudaHostUtils.hpp
	System/CudaDevice/CudaKernelLauncher.cu  
	
	System/Log/Logger.hpp
	
	# Utility
	Utility/PoissonDisk/cyCore.h
	Utility/PoissonDisk/cyHeap.h
	Utility/PoissonDisk/cyPoint.h
	Utility/PoissonDisk/cyPointCloud.h
	Utility/PoissonDisk/cySampleElim.h
	Utility/PoissonDisk/PoissonDisk.h
	Utility/PoissonDisk/SampleGenerator.h

	Utility/Profiler/CudaTimer.hpp
	Utility/Profiler/Timer.hpp
    )
target_link_libraries(p2g_cuda
    mpm_cpp
    cudart
    cudadevrt
	partio
    )

