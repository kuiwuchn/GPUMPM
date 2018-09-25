#	Required version of CMake
cmake_minimum_required(VERSION 3.0)

#	Project title
project(MPM)

include_directories(../ ${PROJECT_SOURCE_DIR}/)

#	C++
set(CMAKE_C_COMPILER /usr/bin/gcc-6)
set(CMAKE_CXX_COMPILER /usr/bin/g++-6)
SET(CMAKE_CXX_FLAGS "-std=c++14 -O3 -Wall")

#	CUDA

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("-- _host_name variable is set to: " ${_host_name})
if(${_host_name} STREQUAL "neumann")

    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-8.0/include)
    include_directories(/usr/local/cuda-8.0/samples/common/inc)
    link_directories(/usr/local/cuda-8.0/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_52,code=sm_52)
    message("arch=52")
elseif(${_host_name} STREQUAL "amax")

    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-8.0/include)
    include_directories(/usr/local/cuda-8.0/samples/common/inc)
    link_directories(/usr/local/cuda-8.0/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")

    include_directories(~/SDKs/partio/build/Linux-4.4.0-x86_64-optimize/include)
    link_directories(~/SDKs/partio/build/Linux-4.4.0-x86_64-optimize/lib)
elseif(${_host_name} STREQUAL "turing")

    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-9.1/include)
    include_directories(/usr/local/cuda-9.1/samples/common/inc)
    link_directories(/usr/local/cuda-9.1/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")

    include_directories(/nobackup/mingg13/partio/build/Linux-4.4.0-x86_64-optimize/include)
    link_directories(/nobackup/mingg13/partio/build/Linux-4.4.0-x86_64-optimize/lib)
elseif(${_host_name} STREQUAL "kratos")

	find_package(CUDA)
    find_package(CUDA REQUIRED)
    include_directories(/usr/include/thrust/system/cuda/)
    link_directories(/usr/lib/cuda)
    include_directories(/usr/local/cuda-9.2/samples/common/inc)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")

    include_directories(/home/fanfu/partio/build/Linux-4.15.0-x86_64-optimize/include)
    link_directories(   /home/fanfu/partio/build/Linux-4.15.0-x86_64-optimize/lib)
elseif(${_host_name} STREQUAL "metalgear")

	find_package(CUDA)
    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-9.0/include)
    include_directories(/usr/local/cuda-9.0/samples/common/inc)
    link_directories(/usr/local/cuda-9.0/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")

    include_directories(/home/fanfu/partio/partio/build/Linux-4.13.0-x86_64-optimize/include)
    link_directories(   /home/fanfu/partio/partio/build/Linux-4.13.0-x86_64-optimize/lib)
elseif(${_host_name} STREQUAL "mine-ubuntu1804")

    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-9.2/include)
    include_directories(/usr/local/cuda-9.2/samples/common/inc)
    link_directories(/usr/local/cuda-9.2/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")
    
    include_directories(/home/mine/SDKs/partio/include)
    link_directories(/home/mine/SDKs/partio/lib)
elseif(${_host_name} STREQUAL "mine-Surface-Book-2")

    find_package(CUDA REQUIRED)
    include_directories(/usr/local/cuda-9.1/include)
    include_directories(/usr/local/cuda-9.1/samples/common/inc)
    link_directories(/usr/local/cuda-9.1/lib64)
    set(CUDA_SEPARABLE_COMPILATION ON)

    set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_61,code=sm_61)
    message("arch=61")
    
    include_directories(/home/mine/SDK/partio/src/lib)
    link_directories(/home/mine/SDK/partio/src/lib)
endif()

cuda_add_library(mpm_cuda
    ../Simulation/MPM/Simulator.cu
    ../Simulation/MPM/SimulatorKernels.cu
    ../Simulation/Particle/ParticleDomain.cu
    ../Simulation/Particle/ParticleKernels.cu
    ../Simulation/DomainTransform/DomainTransformer.cu
    ../Simulation/DomainTransform/DomainTransformKernels.cu
    ../Simulation/Grid/GridDomain.cu
    ../Simulation/Grid/GridKernels.cu
    ../Simulation/TimeIntegrator/MPMTimeIntegrator.cu
    ../Simulation/TimeIntegrator/ExplicitTimeIntegrator.cu
    ../Simulation/TimeIntegrator/ImplicitTimeIntegrator.cu
    ../Simulation/TimeIntegrator/P2GKernels.cu
    ../Simulation/TimeIntegrator/GridUpdateKernels.cu
    ../Simulation/TimeIntegrator/G2PKernels.cu
    ../Simulation/TimeIntegrator/MPMComputationKernels.cu

    ../Benchmarks.cu

    ../MnBase/Math/Matrix/MatrixKernels.cu
    
    ../System/CudaDevice/CudaDeviceUtils.cu
    ../System/CudaDevice/CudaDevice.cpp
    ../SPGrid_Geometry.cpp
    ../SPGrid_Utilities.cpp
    ../System/CudaDevice/CudaHostUtils.hpp
    )
target_link_libraries(mpm_cuda
    mpm_cpp
    cudart
    cudadevrt
    )


