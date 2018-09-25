#	Required version of CMake
cmake_minimum_required(VERSION 3.0)

#	Project title
project(MPM)

include_directories(../ ${PROJECT_SOURCE_DIR}/)

#	C++


cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("-- _host_name variable is set to: " ${_host_name})
if(${_host_name} STREQUAL "neumann")
    set(CMAKE_C_COMPILER /s/gcc-6.1/amd64_ubu14/bin/gcc)
    set(CMAKE_CXX_COMPILER /s/gcc-6.1/amd64_ubu14/bin/g++)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -c")

    include_directories(/usr/local/cuda-8.0/include)
    include_directories(/usr/local/cuda-8.0/samples/common/inc)
    link_directories(/usr/local/cuda-8.0/lib64)
elseif(${_host_name} STREQUAL "amax")
    set(CMAKE_C_COMPILER /usr/local/bin/gcc)
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -c")

    include_directories(/usr/local/cuda-8.0/include)
    include_directories(/usr/local/cuda-8.0/samples/common/inc)
    link_directories(/usr/local/cuda-8.0/lib64)

    include_directories(~/SDKs/partio/build/Linux-4.4.0-x86_64-optimize/include)
    link_directories(~/SDKs/partio/build/Linux-4.4.0-x86_64-optimize/lib)
elseif(${_host_name} STREQUAL "turing")
    set(CMAKE_C_COMPILER /usr/bin/gcc-7)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-7)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -c")
    #set(CMAKE_C_COMPILER /s/gcc-6.1/amd64_ubu14/bin/gcc)
    #set(CMAKE_CXX_COMPILER /s/gcc-6.1/amd64_ubu14/bin/g++)
    #SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -c")

    include_directories(/usr/local/cuda-9.1/include)
    include_directories(/usr/local/cuda-9.1/samples/common/inc)
    link_directories(/usr/local/cuda-9.1/lib64)

    include_directories(/nobackup/mingg13/partio/build/Linux-4.4.0-x86_64-optimize/include)
    link_directories(/nobackup/mingg13/partio/build/Linux-4.4.0-x86_64-optimize/lib)
    # for libz 
    link_directories(/usr/lib/x86_64-linux-gnu/)
    link_directories(/lib/i386-linux-gnu/)
    link_directories(/lib/x86_64-linux-gnu/)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARIES})
elseif(${_host_name} STREQUAL "metalgear")
    set(CMAKE_C_COMPILER /usr/bin/gcc-7)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-7)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -c")
    #set(CMAKE_C_COMPILER /s/gcc-6.1/amd64_ubu14/bin/gcc)
    #set(CMAKE_CXX_COMPILER /s/gcc-6.1/amd64_ubu14/bin/g++)
    #SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -c")

    include_directories(/usr/local/cuda-9.0/include)
    include_directories(/usr/local/cuda-9.0/samples/common/inc)
    link_directories(/usr/local/cuda-9.0/lib64)

    include_directories(/home/fanfu/partio/partio/build/Linux-4.13.0-x86_64-optimize/include)
    link_directories(   /home/fanfu/partio/partio/build/Linux-4.13.0-x86_64-optimize/lib)
    # for libz 
    link_directories(/usr/lib/x86_64-linux-gnu/)
    link_directories(/lib/i386-linux-gnu/)
    link_directories(/lib/x86_64-linux-gnu/)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARIES})
elseif(${_host_name} STREQUAL "kratos")
    set(CMAKE_C_COMPILER /usr/bin/gcc-7)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-7)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -c")
    #set(CMAKE_C_COMPILER /s/gcc-6.1/amd64_ubu14/bin/gcc)
    #set(CMAKE_CXX_COMPILER /s/gcc-6.1/amd64_ubu14/bin/g++)
    #SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -Wall -c")

    include_directories(/usr/include/thrust/system/cuda/)
    link_directories(/usr/lib/cuda)
    include_directories(/usr/local/cuda-9.2/samples/common/inc)

    include_directories(/home/fanfu/partio/build/Linux-4.15.0-x86_64-optimize/include)
    link_directories(   /home/fanfu/partio/build/Linux-4.15.0-x86_64-optimize/lib)

    # for libz 
    link_directories(/usr/lib/x86_64-linux-gnu/)
    link_directories(/lib/i386-linux-gnu/)
    link_directories(/lib/x86_64-linux-gnu/)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARIES})
elseif(${_host_name} STREQUAL "mine-ubuntu1804")
    set(CMAKE_C_COMPILER /usr/bin/gcc-6)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-6)
    SET(CMAKE_CXX_FLAGS "-std=c++14 -O3")

    include_directories(/usr/local/cuda-9.2/include)
    include_directories(/usr/local/cuda-9.2/samples/common/inc)
    link_directories(/usr/local/cuda-9.2/lib64)

    include_directories(/home/mine/SDKs/partio/include)
    link_directories(/home/mine/SDKs/partio/lib)
    # for libz 
    link_directories(/usr/lib/x86_64-linux-gnu/)
    link_directories(/lib/x86_64-linux-gnu/)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARIES})
elseif(${_host_name} STREQUAL "mine-Surface-Book-2")
    set(CMAKE_C_COMPILER /usr/bin/gcc-7)
    set(CMAKE_CXX_COMPILER /usr/bin/g++-7)
    SET(CMAKE_CXX_FLAGS "-std=c++17 -O3 -c")

    include_directories(/usr/local/cuda-9.1/include)
    include_directories(/usr/local/cuda-9.1/samples/common/inc)
    link_directories(/usr/local/cuda-9.1/lib64)
    
    include_directories(/home/mine/SDK/partio/src/lib)
    link_directories(/home/mine/SDK/partio/src/lib)
    # for libz 
    link_directories(/usr/lib/x86_64-linux-gnu/)
    link_directories(/lib/x86_64-linux-gnu/)
    include_directories(${ZLIB_INCLUDE_DIRS})
    link_directories(${ZLIB_LIBRARIES})
endif()

add_library(mpm_cpp
	../Simulation/Material/Material.cpp
	../Simulation/Particle/ParticleDomain.cpp
	../Simulation/Grid/GridDomain.cpp
	../Simulation/DomainTransform/DomainTransformer.cpp
	../Simulation/MPM/Simulator.cpp
	../Simulation/MPM/SimulatorBuilder.cpp
	../MnBase/Math/Probability/Probability.cpp
)

target_link_libraries(mpm_cpp
                      cudart
                      cudadevrt
                      partio
                      z
                      -L/usr/lib/x86_64-linux-gnu/libz.so
                     )
