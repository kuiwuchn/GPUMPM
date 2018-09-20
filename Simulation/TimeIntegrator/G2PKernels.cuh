#ifndef __G2P_KERNELS_CUH_
#define __G2P_KERNELS_CUH_
#include <cuda_runtime.h>
#include <stdint.h>
#include <Setting.h>

namespace mn {

    __global__ void G2P_FLIP( 
		const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, const int** smallest_nodes, 
        int* d_block_offsets, int* d_cellids, int* d_indices, int* d_indexTrans, T** d_sorted_positions, T** d_sorted_velocities,
        T** d_channels, T* d_sorted_F, T* d_tmp, T dt, int** d_adjPage);

    __global__ void G2P_APIC( 
		const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, const int** smallest_nodes, 
        int* d_block_offsets, int* d_cellids, int* d_indices, int* d_indexTrans, T** d_sorted_positions, T** d_sorted_velocities, 
		T** d_channels, T* d_sorted_F, T* d_sorted_B, T* d_tmp, T dt, int** d_adjPage);

    __global__ void G2P_MLS( 
		const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, const int** smallest_nodes, 
        int* d_block_offsets, int* d_cellids, int* d_indices, int* d_indexTrans, T** d_sorted_positions, T** d_sorted_velocities, 
		T** d_channels, T* d_sorted_F, T* d_sorted_B, T* d_tmp, T dt, int** d_adjPage);

    __global__ void G2P_Implicit( 
		const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, const int** smallest_nodes, 
        int* d_block_offsets, int* d_cellids, T** d_sorted_positions, T** d_sorted_velocities,
        T** d_channels, T* d_sorted_F, T dt, int** d_adjPage, T** d_implicit_x, T mu, T lambda, T volume, T* d_tmp_matrix);

    __global__ void G2P_Implicit_MLS( 
		const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, const int** smallest_nodes, 
        int* d_block_offsets, int* d_cellids, T** d_sorted_positions, T** d_sorted_velocities,
        T** d_channels, T* d_sorted_F, T dt, int** d_adjPage, T** d_implicit_x, T mu, T lambda, T volume, T* d_tmp_matrix);

    __global__ void G2P_Implicit_Compute_dP_dF(const int numParticle, int* d_indexTrans, T mu, T lambda, T volume,const T* d_sorted_F, T* d_tmp_matrix);
}

#endif
