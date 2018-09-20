#include "G2PKernels.cuh"
#include <MnBase/Math/Matrix/MatrixKernels.cuh>
#include <System/CudaDevice/CudaDeviceUtils.cuh>
#include <Simulation/ConstitutiveModel/ConstitutiveModelKernels.cuh>
#include <cstdio>

namespace mn {

    __global__ void G2P_MLS(
        const int numParticle, 
        const int* d_targetPages, 
        const int* d_virtualPageOffsets, 
        const int** smallest_nodes,
        int* d_block_offsets,
        int* d_cellids,
        int* d_indices,
        int* d_indexTrans,
        T** d_sorted_positions,
        T** d_sorted_velocities,
        T** d_channels,
        T* d_sorted_F,
        T* d_B,
        T* d_tmp,
        T dt,
        int** d_adjPage
    ) {
        __shared__ T buffer[3][8][8][8];
    
	int pageid = d_targetPages[blockIdx.x] - 1; // from virtual to physical page
    	int cellid = d_block_offsets[pageid]; // 
    	int relParid = 512 * (blockIdx.x - d_virtualPageOffsets[pageid]) + threadIdx.x;
    	int parid = cellid + relParid;
    
        int block = threadIdx.x & 0x3f;
        int ci=block >> 4;
        int cj=(block & 0xc) >> 2;
        int ck=block & 3;
    
        block=threadIdx.x >> 6;
        int bi=block >> 2;
        int bj=(block & 2) >> 1;
        int bk=block & 1;

        int page_idx = block ? d_adjPage[block - 1][pageid] : pageid;

        // vel
        for(int v=0;v<3;++v)
                buffer[v][bi*4+ci][bj*4+cj][bk*4+ck] = 
                    *((T*)((uint64_t)d_channels[1+v]+(int)page_idx*4096)+(ci*16+cj*4+ck));

        __syncthreads();

        int smallest_node[3];
        if(relParid < d_block_offsets[pageid + 1] - d_block_offsets[pageid]) {
            cellid = d_cellids[parid] - 1;
            T wOneD[3][3]; 
    
            smallest_node[0] = smallest_nodes[0][cellid];
            smallest_node[1] = smallest_nodes[1][cellid];
            smallest_node[2] = smallest_nodes[2][cellid];
    
			T xp[3];
			xp[0] = d_sorted_positions[0][parid] - smallest_node[0] * dx;
			xp[1] = d_sorted_positions[1][parid] - smallest_node[1] * dx;
			xp[2] = d_sorted_positions[2][parid] - smallest_node[2] * dx;

			for (int v = 0; v < 3; ++v) {
				T d0 = xp[v]*one_over_dx;
				T z = ((T)1.5 - d0);
				wOneD[v][0] = (T)0.5 * z * z;
				d0 = d0 - 1.0f;
				wOneD[v][1] = (T)0.75 - d0 * d0;	
				z = (T)1.5 - (1.0f - d0);
				wOneD[v][2] = (T)0.5 * z * z;
			}

			int c = 0;
			float tmp[27];
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < 3; ++k) {
						tmp[c++] = wOneD[0][i] * wOneD[1][j] * wOneD[2][k];
					}
				}
			}

            for(int v=0;v<3;++v) 
                smallest_node[v]=smallest_node[v] & 0x3;

            T val[9]; for(int i=0;i<3;++i) val[i]=0.f;

			c = 0;
            for(int i=0;i<3;++i) {
            for(int j=0;j<3;++j) {
            for(int k=0;k<3;++k) {
                // v_pic

                val[0] += tmp[c] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
				val[1] += tmp[c] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
				val[2] += tmp[c++] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
            }
            }
            }

            d_tmp[parid                 ] = val[0];
            d_tmp[parid + numParticle   ] = val[1];
            d_tmp[parid + numParticle *2] = val[2];

            d_sorted_positions[0][parid]+=val[0]*dt;
            d_sorted_positions[1][parid]+=val[1]*dt;
            d_sorted_positions[2][parid]+=val[2]*dt;

			for (int i = 0; i<9; ++i) val[i] = 0.f;

			c = 0;
			for (int i = 0; i<3; ++i) {
				for (int j = 0; j<3; ++j) {
					for (int k = 0; k<3; ++k) {
						// B
						val[0] += tmp[c] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[1] += tmp[c] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[2] += tmp[c] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[3] += tmp[c] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[4] += tmp[c] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[5] += tmp[c] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[6] += tmp[c] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);
						val[7] += tmp[c] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);
						val[8] += tmp[c++] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);
					}
				}
			}

            for(int i=0;i<9;++i) d_tmp[parid+(i+3)*numParticle]=val[i];

            for(int i=0;i<9;++i) val[i]=val[i]*dt*D_inverse;
            val[0]+=1.f; val[4]+=1.f; val[8]+=1.f;

            T F[9];
            int parid_trans = d_indexTrans[parid];
            for(int i=0;i<9;++i) F[i]=d_sorted_F[parid_trans+i*numParticle];

            T result[9];
            matrixMatrixMultiplication(&(val[0]),F,result);

            for(int i=0;i<9;++i) d_tmp[parid+(i+12)*numParticle]=result[i];
        }
    
    }

    __global__ void G2P_APIC(
        const int numParticle, 
        const int* d_targetPages, 
        const int* d_virtualPageOffsets, 
        const int** smallest_nodes,
        int* d_block_offsets,
        int* d_cellids,
        int* d_indices,
        int* d_indexTrans,
        T** d_sorted_positions,
        T** d_sorted_velocities,
        T** d_channels,
        T* d_sorted_F,
        T* d_B,
        T* d_tmp,
        T dt,
        int** d_adjPage
    ) {
		__shared__ T buffer[3][8][8][8];

		int pageid = d_targetPages[blockIdx.x] - 1; // from virtual to physical page
		int cellid = d_block_offsets[pageid]; // 
		int relParid = 512 * (blockIdx.x - d_virtualPageOffsets[pageid]) + threadIdx.x;
		int parid = cellid + relParid;

		int block = threadIdx.x & 0x3f;
		int ci = block >> 4;
		int cj = (block & 0xc) >> 2;
		int ck = block & 3;

		block = threadIdx.x >> 6;
		int bi = block >> 2;
		int bj = (block & 2) >> 1;
		int bk = block & 1;

		int page_idx = block ? d_adjPage[block - 1][pageid] : pageid;

		// vel
		for (int v = 0; v<3; ++v)
			buffer[v][bi * 4 + ci][bj * 4 + cj][bk * 4 + ck] =
			*((T*)((uint64_t)d_channels[1 + v] + (int)page_idx * 4096) + (ci * 16 + cj * 4 + ck));

		__syncthreads();

		int smallest_node[3];
		if (relParid < d_block_offsets[pageid + 1] - d_block_offsets[pageid]) {
			cellid = d_cellids[parid] - 1;
			// quadratic B spline weights

			T wOneD[3][3], wgOneD[3][3];

			smallest_node[0] = smallest_nodes[0][cellid];
			smallest_node[1] = smallest_nodes[1][cellid];
			smallest_node[2] = smallest_nodes[2][cellid];

			T xp[3];
			xp[0] = d_sorted_positions[0][parid] - smallest_node[0] * dx;
			xp[1] = d_sorted_positions[1][parid] - smallest_node[1] * dx;
			xp[2] = d_sorted_positions[2][parid] - smallest_node[2] * dx;

			for (int v = 0; v < 3; ++v) {
				T d0 = xp[v]*one_over_dx;
				T z = ((T)1.5 - d0);
				wOneD[v][0] = (T)0.5 * z * z;
				wgOneD[v][0] = -z;
				d0 = d0 - 1;
				wOneD[v][1] = (T)0.75 - d0 * d0;
				wgOneD[v][1] = -d0 * 2;
				z = (T)1.5 - (1 - d0);
				wOneD[v][2] = (T)0.5 * z * z;
				wgOneD[v][2] = z;
			}

			wgOneD[0][0] *= one_over_dx;
			wgOneD[0][1] *= one_over_dx;
			wgOneD[0][2] *= one_over_dx;
			wgOneD[1][0] *= one_over_dx;
			wgOneD[1][1] *= one_over_dx;
			wgOneD[1][2] *= one_over_dx;
			wgOneD[2][0] *= one_over_dx;
			wgOneD[2][1] *= one_over_dx;
			wgOneD[2][2] *= one_over_dx;

			for (int v = 0; v<3; ++v) smallest_node[v] = smallest_node[v] & 0x3;

			T val[9]; for (int i = 0; i<9; ++i) val[i] = 0.f;

			for (int i = 0; i<3; ++i) {
				for (int j = 0; j<3; ++j) {
					for (int k = 0; k<3; ++k) {
						// v_pic
						val[0] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k]*buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[1] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k]*buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[2] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k]*buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
					}
				}
			}

			__syncthreads();
            d_tmp[parid                 ] = val[0];
            d_tmp[parid + numParticle   ] = val[1];
            d_tmp[parid + numParticle *2] = val[2];
			d_sorted_positions[0][parid] += val[0] * dt;
			d_sorted_positions[1][parid] += val[1] * dt;
			d_sorted_positions[2][parid] += val[2] * dt;

			for (int i = 0; i<9; ++i) val[i] = 0.f;
			for (int i = 0; i<3; ++i) {
				for (int j = 0; j<3; ++j) {
					for (int k = 0; k<3; ++k) {
						// F
						val[0] += buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[1] += buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[2] += buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[3] += buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[4] += buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[5] += buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[6] += buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
						val[7] += buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
						val[8] += buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
					}
				}
			}

			for (int i = 0; i<9; ++i) val[i] = val[i] * dt;
			val[0] += 1.f;
			val[4] += 1.f;
			val[8] += 1.f;

			T F[9];
			__syncthreads();
            int parid_trans=d_indexTrans[parid];
			for (int i = 0; i<9; ++i) F[i] = d_sorted_F[parid_trans + i*numParticle];

			T result[9];
			matrixMatrixMultiplication(&(val[0]), F, result);

			for (int i = 0; i<9; ++i) d_tmp[parid + (i+3)*numParticle] = result[i];

			for (int i = 0; i<9; ++i) val[i] = 0.f;
			for (int i = 0; i<3; ++i) {
				for (int j = 0; j<3; ++j) {
					for (int k = 0; k<3; ++k) {
						// B
						val[0] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[1] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[2] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (i*dx - xp[0]);
						val[3] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[4] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[5] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (j*dx - xp[1]);
						val[6] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);
						val[7] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);
						val[8] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * (k*dx - xp[2]);

					}
				}
			}

			for (int i = 0; i<9; ++i) d_tmp[parid + (i+12)*numParticle] = val[i];
		}
    }

    __global__ void G2P_FLIP(
        const int numParticle, 
        const int* d_targetPages, 
        const int* d_virtualPageOffsets, 
        const int** smallest_nodes,
        int* d_block_offsets,
        int* d_cellids,
        int* d_indices,
        int* d_indexTrans,
        T** d_sorted_positions,
        T** d_sorted_velocities,
        T** d_channels,
        T* d_sorted_F,
        T* d_tmp,
        T dt,
        int** d_adjPage
    ) {
		const static T flip = 0.95f;

		__shared__ T buffer[6][8][8][8];

		int pageid = d_targetPages[blockIdx.x] - 1; // from virtual to physical page
		int cellid = d_block_offsets[pageid]; // 
		int relParid = 512 * (blockIdx.x - d_virtualPageOffsets[pageid]) + threadIdx.x;
		int parid = cellid + relParid;

		int block = threadIdx.x & 0x3f;
		int ci = block >> 4;
		int cj = (block & 0xc) >> 2;
		int ck = block & 3;

		block = threadIdx.x >> 6;
		int bi = block >> 2;
		int bj = (block & 2) >> 1;
		int bk = block & 1;

		int page_idx = block ? d_adjPage[block - 1][pageid] : pageid;

		// vel0
		for (int v = 0; v < 3; ++v)
			buffer[v][bi * 4 + ci][bj * 4 + cj][bk * 4 + ck] =
			*((T*)((uint64_t)d_channels[7 + v] + (int)page_idx * 4096) + (ci * 16 + cj * 4 + ck));
		// vel
		for (int v = 0; v < 3; ++v)
			buffer[v + 3][bi * 4 + ci][bj * 4 + cj][bk * 4 + ck] =
			*((T*)((uint64_t)d_channels[1 + v] + (int)page_idx * 4096) + (ci * 16 + cj * 4 + ck));

		__syncthreads();

		int smallest_node[3];
		if (relParid < d_block_offsets[pageid + 1] - d_block_offsets[pageid]) {
			cellid = d_cellids[parid] - 1;
			// quadratic B spline weights

			T wOneD[3][3], wgOneD[3][3];
			smallest_node[0] = smallest_nodes[0][cellid];
			smallest_node[1] = smallest_nodes[1][cellid];
			smallest_node[2] = smallest_nodes[2][cellid];

			for (int v = 0; v < 3; ++v) {
				T d0 = (d_sorted_positions[v][parid] - (T)smallest_node[v] * dx)*one_over_dx;
				T z = ((T)1.5 - d0);
				wOneD[v][0] = (T)0.5 * z * z;
				wgOneD[v][0] = -z;
				d0 = d0 - 1;
				wOneD[v][1] = (T)0.75 - d0 * d0;
				wgOneD[v][1] = -d0 * 2;
				z = (T)1.5 - (1 - d0);
				wOneD[v][2] = (T)0.5 * z * z;
				wgOneD[v][2] = z;
			}

			wgOneD[0][0] *= one_over_dx;
			wgOneD[0][1] *= one_over_dx;
			wgOneD[0][2] *= one_over_dx;
			wgOneD[1][0] *= one_over_dx;
			wgOneD[1][1] *= one_over_dx;
			wgOneD[1][2] *= one_over_dx;
			wgOneD[2][0] *= one_over_dx;
			wgOneD[2][1] *= one_over_dx;
			wgOneD[2][2] *= one_over_dx;

			for (int v = 0; v < 3; ++v) smallest_node[v] = smallest_node[v] & 0x3;

			T val[9];
			val[0] = 0.0f; val[1] = 0.0f; val[2] = 0.0f;
			val[3] = 0.0f; val[4] = 0.0f; val[5] = 0.0f;

			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < 3; ++k) {
						// v_diff
						val[0] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[0][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[1] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[1][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[2] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[2][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						// v_pic
						val[3] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[3][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[4] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[4][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
						val[5] += wOneD[0][i] * wOneD[1][j] * wOneD[2][k] * buffer[5][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k];
					}
				}
			}
			float local_dt = dt;

			__syncthreads();
            int parid_mapped=d_indexTrans[parid];
            d_tmp[parid                 ] = (val[3] * (1.0f - flip) + (val[0] + d_sorted_velocities[0][parid_mapped]) * flip);
            d_tmp[parid + numParticle   ] = (val[4] * (1.0f - flip) + (val[1] + d_sorted_velocities[1][parid_mapped]) * flip);
            d_tmp[parid + numParticle *2] = (val[5] * (1.0f - flip) + (val[2] + d_sorted_velocities[2][parid_mapped]) * flip);

			d_sorted_positions[0][parid] += val[3] * local_dt;
			d_sorted_positions[1][parid] += val[4] * local_dt;
			d_sorted_positions[2][parid] += val[5] * local_dt;

			for (int i = 0; i < 9; ++i) val[i] = 0.f;

			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					for (int k = 0; k < 3; ++k) {
						// matrix : again column major
						val[0] += buffer[3][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[1] += buffer[4][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[2] += buffer[5][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wgOneD[0][i] * wOneD[1][j] * wOneD[2][k];
						val[3] += buffer[3][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[4] += buffer[4][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[5] += buffer[5][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wgOneD[1][j] * wOneD[2][k];
						val[6] += buffer[3][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
						val[7] += buffer[4][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
						val[8] += buffer[5][smallest_node[0] + i][smallest_node[1] + j][smallest_node[2] + k] * wOneD[0][i] * wOneD[1][j] * wgOneD[2][k];
					}
				}
			}

			for (int i = 0; i < 9; ++i) val[i] = val[i] * local_dt;
			val[0] += 1.f; val[4] += 1.f; val[8] += 1.f;

			__syncthreads();

			T F[9], result[9];
            for (int i = 0; i < 9; ++i) F[i] = d_sorted_F[parid_mapped + i*numParticle];

            matrixMatrixMultiplication(&(val[0]), F, result); 

            for (int i = 0; i < 9; ++i) d_tmp[parid + numParticle * (3 + i)] = result[i];
		}
    }

    __global__ void G2P_Implicit_MLS( const int numParticle, const int* d_targetPages, const int* d_virtualPageOffsets, 
                                      const int** smallest_nodes, int* d_block_offsets, int* d_cellids, T** d_sorted_positions,
                                      T** d_sorted_velocities, T** d_channels, T* d_sorted_F, T dt, int** d_adjPage, T** d_implicit_x,
                                      T mu, T lambda, T volume, T* d_tmp_matrix
    ) {
        __shared__ T buffer[3][8][8][8];
    
        int pageid = d_targetPages[blockIdx.x] - 1; // from virtual to physical page
    	int cellid = d_block_offsets[pageid]; // 
    	int relParid = 512 * (blockIdx.x - d_virtualPageOffsets[pageid]) + threadIdx.x;
    	int parid = cellid + relParid;
    
        int block = threadIdx.x & 0x3f;
        int ci=block >> 4;
        int cj=(block & 0xc) >> 2;
        int ck=block & 3;
    
        block=threadIdx.x >> 6;
        int bi=block >> 2;
        int bj=(block & 2) >> 1;
        int bk=block & 1;

        int page_idx = block ? d_adjPage[block - 1][pageid] : pageid;

        for(int v=0;v<3;++v)
            buffer[v][bi*4+ci][bj*4+cj][bk*4+ck] = 
                *((T*)((uint64_t)d_implicit_x[v]+(int)page_idx*4096)+(ci*16+cj*4+ck));
        __syncthreads();

        int smallest_node[3];
        if(relParid < d_block_offsets[pageid + 1] - d_block_offsets[pageid]) {
            cellid = d_cellids[parid] - 1;

            T wOneD[3][3];

            smallest_node[0] = smallest_nodes[0][cellid];
            smallest_node[1] = smallest_nodes[1][cellid];
            smallest_node[2] = smallest_nodes[2][cellid];

            T xp[3];
            xp[0] = d_sorted_positions[0][parid] - smallest_node[0] * dx;
            xp[1] = d_sorted_positions[1][parid] - smallest_node[1] * dx;
            xp[2] = d_sorted_positions[2][parid] - smallest_node[2] * dx;
            for (int v = 0; v < 3; ++v) {
                T d0 = xp[v] * one_over_dx;
                T z = ((T)1.5 - d0);
                wOneD[v][0] = (T)0.5 * z * z;
                d0 = d0 - 1;
                wOneD[v][1] = (T)0.75 - d0 * d0;
                z = (T)1.5 - (1 - d0);
                wOneD[v][2] = (T)0.5 * z * z;
            }

            for(int v=0;v<3;++v) smallest_node[v] = smallest_node[v] & 0x3;
    
            T val[9];
            for(int i=0;i<9;++i)  val[i]=0.f;

            T weight;
            T vel[3];
            T xi_minus_xp[3];
         
            for (int j = 0; j<3; ++j) {
                for (int k = 0; k<3; ++k) {
                    weight = wOneD[0][0] * wOneD[1][j] * wOneD[2][k];
                    /*weight *= D_inverse;*/

                    vel[0] = buffer[0][smallest_node[0]][smallest_node[1] + j][smallest_node[2] + k];
                    vel[1] = buffer[1][smallest_node[0]][smallest_node[1] + j][smallest_node[2] + k];
                    vel[2] = buffer[2][smallest_node[0]][smallest_node[1] + j][smallest_node[2] + k];

                    xi_minus_xp[0] = - xp[0]; xi_minus_xp[1] = j*dx - xp[1]; xi_minus_xp[2] = k*dx - xp[2];

                    val[0] += weight*vel[0] * xi_minus_xp[0];
                    val[1] += weight*vel[1] * xi_minus_xp[0];
                    val[2] += weight*vel[2] * xi_minus_xp[0];
                    val[3] += weight*vel[0] * xi_minus_xp[1];
                    val[4] += weight*vel[1] * xi_minus_xp[1];
                    val[5] += weight*vel[2] * xi_minus_xp[1];
                    val[6] += weight*vel[0] * xi_minus_xp[2];
                    val[7] += weight*vel[1] * xi_minus_xp[2];
                    val[8] += weight*vel[2] * xi_minus_xp[2];

                }
            }

            for (int j = 0; j<3; ++j) {
                for (int k = 0; k<3; ++k) {
                    weight = wOneD[0][1] * wOneD[1][j] * wOneD[2][k];
                    /*weight *= D_inverse;*/

                    vel[0] = buffer[0][smallest_node[0] + 1][smallest_node[1] + j][smallest_node[2] + k];
                    vel[1] = buffer[1][smallest_node[0] + 1][smallest_node[1] + j][smallest_node[2] + k];
                    vel[2] = buffer[2][smallest_node[0] + 1][smallest_node[1] + j][smallest_node[2] + k];

                    xi_minus_xp[0] = dx - xp[0]; xi_minus_xp[1] = j*dx - xp[1]; xi_minus_xp[2] = k*dx - xp[2];

                    val[0] += weight*vel[0] * xi_minus_xp[0];
                    val[1] += weight*vel[1] * xi_minus_xp[0];
                    val[2] += weight*vel[2] * xi_minus_xp[0];
                    val[3] += weight*vel[0] * xi_minus_xp[1];
                    val[4] += weight*vel[1] * xi_minus_xp[1];
                    val[5] += weight*vel[2] * xi_minus_xp[1];
                    val[6] += weight*vel[0] * xi_minus_xp[2];
                    val[7] += weight*vel[1] * xi_minus_xp[2];
                    val[8] += weight*vel[2] * xi_minus_xp[2];

                }
            }

            for (int j = 0; j<3; ++j) {
                for (int k = 0; k<3; ++k) {
                    weight = wOneD[0][2] * wOneD[1][j] * wOneD[2][k];
                    /*weight *= D_inverse;*/

                    vel[0] = buffer[0][smallest_node[0] + 2][smallest_node[1] + j][smallest_node[2] + k];
                    vel[1] = buffer[1][smallest_node[0] + 2][smallest_node[1] + j][smallest_node[2] + k];
                    vel[2] = buffer[2][smallest_node[0] + 2][smallest_node[1] + j][smallest_node[2] + k];

                    xi_minus_xp[0] = 2*dx - xp[0]; xi_minus_xp[1] = j*dx - xp[1]; xi_minus_xp[2] = k*dx - xp[2];

                    val[0] += weight*vel[0] * xi_minus_xp[0];
                    val[1] += weight*vel[1] * xi_minus_xp[0];
                    val[2] += weight*vel[2] * xi_minus_xp[0];
                    val[3] += weight*vel[0] * xi_minus_xp[1];
                    val[4] += weight*vel[1] * xi_minus_xp[1];
                    val[5] += weight*vel[2] * xi_minus_xp[1];
                    val[6] += weight*vel[0] * xi_minus_xp[2];
                    val[7] += weight*vel[1] * xi_minus_xp[2];
                    val[8] += weight*vel[2] * xi_minus_xp[2];

                }
            }

            for(int i=0;i<9;++i) d_tmp_matrix[parid+i*numParticle]=val[i]*D_inverse;
        }
    }

    __global__ void G2P_Implicit(
            const int numParticle, 
            const int* d_targetPages, 
            const int* d_virtualPageOffsets, 
            const int** smallest_nodes,
            int* d_block_offsets,
            int* d_cellids,
            T** d_sorted_positions,
            T** d_sorted_velocities,
            T** d_channels,
            T* d_sorted_F,
            T dt,
            int** d_adjPage,
            T** d_implicit_x,
            T mu,
            T lambda,
            T volume,
            T* d_tmp_matrix
    ) {
        __shared__ T buffer[3][8][8][8];
    
	    int pageid = d_targetPages[blockIdx.x] - 1; // from virtual to physical page
    	int cellid = d_block_offsets[pageid]; // 
    	int relParid = 512 * (blockIdx.x - d_virtualPageOffsets[pageid]) + threadIdx.x;
    	int parid = cellid + relParid;
    
        int block = threadIdx.x & 0x3f;
        int ci=block >> 4;
        int cj=(block & 0xc) >> 2;
        int ck=block & 3;
    
        block=threadIdx.x >> 6;
        int bi=block >> 2;
        int bj=(block & 2) >> 1;
        int bk=block & 1;

        int page_idx = block ? d_adjPage[block - 1][pageid] : pageid;

        for(int v=0;v<3;++v)
            buffer[v][bi*4+ci][bj*4+cj][bk*4+ck] = *((T*)((uint64_t)d_implicit_x[v]+(int)page_idx*4096)+(ci*16+cj*4+ck));
        __syncthreads();

        int smallest_node[3];
        if(relParid < d_block_offsets[pageid + 1] - d_block_offsets[pageid]) {
            cellid = d_cellids[parid] - 1;

            T wOneD[3][3], wgOneD[3][3];
    
            smallest_node[0] = smallest_nodes[0][cellid];
            smallest_node[1] = smallest_nodes[1][cellid];
            smallest_node[2] = smallest_nodes[2][cellid];

            for (int v = 0; v < 3; ++v) {
                T d0 = (d_sorted_positions[v][parid] - smallest_node[v] * dx) * one_over_dx;
                T z = ((T)1.5 - d0);
                wOneD[v][0] = (T)0.5 * z * z;
                wgOneD[v][0] = -z;
                d0 = d0 - 1;
                wOneD[v][1] = (T)0.75 - d0 * d0;
                wgOneD[v][1] = -d0 * 2;
                z = (T)1.5 - (1 - d0);
                wOneD[v][2] = (T)0.5 * z * z;
                wgOneD[v][2] = z;
            }

            wgOneD[0][0] *= one_over_dx;
            wgOneD[0][1] *= one_over_dx;
            wgOneD[0][2] *= one_over_dx;
            wgOneD[1][0] *= one_over_dx;
            wgOneD[1][1] *= one_over_dx;
            wgOneD[1][2] *= one_over_dx;
            wgOneD[2][0] *= one_over_dx;
            wgOneD[2][1] *= one_over_dx;
            wgOneD[2][2] *= one_over_dx; 

            for(int v=0;v<3;++v) smallest_node[v] = smallest_node[v] & 0x3;
    
            T val[9];
            for(int i=0;i<9;++i) val[i]=0.f;

            for(int i=0;i<3;++i) {
            for(int j=0;j<3;++j) {
            for(int k=0;k<3;++k) {
    
                T wg[3];
                wg[0]=wgOneD[0][i]*wOneD[1][j]*wOneD[2][k];
                wg[1]=wOneD[0][i]*wgOneD[1][j]*wOneD[2][k];
                wg[2]=wOneD[0][i]*wOneD[1][j]*wgOneD[2][k];
    
                T vel[3];
                vel[0]=buffer[0][smallest_node[0]+i][smallest_node[1]+j][smallest_node[2]+k];
                vel[1]=buffer[1][smallest_node[0]+i][smallest_node[1]+j][smallest_node[2]+k];
                vel[2]=buffer[2][smallest_node[0]+i][smallest_node[1]+j][smallest_node[2]+k];

                val[0]+=vel[0]*wg[0];
                val[1]+=vel[1]*wg[0];
                val[2]+=vel[2]*wg[0];
                val[3]+=vel[0]*wg[1];
                val[4]+=vel[1]*wg[1];
                val[5]+=vel[2]*wg[1];
                val[6]+=vel[0]*wg[2];
                val[7]+=vel[1]*wg[2];
                val[8]+=vel[2]*wg[2];
    
            }
            }
            }

            for(int i=0;i<9;++i) d_tmp_matrix[parid+i*numParticle]=val[i];
        }
    }

    __global__ void G2P_Implicit_Compute_dP_dF(const int numParticle, int* d_indexTrans, T mu, T lambda, T volume, const T* d_sorted_F, T* d_tmp_matrix)
    {
        int parid = blockDim.x * blockIdx.x + threadIdx.x;
        if (parid >= numParticle) return;
   
        int parid_trans = d_indexTrans[parid];
        T val[9]; for(int i=0;i<9;++i) val[i]=d_tmp_matrix[parid+i*numParticle];
        T F[9]; for(int i=0;i<9;++i) F[i]=d_sorted_F[parid_trans+i*numParticle];
        T dF[9]; matrixMatrixMultiplication(val, F, dF);
        T dP[9];
        Times_Rotated_dP_dF_FixedCorotated(mu, lambda, F, parid, dF, dP);

		matrixMatrixTransposeMultiplication(dP, F, dF);

        for(int i=0;i<9;++i) d_tmp_matrix[parid+i*numParticle]=dF[i]*volume;
    }

}
