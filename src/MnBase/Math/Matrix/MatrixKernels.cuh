#ifndef __MATRIX_KERNELS_CUH__
#define __MATRIX_KERNELS_CUH__

#include <Setting.h>
#include <cuda_runtime.h>

namespace mn {

    __forceinline__
    __device__ void matrixInverse(const T* x, T* inv) {
        T cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
        T determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;
        T s=1/determinant;
        inv[0]=s*cofactor11; inv[1]=s*cofactor12; inv[2]=s*cofactor13;
        inv[3]=s*x[6]*x[5]-s*x[3]*x[8]; inv[4]=s*x[0]*x[8]-s*x[6]*x[2]; inv[5]=s*x[3]*x[2]-s*x[0]*x[5];
        inv[6]=s*x[3]*x[7]-s*x[6]*x[4]; inv[7]=s*x[6]*x[1]-s*x[0]*x[7]; inv[8]=s*x[0]*x[4]-s*x[3]*x[1];
    }

    __forceinline__
    __device__ T matrixDeterminant(const T* x) {
        return x[0]*(x[4]*x[8]-x[7]*x[5])+x[3]*(x[7]*x[2]-x[1]*x[8])+x[6]*(x[1]*x[5]-x[4]*x[2]);
    }

    __forceinline__
    __device__ void matrixTranspose(const T* x, T* transpose) {
        transpose[0]=x[0]; transpose[1]=x[3]; transpose[2]=x[6];
        transpose[3]=x[1]; transpose[4]=x[4]; transpose[5]=x[7];
        transpose[6]=x[2]; transpose[7]=x[5]; transpose[8]=x[8];
    }


    __forceinline__
    __device__ T matrixTrace(const T* x) {
        return x[0]+x[4]+x[8];
    }

    __forceinline__
    __device__ void matrixMatrixMultiplication(const T* a, const T* b, T* c)
    {
        c[0]=a[0]*b[0]+a[3]*b[1]+a[6]*b[2]; 
        c[1]=a[1]*b[0]+a[4]*b[1]+a[7]*b[2];
        c[2]=a[2]*b[0]+a[5]*b[1]+a[8]*b[2];
        c[3]=a[0]*b[3]+a[3]*b[4]+a[6]*b[5];
        c[4]=a[1]*b[3]+a[4]*b[4]+a[7]*b[5];
        c[5]=a[2]*b[3]+a[5]*b[4]+a[8]*b[5];
        c[6]=a[0]*b[6]+a[3]*b[7]+a[6]*b[8];
        c[7]=a[1]*b[6]+a[4]*b[7]+a[7]*b[8];
        c[8]=a[2]*b[6]+a[5]*b[7]+a[8]*b[8];
    }
    
    __forceinline__
    __device__ void matrixDiagonalMatrixMultiplication(const T* a, const T* b, T* c)
    {
        c[0]=a[0]*b[0];
        c[1]=a[1]*b[0];
        c[2]=a[2]*b[0];
        c[3]=a[3]*b[1];
        c[4]=a[4]*b[1];
        c[5]=a[5]*b[1];
        c[6]=a[6]*b[2];
        c[7]=a[7]*b[2];
        c[8]=a[8]*b[2];
    }

    __forceinline__
    __device__ void matrixTransposeMatrixMultiplication(const T* a, const T* b, T* c) {
        c[0]=a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; 
        c[1]=a[3]*b[0]+a[4]*b[1]+a[5]*b[2];
        c[2]=a[6]*b[0]+a[7]*b[1]+a[8]*b[2];
        c[3]=a[0]*b[3]+a[1]*b[4]+a[2]*b[5];
        c[4]=a[3]*b[3]+a[4]*b[4]+a[5]*b[5];
        c[5]=a[6]*b[3]+a[7]*b[4]+a[8]*b[5];
        c[6]=a[0]*b[6]+a[1]*b[7]+a[2]*b[8];
        c[7]=a[3]*b[6]+a[4]*b[7]+a[5]*b[8];
        c[8]=a[6]*b[6]+a[7]*b[7]+a[8]*b[8];
    }

    __forceinline__
    __device__ __host__ void matrixVectorMultiplication(const T* x,const T* v,T* result)
    {
        result[0]=x[0]*v[0]+x[3]*v[1]+x[6]*v[2];
        result[1]=x[1]*v[0]+x[4]*v[1]+x[7]*v[2];
        result[2]=x[2]*v[0]+x[5]*v[1]+x[8]*v[2];
    }

    __forceinline__
    __device__ void matrixMatrixTransposeMultiplication(const T* a, const T* b, T* c) {
        c[0]=a[0]*b[0]+a[3]*b[3]+a[6]*b[6]; 
        c[1]=a[1]*b[0]+a[4]*b[3]+a[7]*b[6];
        c[2]=a[2]*b[0]+a[5]*b[3]+a[8]*b[6];
        c[3]=a[0]*b[1]+a[3]*b[4]+a[6]*b[7];
        c[4]=a[1]*b[1]+a[4]*b[4]+a[7]*b[7];
        c[5]=a[2]*b[1]+a[5]*b[4]+a[8]*b[7];
        c[6]=a[0]*b[2]+a[3]*b[5]+a[6]*b[8];
        c[7]=a[1]*b[2]+a[4]*b[5]+a[7]*b[8];
        c[8]=a[2]*b[2]+a[5]*b[5]+a[8]*b[8];
    }

    __forceinline__
    __device__ T vectorMaxComponent(const T* x)
    {
        T tmp=x[0];
        if(tmp<x[1]) tmp=x[1];
        if(tmp<x[2]) tmp=x[2];
        return tmp;
    }

    __forceinline__
    __device__ T vectorMagnitude(const T* x)
    {
        return sqrtf(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    }

    __forceinline__
    __device__ void vectorComponentMax(const T* x, const T* y, T* result)
    {
        for(int v=0;v<3;++v)
            result[v]=x[v]>y[v]?x[v]:y[v];
    }

    __forceinline__
    __device__ T signedDistanceOrientedBox(const T* point,const T* box_center,const T* edges,const T* rotation)
    {
        T tmp[3];
        for(int v=0;v<3;++v)
            tmp[v]=point[v]-box_center[v];
        T diff[3];
        matrixVectorMultiplication(rotation,tmp,diff);
        T phi[3];
        for(int v=0;v<3;++v)
            phi[v]=(diff[v]>0?diff[v]:-diff[v])-edges[v]*.5f;

        if(phi[0]<=0 && phi[1]<=0 && phi[2]<=0)
           return vectorMaxComponent(phi); 
        else
        {
            T zeros[3]={0,0,0};
            vectorComponentMax(phi,zeros,diff);
            return vectorMagnitude(diff);
        }
    }
}

#endif
