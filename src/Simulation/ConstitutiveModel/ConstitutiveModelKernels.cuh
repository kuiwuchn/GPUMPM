#ifndef __CONSTITUTIVE_MODEL_KERNELS_CUH__
#define __CONSTITUTIVE_MODEL_KERNELS_CUH__

#include <Setting.h>
#include <cuda_runtime.h>
#include <MnBase/Math/Matrix/svd3.h>

namespace mn {

    __forceinline__ __device__ T Clamp_Small_Magnitude(const T input)
    {
        T magnitude=input>0?input:-input;
        T sign=input>0?1.f:-1.f;
        T output=magnitude>1e-6?magnitude:1e-6;
        return output*sign;
    }

    __forceinline__  __device__ T log_1px_over_x(const T x, const T eps)
    {
#ifdef DEBUG_INFO
        if(eps<0)
        {
            printf("eps of log_1px_over_x cannot be negative\n");
        }
#endif
        T absX=x>0?x:-x;
        if (absX < eps)
            return 1.f;
        else
            return log1pf(x) / x; // TODO replace with more accrucate version
    }

    __forceinline__ __device__ T diff_log_over_diff(const T x, const T y, const T eps)
    {
#ifdef DEBUG_INFO
        if(eps<0)
        {
            printf("eps of diff_log_over_diff cannot be negative\n");
        }
#endif
        T p = x / y - 1;
        return log_1px_over_x(p, eps) / y;
    }

    __forceinline__ __device__ T diff_interlock_log_over_diff(const T x, const T y, const T logy, const T eps)
    {
#ifdef DEBUG_INFO
        if(eps<0)
        {
            printf("eps of diff_interlock_log_over_diff cannot be negative\n");
        }
#endif
        return logy - y * diff_log_over_diff(x, y, eps);
    }

    __forceinline__ __device__ bool too_close(const T x, const T y, const T eps)
    {
#ifdef DEBUG_INFO
        if(eps<0)
        {
            printf("eps of too_close cannot be negative\n");
        }
#endif
        T p=x/y-1.f;
        T absP=p>0?p:-p;
        if (absP < eps)
            return true;
        else 
            return false;

    }
    __device__ void Times_Rotated_dP_dF_FixedCorotated(const T mu, const T lambda, const T* F, const int parid, const T* dF, T* dP) {

        T U[9]; T S[3]; T V[9];
        svd(F[0],F[3],F[6],F[1],F[4],F[7],F[2],F[5],F[8], U[0],U[3],U[6],U[1],U[4],U[7],U[2],U[5],U[8], S[0],S[1],S[2], V[0],V[3],V[6],V[1],V[4],V[7],V[2],V[5],V[8]);

        // 
        T J=S[0]*S[1]*S[2]; T scaled_mu=2.f*mu; T scaled_lambda=lambda*(J-1.f);
        T P_hat[3]; P_hat[0]=scaled_mu*(S[0]-1.f)+scaled_lambda*(S[1]*S[2]); P_hat[1]=scaled_mu*(S[1]-1.f)+scaled_lambda*(S[0]*S[2]); P_hat[2]=scaled_mu*(S[2]-1.f)+scaled_lambda*(S[0]*S[1]);

        T dP_hat_dSigma_upper[6]; scaled_lambda=lambda*(2.f*J-1.f)*J;
        for(int i=0;i<3;++i) dP_hat_dSigma_upper[i]=scaled_mu+lambda*J*J/(S[i]*S[i]);
        dP_hat_dSigma_upper[3]=scaled_lambda/(S[0]*S[1]); dP_hat_dSigma_upper[4]=scaled_lambda/(S[0]*S[2]); dP_hat_dSigma_upper[5]=scaled_lambda/(S[1]*S[2]);

        scaled_lambda=-lambda*(J-1.f)*J;
        T M[3]; M[0]=0.5f*(2.f*mu+scaled_lambda/(S[0]*S[1])); M[1]=0.5f*(2.f*mu+scaled_lambda/(S[0]*S[2])); M[2]=0.5f*(2.f*mu+scaled_lambda/(S[1]*S[2]));
        //

        T P[3]; 
		P[0]=0.5*(P_hat[0]+P_hat[1])/Clamp_Small_Magnitude(S[0]+S[1]); 
		P[1]=0.5*(P_hat[0]+P_hat[2])/Clamp_Small_Magnitude(S[0]+S[2]); 
		P[2]=0.5*(P_hat[1]+P_hat[2])/Clamp_Small_Magnitude(S[1]+S[2]);

        T dF_hat[9];
        dF_hat[0]=(dF[0] * U[0] + dF[1] * U[1] + dF[2] * U[2]) * V[0] + (dF[3] * U[0] + dF[4] * U[1] + dF[5] * U[2]) * V[1] + (dF[6] * U[0] + dF[7] * U[1] + dF[8] * U[2]) * V[2]; 
        dF_hat[1]=(dF[0] * U[3] + dF[1] * U[4] + dF[2] * U[5]) * V[0] + (dF[3] * U[3] + dF[4] * U[4] + dF[5] * U[5]) * V[1] + (dF[6] * U[3] + dF[7] * U[4] + dF[8] * U[5]) * V[2];
        dF_hat[2]=(dF[0] * U[6] + dF[1] * U[7] + dF[2] * U[8]) * V[0] + (dF[3] * U[6] + dF[4] * U[7] + dF[5] * U[8]) * V[1] + (dF[6] * U[6] + dF[7] * U[7] + dF[8] * U[8]) * V[2];
        dF_hat[3]=(dF[0] * U[0] + dF[1] * U[1] + dF[2] * U[2]) * V[3] + (dF[3] * U[0] + dF[4] * U[1] + dF[5] * U[2]) * V[4] + (dF[6] * U[0] + dF[7] * U[1] + dF[8] * U[2]) * V[5];
        dF_hat[4]=(dF[0] * U[3] + dF[1] * U[4] + dF[2] * U[5]) * V[3] + (dF[3] * U[3] + dF[4] * U[4] + dF[5] * U[5]) * V[4] + (dF[6] * U[3] + dF[7] * U[4] + dF[8] * U[5]) * V[5];
        dF_hat[5]=(dF[0] * U[6] + dF[1] * U[7] + dF[2] * U[8]) * V[3] + (dF[3] * U[6] + dF[4] * U[7] + dF[5] * U[8]) * V[4] + (dF[6] * U[6] + dF[7] * U[7] + dF[8] * U[8]) * V[5];
        dF_hat[6]=(dF[0] * U[0] + dF[1] * U[1] + dF[2] * U[2]) * V[6] + (dF[3] * U[0] + dF[4] * U[1] + dF[5] * U[2]) * V[7] + (dF[6] * U[0] + dF[7] * U[1] + dF[8] * U[2]) * V[8];
        dF_hat[7]=(dF[0] * U[3] + dF[1] * U[4] + dF[2] * U[5]) * V[6] + (dF[3] * U[3] + dF[4] * U[4] + dF[5] * U[5]) * V[7] + (dF[6] * U[3] + dF[7] * U[4] + dF[8] * U[5]) * V[8];
        dF_hat[8]=(dF[0] * U[6] + dF[1] * U[7] + dF[2] * U[8]) * V[6] + (dF[3] * U[6] + dF[4] * U[7] + dF[5] * U[8]) * V[7] + (dF[6] * U[6] + dF[7] * U[7] + dF[8] * U[8]) * V[8];

        T dP_hat[9];
        dP_hat[0] = dP_hat_dSigma_upper[0] * dF_hat[0] + dP_hat_dSigma_upper[3] * dF_hat[4] + dP_hat_dSigma_upper[4] * dF_hat[8];
        dP_hat[4] = dP_hat_dSigma_upper[3] * dF_hat[0] + dP_hat_dSigma_upper[1] * dF_hat[4] + dP_hat_dSigma_upper[5] * dF_hat[8];
        dP_hat[8] = dP_hat_dSigma_upper[4] * dF_hat[0] + dP_hat_dSigma_upper[5] * dF_hat[4] + dP_hat_dSigma_upper[2] * dF_hat[8];
        dP_hat[3] = ((M[0] + P[0]) * dF_hat[3] + (M[0] - P[0]) * dF_hat[1]) ;
        dP_hat[1] = ((M[0] - P[0]) * dF_hat[3] + (M[0] + P[0]) * dF_hat[1]) ;
        dP_hat[6] = ((M[1] + P[1]) * dF_hat[6] + (M[1] - P[1]) * dF_hat[2]) ;
        dP_hat[2] = ((M[1] - P[1]) * dF_hat[6] + (M[1] + P[1]) * dF_hat[2]) ;
        dP_hat[7] = ((M[2] + P[2]) * dF_hat[7] + (M[2] - P[2]) * dF_hat[5]) ;
        dP_hat[5] = ((M[2] - P[2]) * dF_hat[7] + (M[2] + P[2]) * dF_hat[5]) ;

        dP[0]=(dP_hat[0] * U[0] + dP_hat[1] * U[3] + dP_hat[2] * U[6]) * V[0] + (dP_hat[3] * U[0] + dP_hat[4] * U[3] + dP_hat[5] * U[6]) * V[3] + (dP_hat[6] * U[0] + dP_hat[7] * U[3] + dP_hat[8] * U[6]) * V[6]; 
        dP[1]=(dP_hat[0] * U[1] + dP_hat[1] * U[4] + dP_hat[2] * U[7]) * V[0] + (dP_hat[3] * U[1] + dP_hat[4] * U[4] + dP_hat[5] * U[7]) * V[3] + (dP_hat[6] * U[1] + dP_hat[7] * U[4] + dP_hat[8] * U[7]) * V[6];
        dP[2]=(dP_hat[0] * U[2] + dP_hat[1] * U[5] + dP_hat[2] * U[8]) * V[0] + (dP_hat[3] * U[2] + dP_hat[4] * U[5] + dP_hat[5] * U[8]) * V[3] + (dP_hat[6] * U[2] + dP_hat[7] * U[5] + dP_hat[8] * U[8]) * V[6];
        dP[3]=(dP_hat[0] * U[0] + dP_hat[1] * U[3] + dP_hat[2] * U[6]) * V[1] + (dP_hat[3] * U[0] + dP_hat[4] * U[3] + dP_hat[5] * U[6]) * V[4] + (dP_hat[6] * U[0] + dP_hat[7] * U[3] + dP_hat[8] * U[6]) * V[7];
        dP[4]=(dP_hat[0] * U[1] + dP_hat[1] * U[4] + dP_hat[2] * U[7]) * V[1] + (dP_hat[3] * U[1] + dP_hat[4] * U[4] + dP_hat[5] * U[7]) * V[4] + (dP_hat[6] * U[1] + dP_hat[7] * U[4] + dP_hat[8] * U[7]) * V[7];
        dP[5]=(dP_hat[0] * U[2] + dP_hat[1] * U[5] + dP_hat[2] * U[8]) * V[1] + (dP_hat[3] * U[2] + dP_hat[4] * U[5] + dP_hat[5] * U[8]) * V[4] + (dP_hat[6] * U[2] + dP_hat[7] * U[5] + dP_hat[8] * U[8]) * V[7];
        dP[6]=(dP_hat[0] * U[0] + dP_hat[1] * U[3] + dP_hat[2] * U[6]) * V[2] + (dP_hat[3] * U[0] + dP_hat[4] * U[3] + dP_hat[5] * U[6]) * V[5] + (dP_hat[6] * U[0] + dP_hat[7] * U[3] + dP_hat[8] * U[6]) * V[8];
        dP[7]=(dP_hat[0] * U[1] + dP_hat[1] * U[4] + dP_hat[2] * U[7]) * V[2] + (dP_hat[3] * U[1] + dP_hat[4] * U[4] + dP_hat[5] * U[7]) * V[5] + (dP_hat[6] * U[1] + dP_hat[7] * U[4] + dP_hat[8] * U[7]) * V[8];
        dP[8]=(dP_hat[0] * U[2] + dP_hat[1] * U[5] + dP_hat[2] * U[8]) * V[2] + (dP_hat[3] * U[2] + dP_hat[4] * U[5] + dP_hat[5] * U[8]) * V[5] + (dP_hat[6] * U[2] + dP_hat[7] * U[5] + dP_hat[8] * U[8]) * V[8];
    };
}

#endif
