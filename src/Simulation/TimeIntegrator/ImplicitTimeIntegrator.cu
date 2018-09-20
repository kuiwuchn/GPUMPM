#include "TimeIntegrator.cuh"

#include <MnBase/Math/Matrix/MatrixKernels.cuh>
#include "MPMComputationKernels.cuh"
#include "GridUpdateKernels.cuh"
#include "P2GKernels.cuh"
#include "G2PKernels.cuh"

#include <System/CudaDevice/CudaDeviceUtils.cuh>
#include <System/CudaDevice/CudaKernelLauncher.cu>

namespace mn {

    ImplicitTimeIntegrator::ImplicitTimeIntegrator(int transferScheme, int numParticle, T* dMemTrunk) :
        MPMTimeIntegrator(transferScheme, numParticle, dMemTrunk) {
        checkCudaErrors(cudaMalloc((void**)&d_tmpMatrix, sizeof(T) * numParticle * 9 + 2));
        checkCudaErrors(cudaMalloc((void**)&d_innerProduct, sizeof(double)));
        checkCudaErrors(cudaMalloc((void**)&d_normSquared, sizeof(T)));
    }

    ImplicitTimeIntegrator::~ImplicitTimeIntegrator() {
        checkCudaErrors(cudaFree(d_tmpMatrix));
        checkCudaErrors(cudaFree(d_innerProduct));
        checkCudaErrors(cudaFree(d_normSquared));
    }

    void ImplicitTimeIntegrator::integrate(
        const T dt,
        Model& model, 
        std::unique_ptr<SPGrid>& grid,
        std::unique_ptr<DomainTransformer<TEST_STRUCT<T>>>& trans) 
    {

        auto& geometry = model.refGeometryPtr();

        trans->rebuild();
        Logger::recordSection<TimerType::GPU>("build_particle_grid_mapping");

        configuredLaunch({"CalcIndex", trans->_numCell}, calcIndex,
                         trans->_numCell, one_over_dx, (const int*)(trans->d_cell2particle), 
                         (const T**)geometry->d_orderedPos, geometry->d_smallestNodeIndex);
        Logger::recordSection<TimerType::GPU>("calc_smallest_index");

        computeForceCoefficient(model);
        Logger::recordSection<TimerType::GPU>("contribution_calculation(include_svd)");

        transferP2G(dt,geometry, grid, trans);
        Logger::recordSection<TimerType::GPU>("p2g_calculation");
        call_postP2G(grid, trans);

        applyExternalForce(dt, grid, trans);

        //updateGridVelocity(dt, grid, trans);  ///< explicit
        implicitSolve(dt, model, grid, trans);
        recordLaunch(std::string("ImplicitCopy"), trans-> _numTotalPage, 64, (size_t)0, implicitCopy,
                        0, (const T**)grid->d_implicit_x, 1, grid->d_channels);

        call_preG2P(grid, trans);
        Logger::recordSection<TimerType::GPU>("implicit_grid_update");

        transferG2P(dt, geometry, grid, trans);
        Logger::recordSection<TimerType::GPU>("g2p_calculation");

        geometry->reorder();
        Logger::recordSection<TimerType::GPU>("particle_reorder");
        Logger::blankLine<TimerType::GPU>();
    }

    void ImplicitTimeIntegrator::computeForceCoefficient(Model& model) {
        auto& geometry = model.refGeometryPtr();
        auto& refMaterialPtr = model.refMaterialDynamicsPtr();
        auto material = (ElasticMaterialDynamics*)refMaterialPtr.get();
        configuredLaunch({"ComputeContributionFixedCorotated", geometry->_numParticle}, computeContributionFixedCorotated,
                         geometry->_numParticle, (const T*)geometry->d_F,
                         material->_lambda, material->_mu, material->_volume, d_contribution);
    }

    /// implicit grid update
    void ImplicitTimeIntegrator::multiply(const T dt, 
        Model& model, 
        std::unique_ptr<SPGrid>& grid,
        std::unique_ptr<DomainTransformer<TEST_STRUCT<T>>>& trans,
        T** x, T** result) {

        auto& geometry = model.refGeometryPtr();
        auto material = (ElasticMaterialDynamics*)model.refMaterialDynamicsPtr().get();

#if TRANSFER_SCHEME != 2
        recordLaunch(std::string("G2P_Implicit"), (int)trans->_numVirtualPage, 512, (size_t)0, G2P_Implicit,
                    geometry->_numParticle, 
                    (const int*)trans->d_targetPage,
                    (const int*)trans->d_virtualPageOffset,
                    (const int**)geometry->d_smallestNodeIndex, 
                    trans->d_page2particle, 
                    trans->d_particle2cell, 
                    geometry->d_orderedPos, 
                    geometry->d_orderedVel, 
                    grid->d_channels,
                    geometry->d_F,
                    dt,
                    trans->d_adjPage,
                    // implicit
                    x,
                    material->_mu,
                    material->_lambda,
                    material->_volume,
                    d_tmpMatrix
                    );
#else
        recordLaunch(std::string("G2P_Implicit_MLS"), (int)trans->_numVirtualPage, 512, (size_t)0, G2P_Implicit_MLS,
                    geometry->_numParticle, 
                    (const int*)trans->d_targetPage,
                    (const int*)trans->d_virtualPageOffset,
                    (const int**)geometry->d_smallestNodeIndex, 
                    trans->d_page2particle, 
                    trans->d_particle2cell, 
                    geometry->d_orderedPos,  
                    geometry->d_orderedVel, 
                    grid->d_channels,
                    geometry->d_F,
                    dt,
                    trans->d_adjPage,
                    // implicit
                    x,
                    material->_mu,
                    material->_lambda,
                    material->_volume,
                    d_tmpMatrix
                    );
#endif

        configuredLaunch({"G2P_Implicit_Compute_dP_dF", geometry->_numParticle}, G2P_Implicit_Compute_dP_dF,
                             geometry->_numParticle, geometry->d_indexTrans, material->_mu, material->_lambda, material->_volume, (const T*)geometry->d_F, d_tmpMatrix);

        recordLaunch(std::string("ImplicitClear"), trans-> _numTotalPage, 64, (size_t)0, implicitClear, result);

#if TRANSFER_SCHEME != 2
        recordLaunch(std::string("P2G_Implicit"), (int)trans->_numVirtualPage, 512, (size_t)0, P2G_Implicit,
                     geometry->_numParticle, 
                     (const int*)trans->d_targetPage,
                     (const int*)trans->d_virtualPageOffset,
                     (const int**)geometry->d_smallestNodeIndex, 
                     trans->d_page2particle, 
                     trans->d_particle2cell, 
                     geometry->d_orderedPos, 
                     geometry->d_orderedMass, 
                     geometry->d_orderedVel, 
                     grid->d_channels,
                     trans->d_adjPage,
                     d_tmpMatrix,
                     result
                    );
#else
        recordLaunch(std::string("P2G_Implicit_MLS"), (int)trans->_numVirtualPage, 512, (size_t)0, P2G_Implicit_MLS,
                     geometry->_numParticle, 
                     (const int*)trans->d_targetPage,
                     (const int*)trans->d_virtualPageOffset,
                     (const int**)geometry->d_smallestNodeIndex, 
                     trans->d_page2particle, 
                     trans->d_particle2cell, 
                     geometry->d_orderedPos, 
                     geometry->d_orderedMass, 
                     geometry->d_orderedVel, 
                     grid->d_channels,
                     trans->d_adjPage,
                     d_tmpMatrix,
                     result
                    );
#endif

        bool trapezoidal=false;
        const T scaled_dt_squared=dt*dt/(1+trapezoidal);

        recordLaunch(std::string("ImplicitSystemMatrix"), trans-> _numTotalPage, 64, (size_t)0, implicitSystemMatrix, scaled_dt_squared, grid->d_channels, x, result);
    }

    bool ImplicitTimeIntegrator::implicitSolve(const T dt,
        Model& model, 
        std::unique_ptr<SPGrid>& grid,
        std::unique_ptr<DomainTransformer<TEST_STRUCT<T>>>& trans) {

        auto& geometry = model.refGeometryPtr();
        auto material = (ElasticMaterialDynamics*)model.refMaterialDynamicsPtr().get();

        T nullspace_tolerance=(T)1e-5;
        T residual_magnitude_squared=0;
        T nullspace_measure=0; // extra convergence information
        int restart_iterations=0;

        // some inputs
        T tolerance=5e-5;
        int min_iterations=0;
        int max_iterations=20;

        // compute rhs
        recordLaunch(std::string("UpdateVelocity"), trans-> _numTotalPage, 64, (size_t)0, updateVelocity, dt, grid->d_channels);

		printf("here %d\n", grid->h_masks[0]);

        // set flags for implicit grid flags
        recordLaunch(std::string("SetFlags"), trans-> _numTotalPage, 64, (size_t)0, setFlags, 
                     make_ulonglong3(grid->h_masks[0], grid->h_masks[1], grid->h_masks[2]), trans->d_pageOffset, grid->d_implicit_flags);

        // conjugate residual 
        static const T small_number=1e-7;

        T rho_old=0;T convergence_norm=0; T initial_residual_norm=0.f;
        double h_innerProduct=0.;
        int iterations;
        // precalc weight & weight gradient
        for (iterations=0; ; iterations++){
            bool restart = !iterations || (restart_iterations && iterations%restart_iterations == 0);
            if (restart) {
                //r=b;
                recordLaunch(std::string("ImplicitCopy"), trans->_numTotalPage, 64, (size_t)0, implicitCopy, 1, (const T**)grid->d_channels, 0, grid->d_implicit_r);

                checkCudaErrors(cudaMemset(d_normSquared, 0, sizeof(T)));
                recordLaunch(std::string("ImplicitConvergenceNorm"), trans-> _numTotalPage, 64, (size_t)0, implicitConvergenceNorm, 
                             (const T**)grid->d_implicit_r, d_normSquared);
                checkCudaErrors(cudaMemcpy((void*)&convergence_norm, d_normSquared, sizeof(T), cudaMemcpyDeviceToHost));
                initial_residual_norm = std::sqrt(convergence_norm);
                printf("initial convergence_norm is %f\n",initial_residual_norm);
                
                //system.Multiply(x,p);
                multiply(dt,model, grid, trans, grid->d_implicit_x,grid->d_implicit_p);
                //r-=p;
                recordLaunch(std::string("ImplicitMinus"), trans->_numTotalPage, 64, (size_t)0, implicitMinus, (const T**)grid->d_implicit_p, grid->d_implicit_r);
                //system.Project(r);
                recordLaunch(std::string("ImplicitProject"), trans->_numTotalPage, 64, (size_t)0, implicitProject, (const unsigned*)grid->d_implicit_flags, grid->d_implicit_r);
            }

            //convergence_norm=system.Convergence_Norm(r);
            checkCudaErrors(cudaMemset(d_normSquared, 0, sizeof(T)));
            recordLaunch(std::string("ImplicitConvergenceNorm"), trans->_numTotalPage, 64, (size_t)0, implicitConvergenceNorm, (const T**)grid->d_implicit_r, d_normSquared);
            checkCudaErrors(cudaMemcpy((void*)&convergence_norm, d_normSquared, sizeof(T), cudaMemcpyDeviceToHost));
            convergence_norm = std::sqrt(convergence_norm);
            printf("convergence_norm is %f\n",convergence_norm);

            //residual_magnitude_squared=(T)system.Inner_Product(r,r);
            checkCudaErrors(cudaMemset(d_innerProduct, 0, sizeof(double)));
            recordLaunch(std::string("ImplicitInnerProduct"), (trans->_numTotalPage * 64 / 512) + 1, 512, (size_t)0, implicitInnerProduct,
			//recordLaunch(std::string("ImplicitInnerProduct"), trans->_numTotalPage, 64, (size_t)0, implicitInnerProduct,
                (const T**)grid->d_channels, (const T**)grid->d_implicit_r, (const T**)grid->d_implicit_r, d_innerProduct, trans->_numTotalPage * 64);
            checkCudaErrors(cudaMemcpy((void*)&h_innerProduct, d_innerProduct, sizeof(double), cudaMemcpyDeviceToHost));
            residual_magnitude_squared=(T)h_innerProduct;

            nullspace_measure=(residual_magnitude_squared>small_number*small_number*100)?abs(rho_old/residual_magnitude_squared):0;
            if((convergence_norm<=std::fmin(tolerance*initial_residual_norm,tolerance)  || (iterations && nullspace_measure<=nullspace_tolerance)) &&
               (iterations>=min_iterations || convergence_norm<small_number)){ // TODO: get the stopping criterion right
                printf("implicit solve %d iterations used\n",iterations);
                return true;
            }
            if(iterations==max_iterations) break;

            //system.Multiply(mr,ar);
            multiply(dt, model, grid, trans, grid->d_implicit_r,grid->d_implicit_ar);

            //system.Project(ar);
            recordLaunch(std::string("ImplicitProject"), trans->_numTotalPage, 64, (size_t)0, implicitProject, (const unsigned*)grid->d_implicit_flags, grid->d_implicit_ar);

             //T rho=(T)system.Inner_Product(mr,ar);
            checkCudaErrors(cudaMemset(d_innerProduct, 0, sizeof(double)));
            recordLaunch(std::string("ImplicitInnerProduct"), (trans->_numTotalPage * 64 / 512) + 1, 512, (size_t)0, implicitInnerProduct,
			//recordLaunch(std::string("ImplicitInnerProduct"), trans->_numTotalPage, 64, (size_t)0, implicitInnerProduct,
                (const T**)grid->d_channels, (const T**)grid->d_implicit_r, (const T**)grid->d_implicit_ar, d_innerProduct, trans->_numTotalPage * 64);
            T rho=0.f;
            checkCudaErrors(cudaMemcpy((void*)&h_innerProduct, d_innerProduct, sizeof(double), cudaMemcpyDeviceToHost));
            rho=(T)h_innerProduct;
            if(!rho) break;
            if(restart)
            {
                //p=mr;
                recordLaunch(std::string("ImplicitCopy"), trans->_numTotalPage, 64, (size_t)0, implicitCopy, 0, (const T**)grid->d_implicit_r,0,grid->d_implicit_p);
                //ap=ar;
                recordLaunch(std::string("ImplicitCopy"), trans->_numTotalPage, 64, (size_t)0, implicitCopy, 0, (const T**)grid->d_implicit_ar,0,grid->d_implicit_ap);
            }
            else{
                T beta=rho/rho_old;
                //p.Copy(beta,p,mr);
                recordLaunch(std::string("ImplicitScale"), trans->_numTotalPage, 64, (size_t)0, implicitScale, 
                    beta, (const T**)grid->d_implicit_p, (const T**)grid->d_implicit_r, grid->d_implicit_p);
                //ap.Copy(beta,ap,ar);
                recordLaunch(std::string("ImplicitScale"), trans->_numTotalPage, 64, (size_t)0, implicitScale, 
                    beta, (const T**)grid->d_implicit_ap, (const T**)grid->d_implicit_ar, grid->d_implicit_ap);
            }

            //T denom=(T)system.Inner_Product(map,ap);
            checkCudaErrors(cudaMemset(d_innerProduct, 0, sizeof(double)));
            recordLaunch(std::string("ImplicitInnerProduct"), (trans->_numTotalPage * 64 / 512) + 1, 512, (size_t)0, implicitInnerProduct,
			//recordLaunch(std::string("ImplicitInnerProduct"), trans->_numTotalPage, 64, (size_t)0, implicitInnerProduct,
                (const T**)grid->d_channels, (const T**)grid->d_implicit_ap, (const T**)grid->d_implicit_ap, d_innerProduct, trans->_numTotalPage * 64);
            T denom=0.f;
            checkCudaErrors(cudaMemcpy((void*)&h_innerProduct, d_innerProduct, sizeof(double), cudaMemcpyDeviceToHost));
            denom=h_innerProduct;
            if(!denom) break;
            T alpha=rho/denom;
            //x.Copy(alpha,p,x);
            recordLaunch(std::string("ImplicitScale"), trans->_numTotalPage, 64, (size_t)0, implicitScale, 
                alpha, (const T**)grid->d_implicit_p, (const T**)grid->d_implicit_x, grid->d_implicit_x);
            //r.Copy(-alpha,ap,r);
            recordLaunch(std::string("ImplicitScale"), trans->_numTotalPage, 64, (size_t)0, implicitScale, 
                -alpha, (const T**)grid->d_implicit_ap, (const T**)grid->d_implicit_r, grid->d_implicit_r);
            rho_old=rho;
        }
        return false;
    }
}
