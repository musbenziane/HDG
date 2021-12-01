subroutine dgDefine(SET,DG)
    !$ use omp_lib
    implicit none
        INCLUDE "set.h"
        INCLUDE "dg.h"
    
        
        TYPE(dgProblem)                                            :: DG
        TYPE(SETup)                                                :: SET
        INTEGER                                                    :: i, j
        
        CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                    ! GettINg GLL poINts and weights

        DG%mu = SET%v1D**2*SET%rho1D

        !##########################################
        !####### Construct the mass matrix ########
        !##########################################

        do i=1,SET%N+1
            DG%Minv(i,i) = 1. / (SET%wi(i) * SET%Jc)
        end do
    
        !##########################################
        !##### Construct the stiffness matrix #####
        !##########################################
    
        do i=1,SET%N+1
            do j=1,SET%N+1
                   DG%Ke(i,j) = SET%wi(j) * SET%lprime(i,j)
            end do
        end do
    
    
        !##########################################
        !##### Construct the Flux matrices    #####
        !##########################################
    
        DG%Z = SET%rho1D * SET%v1D
    
        !$OMP PARALLEL DO PRIVATE(i) SHARED(SET,DG) SCHEDULE(static)
        do i=1,DG%ne-2
            DG%Ar(i,1,1) =  .5 * SET%v1D(i)
            DG%Ar(i,1,2) = -.5 * DG%Z(i) * SET%v1D(i)
            DG%Ar(i,2,1) = -.5 * SET%v1D(i) / DG%Z(i)
            DG%Ar(i,2,2) =  .5 * SET%v1D(i)
    
            DG%Al(i,1,1) = -.5 * SET%v1D(i)
            DG%Al(i,1,2) = -.5 * DG%Z(i) * SET%v1D(i)
            DG%Al(i,2,1) = -.5 * SET%v1D(i) / DG%Z(i)
            DG%Al(i,2,2) = -.5 * SET%v1D(i);
        end do
        !$OMP END PARALLEL DO
    
    
        
    
    end subroutine