subroutine dgDefine(SET,DG)
    !$ use omp_lib
    implicit none
        INCLUDE "set.h"
        INCLUDE "dg.h"
    
        
        TYPE(dgProblem)                                            :: DG
        TYPE(SETup)                                                :: SET
        INTEGER                                                    :: i, j
        
        CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                    ! GettINg GLL poINts and weights

        DG%mu = SET%v1D(DG%ne:SET%ne)**2*SET%rho1D(DG%ne:SET%ne)

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
        j = DG%ne

        DG%Z = SET%rho1D(DG%ne:SET%ne) * SET%v1D(DG%ne:SET%ne)
    

        do i=1,DG%ne-2
            DG%Ar(i,1,1) =  .5 * SET%v1D(j)
            DG%Ar(i,1,2) = -.5 * DG%Z(i) * SET%v1D(j)
            DG%Ar(i,2,1) = -.5 * SET%v1D(j) / DG%Z(i)
            DG%Ar(i,2,2) =  .5 * SET%v1D(j)
    
            DG%Al(i,1,1) = -.5 * SET%v1D(j)
            DG%Al(i,1,2) = -.5 * DG%Z(i) * SET%v1D(j)
            DG%Al(i,2,1) = -.5 * SET%v1D(j) / DG%Z(i)
            DG%Al(i,2,2) = -.5 * SET%v1D(j);
            j = j + 1
        end do
    
    
        
    
    end subroutine