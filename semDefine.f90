subroutine semDefine(SET, SEM)

    !$ USE OMP_LIB
    implicit none
    INCLUDE "sem.h"
    INCLUDE "set.h"

    TYPE(semProblem), INTENT(INOUT)                            :: SEM
    TYPE(SETup), INTENT(INOUT)                                 :: SET
    REAL(kind=8)                                               :: sum
    INTEGER                                                    :: i, j, k, el

    !##########################################
    !####### Construct the mass matrix ########
    !##########################################
    
    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                                ! GettINg GLL poINts and weights

    do i=1,SEM%ne
        do j=1,SET%N+1
            SEM%Me(j)           =  SEM%rho1Dgll(SET%Cij(j,i)) * SET%wi(j) * SET%Jc   ! Elemental mass matrix construction
            SEM%M(SET%Cij(j,i)) =  SEM%M(SET%Cij(j,i)) + SEM%Me(j)                  ! Filling the global mass matrix
        end do
    end do

    SEM%Minv(:,:) = 0
    ! Invert mass matrix
    do i=1,SEM%ngll
        SEM%Minv(i,i) = 1 / SEM%M(i)
    end do



    !###############################################
    !####### Construct the Stiffness matrix ########
    !###############################################
    SEM%mu1Dgll(:) =  SEM%rho1Dgll(:) * SEM%v1Dgll(:)**2.          ! Shear modulus


    SEM%Kg(:,:) = 0
    do el=1,SEM%ne
        !$omp parallel do private(i,j,k,sum) shared(SEM,SET) schedule(static)

        do i=1,SET%N+1                                    ! Elemental stifness matrix
            do j=1,SET%N+1
                sum = 0.
                do k=1,SET%N+1
                    sum = sum + SEM%mu1Dgll(SET%Cij(k,el)) * SET%wi(k) * SET%lprime(j,k) *SET%lprime(i,k) & 
                                                           * SET%Jc * SET%Jci**2
                end do
                SEM%Ke(i,j) = sum
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i,j) shared(SEM,SET) schedule(static)
        do i=1,SET%N+1                                      ! Global assembly
            do j=1,SET%N+1
                SEM%Kg(SET%Cij(i,el),SET%Cij(j,el)) = SEM%Kg(SET%Cij(i,el),SET%Cij(j,el)) + SEM%Ke(i,j)
            end do
        end do
        !$omp end parallel do

    end do

    SEM%uold(:)     = 0.
    SEM%u(:)        = 0.
    SEM%unew(:)     = 0.
    SEM%uddot(:)    = 0.
    SEM%uddotnew(:) = 0.
    SEM%udot(:)     = 0.
    SEM%udotnew(:)  = 0.
    SEM%F(:)        = 0.
    SEM%tauL        = 0.
    SEM%tauLnew     = 0.
    SEM%tauLold     = 0.
    SEM%T           = 0.
    SEM%sigma       = 0.
    SEM%sigmaold    = 0.
    SEM%sigmanew    = 0.
    
end subroutine