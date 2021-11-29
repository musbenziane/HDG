subroutine semSolve(SET, SEM)

    !$ USE OMP_LIB
    implicit none
    TYPE SETup
        SEQUENCE
        INTEGER                                         :: N, ne, nt, esrc, gsrc, ngll, isnap
        REAL(kind=8)                                    :: h, f0, dt, Jc, Jci
        REAL(kind=8), DIMENSION(:),   ALLOCATABLE       :: v1D, rho1D, src, xi, wi
        REAL(kind=8), DIMENSION(:,:), ALLOCATABLE       :: lprime
        integer, dimension(:,:), allocatable            :: Cij

    END TYPE

    TYPE semProblem
        SEQUENCE
        REAL(kind=8), dimension(:),   ALLOCATABLE :: mu1Dgll, v1Dgll, rho1Dgll, M, xgll, u, uold, unew, F, &
                                                     udot, udotnew, uddot, uddotnew, sigma, recsigma, tauL, Me 
        REAL(kind=8), dimension(:,:), ALLOCATABLE :: Minv, Kg, Ke, Uout, Udotout, sigmaout
    END TYPE    

    TYPE(semProblem)                                           :: SEM
    TYPE(SETup)                                                :: SET
    REAL(kind=8), DIMENSION(SET%N+1,SET%N+1)                   :: lprime
    INTEGER                                                    :: i, j, k, l, el




    


    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,SET%ne
        do j=1,SET%N+1
            SEM%Me(j)           = SEM%rho1Dgll(SET%Cij(j,i)) * SET%wi(j) * SET%Jc   ! Elemental mass matrix construction
            SEM%M(SET%Cij(j,i)) =  SEM%M(SET%Cij(j,i)) + SEM%Me(j)                  ! Filling the global mass matrix
        end do
    end do

    SEM%Minv(:,:) = 0
    ! Invert mass matrix
    do i=1,SET%ngll
        SEM%Minv(i,i) = 1 / SEM%M(i)
    end do



end subroutine