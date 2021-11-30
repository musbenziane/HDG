subroutine dgSolve(SET,DG)
!$ use omp_lib
implicit none
    TYPE SETup
        SEQUENCE
        INTEGER                                         :: N, nt, esrc, gsrc, isnap
        REAL(kind=8)                                    :: h, f0, dt, Jc, Jci
        REAL(kind=8), DIMENSION(:),   ALLOCATABLE       :: v1D, rho1D, src, xi, wi
        REAL(kind=8), DIMENSION(:,:), ALLOCATABLE       :: lprime
        integer, dimension(:,:), allocatable            :: Cij
    END TYPE



    TYPE dgProblem
        SEQUENCE
        REAL(kind=8), DIMENSION(:), ALLOCATABLE                :: Z
        REAL(KIND=8), dimension(:,:), ALLOCATABLE              :: v1D, rho1D, mu1Dgll, v1Dgll, rho1Dgll, xgll, Minv, &
                                                                  Ke, sigma, v
        real(kind=8), dimension(:,:,:), allocatable            :: Al, Ar, u, unew, k1, k2, flux, source
        INTEGER                                                :: ne, ngll
    END TYPE

    
    TYPE(dgProblem)                                            :: DG
    TYPE(SETup)                                                :: SET
    INTEGER                                                    :: i, j, k, l, el, it











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