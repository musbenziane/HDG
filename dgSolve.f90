subroutine dgSolve(SET,DG,SEM,it)
    !$USE OMP_LIB
    implicit none
        INCLUDE "set.h"
        INCLUDE "dg.h"
        INCLUDE "sem.h"


    TYPE(semProblem)                                           :: SEM
    TYPE(dgProblem)                                            :: DG
    TYPE(SETup)                                                :: SET
    INTEGER                                                    :: i, j, el, it

    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                    ! GettINg GLL poINts and weights

    IF (SET%esrc .gt. SEM%ne) THEN
        DG%source(SET%esrc,SET%gsrc,2) =  SET%src(it) * SET%wi(SET%gsrc) * SET%Jc
    END IF

    call compute_flux(DG%ne,DG%u,SET%N,DG%Al,DG%Ar,DG%flux,SEM,SET)

    !$OMP PARALLEL DO PRIVATE(el) SHARED(SET,DG) SCHEDULE(static)
    do el=1,DG%ne
        DG%k1(el,:,1)   = MATMUL(DG%Minv,(-DG%mu(el) * MATMUL(DG%Ke, DG%u(el,:,2))) - DG%flux(el,:,1) + DG%source(el,:,1))
        DG%k1(el,:,2)   = MATMUL(DG%Minv,(-1. / SET%rho1D(el)) * MATMUL(DG%Ke,DG%u(el,:,1)) - DG%flux(el,:,2) &
                                                                                            + DG%source(el,:,2))
    end do
    !$OMP END PARALLEL DO 

    !$OMP PARALLEL DO PRIVATE(el) SHARED(SET,DG) SCHEDULE(static)
    do el=1,DG%ne
        DG%unew(el,:,1) = SET%dt * MATMUL(DG%Minv,(-DG%mu(el) * MATMUL(DG%Ke,DG%u(el,:,2))) - &
                                          DG%flux(el,:,1) + DG%source(el,:,1)) + DG%u(el,:,1)
        DG%unew(el,:,2) = SET%dt * MATMUL(DG%Minv,(-1. / SET%rho1D(el)) * MATMUL(DG%Ke,DG%u(el,:,1)) - &
                                          DG%flux(el,:,2) + DG%source(el,:,2)) + DG%u(el,:,2)
    end do
    !$OMP END PARALLEL DO

    call compute_flux(DG%ne,DG%unew,SET%N,DG%Al,DG%Ar,DG%flux,SEM,SET)

    !$OMP PARALLEL DO PRIVATE(el) SHARED(SET,DG) SCHEDULE(static)
    do el=1,DG%ne
        DG%k2(el,:,1)   = MATMUL(DG%Minv,(-DG%mu(el) * MATMUL(DG%Ke, DG%unew(el,:,2))) -  &
                                                              DG%flux(el,:,1) + DG%source(el,:,1))
        DG%k2(el,:,2)   = MATMUL(DG%Minv,(-1. / SET%rho1D(el)) * MATMUL(DG%Ke,DG%unew(el,:,1)) - &
                                                                        DG%flux(el,:,2) + DG%source(el,:,2))
    end do
    !$OMP END PARALLEL DO

    DG%unew = DG%u + .5 * SET%dt * (DG%k1 + DG%k2)
    DG%u    = DG%unew

        
    
    end subroutine