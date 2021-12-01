subroutine semSolve(SET, SEM, it)
    !$ USE OMP_LIB
    implicit none
    INCLUDE "sem.h"
    INCLUDE "set.h"

    TYPE(semProblem), INTENT(INOUT)               :: SEM
    TYPE(SETup), INTENT(INOUT)                    :: SET
    REAL(kind=8)                                  :: temp1(SEM%ngll), temp2(SEM%ngll), temp3(SEM%ngll), tmp, tmpo, tmpn
    INTEGER                                       :: i, j, k, l, el, it

    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                                ! GettINg GLL poINts and weights

    !SEM%tauL(SET%Cij(SET%N+1,200)) = SEM%recsigma(it)
    IF (SET%esrc .le. SEM%ne) THEN
        SEM%F(SET%Cij(SET%gsrc,SET%esrc)) =  SET%src(it) * SET%wi(SET%gsrc) * SET%Jc
    END if

    SEM%unew(:)     = SEM%u(:) + SET%dt * SEM%udot(:) + (SET%dt**2)/2 * SEM%uddot(:)

    temp1(:) = 0
    !$OMP PARALLEL DO PRIVATE(k) SHARED(SEM,SET) SCHEDULE(static) 
    do l=1, SEM%ngll
        temp1(l) =  DOT_PRODUCT(SEM%Kg(l , : ), SEM%unew(:))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL WORKSHARE
    temp2 = SEM%F(:) - temp1(:) + SEM%T(:)
    !$OMP END PARALLEL WORKSHARE

    !$OMP PARALLEL WORKSHARE
    SEM%uddotnew(:)  = (1/SEM%M(:)) * temp2(:)
    !$OMP END PARALLEL WORKSHARE

    !$OMP PARALLEL WORKSHARE
    SEM%udotnew(:)  = SEM%udot(:) + (SET%dt/2) * (SEM%uddot(:) + SEM%uddotnew(:))
    !$OMP END PARALLEL WORKSHARE

    do el=1,SEM%ne
        do i=1,SET%N+1
            tmp = 0.
            do j=1,SET%N+1
                tmpo = tmpo + SEM%mu1Dgll(SET%Cij(i,el)) * SEM%uold(SET%Cij(j,el)) * SET%lprime(j,i) * SET%Jci    
                tmp  = tmp  + SEM%mu1Dgll(SET%Cij(i,el)) * SEM%u(SET%Cij(j,el))    * SET%lprime(j,i) * SET%Jci
                tmpn = tmpn + SEM%mu1Dgll(SET%Cij(i,el)) * SEM%unew(SET%Cij(j,el)) * SET%lprime(j,i) * SET%Jci        
            end do
            SEM%sigmaold(SET%Cij(i,el)) = tmpo    ! Traction history
            SEM%sigma(SET%Cij(i,el))    = tmp
            SEM%sigmanew(SET%Cij(i,el)) = tmpn

        end do 
    end do

    ! Traction at the boundary
    SEM%tauLold = SEM%sigmaold(SET%Cij(SET%N+1,SEM%ne))
    SEM%tauL    = SEM%sigma(SET%Cij(SET%N+1,SEM%ne))
    SEM%tauLnew = SEM%sigmanew(SET%Cij(SET%N+1,SEM%ne))   

    SEM%uold  = SEM%u        
    SEM%u     = SEM%unew
    SEM%udot  = SEM%udotnew
    SEM%uddot = SEM%uddotnew
    
    !recsigma(t)  = sigma(Cij(N+1,200))

    
end subroutine