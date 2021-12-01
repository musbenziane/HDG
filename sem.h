    TYPE semProblem
        SEQUENCE
        REAL(kind=8), dimension(:),   ALLOCATABLE       :: mu1Dgll, v1Dgll, rho1Dgll, M, xgll, u, uold, unew, F, &
                                                           udot, udotnew, uddot, uddotnew, sigma, recsigma, Me, &
                                                           tauL, tauLnew, tauLold, sigmaold, sigmanew, T 
        REAL(kind=8), dimension(:,:), ALLOCATABLE       :: Minv, Kg, Ke, Uout, Udotout, sigmaout
        INTEGER                                         :: ne,ngll
    END TYPE    