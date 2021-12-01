    TYPE dgProblem
        SEQUENCE
        REAL(kind=8), DIMENSION(:), ALLOCATABLE                :: Z, mu
        REAL(KIND=8), dimension(:,:), ALLOCATABLE              :: v1D, rho1D, mu1Dgll, v1Dgll, rho1Dgll, xgll, Minv, &
                                                                  Ke, sigma, v
        real(kind=8), dimension(:,:,:), allocatable            :: Al, Ar, u, unew, k1, k2, flux, source
        INTEGER                                                :: ne, ngll
    END TYPE