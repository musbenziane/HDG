    TYPE SETup
        SEQUENCE
        INTEGER                                         :: N, nt, esrc, gsrc, isnap, ne
        REAL(kind=8)                                    :: h, f0, dt, Jc, Jci
        REAL(kind=8), DIMENSION(:),   ALLOCATABLE       :: v1D, rho1D, src, xi, wi
        REAL(kind=8), DIMENSION(:,:), ALLOCATABLE       :: lprime
        integer, dimension(:,:), allocatable            :: Cij


    END TYPE