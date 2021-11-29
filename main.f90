PROGRAM HDG
    
    !$ USE OMP_LIB
    IMPLICIT NONE

    INTERFACE connectivity_matrix
        FUNCTION connectivity_matrix(N,ne)
            INTEGER, INTENT(IN)                         :: N, ne
            INTEGER connectivity_matrix(N+1,ne)
        END FUNCTION
    END INTERFACE
    
    TYPE SETup
        SEQUENCE
        INTEGER                                         :: N, ne, nt, esrc, gsrc, ngll, isnap
        REAL(KIND=8)                                    :: h, f0, dt, Jc, Jci
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE       :: v1D, rho1D, src, xi, wi
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE       :: lprime
        INTEGER, dimension(:,:), allocatable            :: Cij

    END TYPE

    TYPE semProblem
        SEQUENCE
        REAL(KIND=8), dimension(:),   ALLOCATABLE :: mu1Dgll, v1Dgll, rho1Dgll, M, xgll, u, uold, unew, F, &
                                                     udot, udotnew, uddot, uddotnew, sigma, recsigma, tauL, Me 
        REAL(KIND=8), dimension(:,:), ALLOCATABLE :: MINv, Kg, Ke, Uout, Udotout, sigmaout
    END TYPE

    TYPE dgProblem
        SEQUENCE
        REAL(KIND=8), dimension(:,:), ALLOCATABLE :: v1D, rho1D, mu1Dgll, v1Dgll, rho1Dgll
    END TYPE

    TYPE(SETup)                                                :: SET
    TYPE(semProblem)                                           :: SEM
    REAL (KIND=8)                                              :: sum, CFL, mINdist, lambdamIN, taper
    REAL (KIND=8)                                              :: time, t_cpu_0, t_cpu_1, t_cpu, tmp, tmpBC, attConst, sd
                                                            
    REAL (KIND=8), dimension(:), allocatable                   ::  g
    REAL (KIND=8), dimension(:), allocatable                   :: temp1, temp2, temp3 
    INTEGER                                                    :: i, j, k, l, it, el, reclsnaps
    INTEGER                                                    :: ir, t0, t1, bc, sbc, gWidth, IC
    character(len=40)                                          :: filename, filecheck, outname, modnameprefix
    !$ INTEGER                                                 :: n_workers
    LOGICAL                                                    :: OMPcheck = .false.



    WRITE(*,*) "##########################################"
    WRITE(*,*) "############### OPENMP     ###############"
    !$ OMPcheck = .true.
    IF (OMPcheck) THEN
        !$OMP PARALLEL
        !$ n_workers = OMP_GET_NUM_THREADS()
        !$OMP END PARALLEL
        !$ PRINT '(3X,"Number of workers ->  ",i2)',n_workers
    ELSE
        WRITE(*,*) "Program has been compiled without OPENMP; Time marchINg will run IN serial"
    END IF


    CALL cpu_time(t_cpu_0)
    CALL system_clock(count=t0, count_rate=ir)

    outname = "OUTPUT/snapshots.bIN"
    filename          = "parameters.IN"


    WRITE(*,*) "##########################################"
    WRITE(*,*) "######## READINg SETmeters file #########"
    WRITE(*,*) "##########################################"

    
    PRINT*,"Is the SETmeters INput file (parameters.IN) [Yes/no]"
    READ(*,*) filecheck

    IF (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") THEN
        WRITE(*,*) "READINg simulation parameters..."

    ELSEIF  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") THEN
        WRITE(*,*) "Enter simulation parameters text file name with extension"
        WRITE(*,*) "40 characters max"
        READ(*,*) filename

    ELSE
        WRITE(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        WRITE(*,*) "The program have been termINated, please star over"
        STOP
    END IF

    OPEN (2, file=filename, status = 'old')
    READ(2,*) modnameprefix
    READ(2,*) SET%N
    READ(2,*) SET%ne
    READ(2,*) SET%h
    READ(2,*) SET%f0
    READ(2,*) SET%dt
    READ(2,*) SET%nt
    READ(2,*) SET%esrc
    READ(2,*) SET%gsrc
    READ(2,*) SET%isnap
    READ(2,*) bc
    READ(2,*) sbc
    READ(2,*) gWidth
    READ(2,*) attConst
    READ(2,*) IC
    READ(2,*) sd
    CLOSE(2)

    PRINT*,"Polynomial order          -> ",SET%N
    PRINT*,"Number of elements        -> ",SET%ne
    PRINT*,"Element size              -> ",SET%h
    PRINT*,"Wavelet's peak frequency  -> ",SET%f0
    PRINT*,"Time step                 -> ",SET%dt
    PRINT*,"Number of time steps      -> ",SET%nt
    PRINT*,"Source location [nel/ngll]-> ",SET%esrc, SET%gsrc
    PRINT*,"Snapshot INterval         -> ",SET%isnap

    
    SET%ngll  = SET%N * SET%ne + 1                                  ! Total GLL poINts
    SET%Jc    = SET%h / 2                                           ! Jacobian for structured 1D mesh
    SET%Jci   = 1 / SET%Jc                                          ! Jacobian INverse


    ALLOCATE(SET%Cij(SET%N+1,SET%ne))                               ! Connectivity matrix
    ALLOCATE(SET%xi(SET%N+1))                                       ! GLL poINts
    ALLOCATE(SET%wi(SET%N+1))                                       ! GLL Quadrature weights
    ALLOCATE(SET%v1D(SET%ne))                                       ! 1D velocity model IN elements
    ALLOCATE(SET%rho1D(SET%ne))                                     ! Density velocity model IN elements
    ALLOCATE(SET%lprime(SET%N+1,SET%N+1))                           ! Dervatives of Lagrange polynomials

    !##########################################
    !############ SEM SOLVER INIT #############
    !##########################################


    ALLOCATE(SEM%rho1Dgll(SET%ngll))                                        ! 1D density model mapped
    ALLOCATE(SEM%v1Dgll(SET%ngll))                                          ! 1D velocity mapped
    ALLOCATE(SEM%M(SET%ngll))                                               ! Global mass matrix IN vector form
    ALLOCATE(SEM%MINv(SET%ngll,SET%ngll))                                   ! INverse of the mass matrix
    ALLOCATE(SEM%Me(SET%N+1))                                               ! Elemental mass matrix
    ALLOCATE(SEM%Kg(SET%ngll,SET%ngll))                                     ! Global stIFness matrix
    ALLOCATE(SEM%Ke(SET%N+1,SET%N+1))                                       ! Elemental stIFness matrix
    ALLOCATE(SEM%u(SET%ngll))                                               ! Displacement vector at time t
    ALLOCATE(SEM%unew(SET%ngll),&                                           ! displacement vecotr at time t - dt
             SEM%uddotnew(SET%ngll), & 
             SEM%udotnew(SET%ngll))
    ALLOCATE(SEM%uddot(SET%ngll),    &
             SEM%udot(SET%ngll), SEM%sigma(SET%ngll))
    ALLOCATE(SEM%uold(SET%ngll))                                            ! displacement vector at time t + dt
    ALLOCATE(SET%src(SET%nt))                                               ! Source time FUNCTION
    ALLOCATE(SEM%F(SET%ngll))                                               ! External force
    ALLOCATE(SEM%Uout(NINT(REAL(SET%nt/SET%isnap)),SET%ngll))               ! Snapshots
    ALLOCATE(SEM%Udotout(NINT(REAL(SET%nt/SET%isnap)),SET%ngll), & 
             SEM%sigmaout(NINT(REAL(SET%nt/SET%isnap)),SET%ngll))
    ALLOCATE(SEM%mu1Dgll(SET%ngll))                                         ! Shear modulus mapped
    ALLOCATE(SEM%xgll(SET%ngll))                                            ! Array for global mappINg
    ALLOCATE(g(SET%ngll))
    ALLOCATE(SEM%tauL(SET%ngll), SEM%recsigma(SET%nt))
   
    
    SET%Cij  = connectivity_matrix(SET%N,SET%ne)

    CALL lagrangeprime(SET%N,SET%lprime)                                    ! Lagrange polynomials derivatives
    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                                ! GettINg GLL poINts and weights
    CALL READmodelfiles1D(SET%v1D, SET%rho1D, SET%ne,modnameprefix)         ! READINg model files
    CALL shapefunc(SET%N,SET%h,SET%ne,SET%Cij,SEM%xgll)                     ! Global domaIN mappINg
    CALL mapmodel(SET%N,SET%ne,SET%rho1D,SET%v1D,SEM%rho1Dgll,SEM%v1Dgll)   ! MappINg models
    CALL ricker(SET%nt,SET%f0,SET%dt,SET%src)                               ! Source time FUNCTION

    CALL semSolve(SET,SEM)

END PROGRAM