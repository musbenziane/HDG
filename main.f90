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
        INTEGER                                         :: N, nt, esrc, gsrc, isnap
        REAL(kind=8)                                    :: h, f0, dt, Jc, Jci
        REAL(kind=8), DIMENSION(:),   ALLOCATABLE       :: v1D, rho1D, src, xi, wi
        REAL(kind=8), DIMENSION(:,:), ALLOCATABLE       :: lprime
        integer, dimension(:,:), allocatable            :: Cij

    END TYPE

    TYPE semProblem
        SEQUENCE
        REAL(kind=8), dimension(:),   ALLOCATABLE       :: mu1Dgll, v1Dgll, rho1Dgll, M, xgll, u, uold, unew, F, &
                                                           udot, udotnew, uddot, uddotnew, sigma, recsigma, tauL, Me 
        REAL(kind=8), dimension(:,:), ALLOCATABLE       :: Minv, Kg, Ke, Uout, Udotout, sigmaout
        INTEGER                                         :: ne, ngll
    END TYPE   

    TYPE dgProblem
        SEQUENCE
        REAL(kind=8), DIMENSION(:), ALLOCATABLE                :: Z
        REAL(KIND=8), dimension(:,:), ALLOCATABLE              :: v1D, rho1D, mu1Dgll, v1Dgll, rho1Dgll, xgll, Minv, &
                                                                  Ke, sigma, v
        real(kind=8), dimension(:,:,:), allocatable            :: Al, Ar, u, unew, k1, k2, flux, source
        INTEGER                                                :: ne, ngll
    END TYPE

    TYPE(SETup)                                                :: SET
    TYPE(semProblem)                                           :: SEM
    TYPE(dgProblem)                                            :: DG
    REAL (KIND=8)                                              :: sum, CFL, mINdist, lambdamIN, taper
    REAL (KIND=8)                                              :: time, t_cpu_0, t_cpu_1, t_cpu, tmp, tmpBC, attConst, sd
                                                            
    REAL (KIND=8), dimension(:), allocatable                   ::  g
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

    outname = "OUTPUT/snapshots.bin"
    filename          = "parameters.in"


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
    READ(2,*) SEM%ne
    READ(2,*) DG%ne
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
    PRINT*,"Number of elements        -> ",SEM%ne
    PRINT*,"Element size              -> ",SET%h
    PRINT*,"Wavelet's peak frequency  -> ",SET%f0
    PRINT*,"Time step                 -> ",SET%dt
    PRINT*,"Number of time steps      -> ",SET%nt
    PRINT*,"Source location [nel/ngll]-> ",SET%esrc, SET%gsrc
    PRINT*,"Snapshot INterval         -> ",SET%isnap

    
    SEM%ngll  = SET%N * SEM%ne + 1                                  ! Total GLL poINts
    DG%ngll   = (SET%N + 1) * DG%ne 
    SET%Jc    = SET%h / 2                                           ! Jacobian for structured 1D mesh
    SET%Jci   = 1 / SET%Jc                                          ! Jacobian INverse


    ALLOCATE(SET%Cij(SET%N+1,SEM%ne))                               ! Connectivity matrix
    ALLOCATE(SET%xi(SET%N+1))                                       ! GLL poINts
    ALLOCATE(SET%wi(SET%N+1))                                       ! GLL Quadrature weights
    ALLOCATE(SET%v1D(SEM%ne))                                       ! 1D velocity model IN elements
    ALLOCATE(SET%rho1D(SEM%ne))                                     ! Density velocity model IN elements
    ALLOCATE(SET%lprime(SET%N+1,SET%N+1))                           ! Dervatives of Lagrange polynomials

    !##########################################
    !############ SEM SOLVER INIT #############
    !##########################################


    ALLOCATE(SEM%rho1Dgll(SEM%ngll))                                        ! 1D density model mapped
    ALLOCATE(SEM%v1Dgll(SEM%ngll))                                          ! 1D velocity mapped
    ALLOCATE(SEM%M(SEM%ngll))                                               ! Global mass matrix IN vector form
    ALLOCATE(SEM%MINv(SEM%ngll,SEM%ngll))                                   ! INverse of the mass matrix
    ALLOCATE(SEM%Me(SET%N+1))                                               ! Elemental mass matrix
    ALLOCATE(SEM%Kg(SEM%ngll,SEM%ngll))                                     ! Global stIFness matrix
    ALLOCATE(SEM%Ke(SET%N+1,SET%N+1))                                       ! Elemental stIFness matrix
    ALLOCATE(SEM%u(SEM%ngll))                                               ! Displacement vector at time t
    ALLOCATE(SEM%unew(SEM%ngll),&                                           ! displacement vecotr at time t - dt
             SEM%uddotnew(SEM%ngll), & 
             SEM%udotnew(SEM%ngll))
    ALLOCATE(SEM%uddot(SEM%ngll),    &
             SEM%udot(SEM%ngll), SEM%sigma(SEM%ngll))
    ALLOCATE(SEM%uold(SEM%ngll))                                            ! displacement vector at time t + dt
    ALLOCATE(SET%src(SET%nt))                                               ! Source time FUNCTION
    ALLOCATE(SEM%F(SEM%ngll))                                               ! External force
    ALLOCATE(SEM%Uout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll))               ! Snapshots
    ALLOCATE(SEM%Udotout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll), & 
             SEM%sigmaout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll))
    ALLOCATE(SEM%mu1Dgll(SEM%ngll))                                         ! Shear modulus mapped
    ALLOCATE(SEM%xgll(SEM%ngll))                                            ! Array for global mappINg
    ALLOCATE(g(SEM%ngll))
    ALLOCATE(SEM%tauL(SEM%ngll), SEM%recsigma(SET%nt))
   

    !##########################################
    !############ DG SOLVER INIT  #############
    !##########################################
    ALLOCATE(DG%Z(DG%ne))                                                                 ! Impepdences   
    ALLOCATE(DG%xgll(DG%ne,SET%N+1))                                                      ! Array for global mapping
    ALLOCATE(DG%Minv(SET%N+1,SET%N+1))                                                    ! Elemental mass matrix
    ALLOCATE(DG%Ke(SET%N+1,SET%N+1))                                                      ! Elemental stifness matrix
    ALLOCATE(DG%Ar(DG%ne,2,2),DG%Al(DG%ne,2,2))                                           ! Wave PDE Coefficient Matrix
    ALLOCATE(DG%u(DG%ne,SET%N+1,2),DG%unew(DG%ne,SET%N+1,2))                              ! Solution fields
    ALLOCATE(DG%k1(DG%ne,SET%N+1,2),DG%k2(DG%ne,SET%N+1,2),DG%source(DG%ne,SET%N+1,2))    ! RK
    ALLOCATE(DG%flux(DG%ne,SET%N+1,2))                                                    ! Flux matrix
    ALLOCATE(DG%sigma(NINT(SET%nt/REAL(SET%isnap)),DG%ngll),  &                           ! Stretched solution fields
             DG%v(NINT(SET%nt/REAL(SET%isnap)),DG%ngll))

    
    SET%Cij  = connectivity_matrix(SET%N,SEM%ne)

    CALL lagrangeprime(SET%N,SET%lprime)                                    ! Lagrange polynomials derivatives
    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                                ! GettINg GLL poINts and weights
    CALL READmodelfiles1D(SET%v1D, SET%rho1D, SEM%ne,modnameprefix)         ! READINg model files
    CALL shapefunc(SET%N,SET%h,SEM%ne,SET%Cij,SEM%xgll)                     ! Global domaIN mappINg
    CALL shapefuncDG(SET%N,SET%h,DG%ne, DG%xgll)                            ! Global domain mapping
    CALL mapmodel(SET%N,SEM%ne,SET%rho1D,SET%v1D,SEM%rho1Dgll,SEM%v1Dgll)   ! MappINg models
    CALL ricker(SET%nt,SET%f0,SET%dt,SET%src)                               ! Source time FUNCTION

    CALL semSolve(SET,SEM)


    write(*,*) "##########################################"
    write(*,*) "######### Write solution binary ##########"
    write(*,*) "######### Solution in OUTPUT/   ##########"
    write(*,*) "##########################################"

    outname = "OUTPUT/SEM_snapshots_U.bin"

    inquire(iolength=reclsnaps) SEM%Uout
    open(15,file=outname,access="direct",recl=reclsnaps)
    write(15,rec=1) SEM%Uout
    close(15)


    outname = "OUTPUT/SEM_snapshots_V.bin"

    open(17,file=outname,access="direct",recl=reclsnaps)
    write(17,rec=1) SEM%Udotout
    close(17)


    outname = "OUTPUT/SEM_snapshots_Sigma.bin"

    open(18,file=outname,access="direct",recl=reclsnaps)
    write(18,rec=1) SEM%sigmaout
    close(18)

   
    ! Temps elapsed final.
    call system_clock(count=t1, count_rate=ir)
    time = real(t1 - t0,kind=8) / real(ir,kind=8)

    call cpu_time(t_cpu_1)
    t_cpu = t_cpu_1 - t_cpu_0

    write(*,*) "##########################################"
    write(*,*) "######### TIME:                 ##########"
    print '(//3X,"Elapsed Time        : ",1PE10.3," [s]",/ &
            &,3X,"CPU Time            : ",1PE10.3," [s]",//)', &
            & time,t_cpu
    write(*,*) "##########################################"


END PROGRAM