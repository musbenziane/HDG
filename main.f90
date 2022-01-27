    !###################################################################################################################
    ! HDG (SEM & DG) Hybrod Solver with a regular mesh for the 1D elastic wave equation 
    ! NOTE: This program is yet to be tested against analytical solutions.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Sources:  Igel 2017, Leveque 2002, Hesthaven & Warburton 2008
    ! 
    ! The code used for arbitrary GLL points and weights was created by M.I.T departement of engineering. Link is hereafter
    ! https://geodynamics.org/cig/doxygen/release/specfem3d/    | file name: gll_library.f90
    ! 
    ! This is part of the Numerical Modelling Workshop.
    ! 
    ! Supervisor: Pr. Emmanual Chaljub
    ! Author    : Mus Benziane
    !
    ! Input file: [Example] [Last BC options are yet to be impelented...]
    ! testing              ! Model name prefix (model names: prefix_vp, prefix_rho)
    ! 6                    ! Polynomial order
    ! 200                  ! Number of elements SEM
    ! 200                  ! Number of elements DG
    ! 30.                  ! Element size
    ! 20.                  ! Wavelet's peak frequency
    ! 0.0002               ! Time step
    ! 20000                ! Number of time steps
    ! 150                  ! Source location (element number)
    ! 4                    ! Source location (gll point)
    ! 50                   ! Snapshot interval
    ! 1                    ! [1/2/3/4] 1: Free surface, 2: Rigid wall, 3: Periodic, 4: Sponge 
    ! 0                    ! Boundary condition only on the left side [not for periodic BC | on the right side for absorbing]
    ! 20                   ! Sponge layer width in gll points
    ! .65                  ! Att constant for sponge layer
    ! 0                    ! Use intial conditions instead of point source, and produce analytical solution
    ! 100                  ! Width of Gaussian used as intial condition


    ! -> Model files in C-Style binary floats [doubles]: Vs, Rho files are needed.
    !                                                  : For simple models, use create1Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    !    Output files in OUTPUT directory:
    !
    !####################################################################################################################

PROGRAM HDG
    
    !$ USE OMP_LIB
    IMPLICIT NONE
        
    INCLUDE "sem.h"
    INCLUDE "set.h"
    INCLUDE "dg.h"
    
    INTERFACE connectivity_matrix
        FUNCTION connectivity_matrix(N,ne)
            INTEGER, INTENT(IN)                         :: N, ne
            INTEGER connectivity_matrix(N+1,ne)
        END FUNCTION
    END INTERFACE


    TYPE(SETup)                                                :: SET
    TYPE(semProblem)                                           :: SEM
    TYPE(dgProblem)                                            :: DG
    REAL (KIND=8)                                              :: CFL, mINdist, lambdamIN
    REAL (KIND=8)                                              :: time, t_cpu_0, t_cpu_1, t_cpu, attConst, sd
                                                            
    REAL (KIND=8), dimension(:), allocatable                   :: g
    REAL (KIND=8), DIMENSION(:,:), ALLOCATABLE                 :: globalsol_s, globalsol_v
    INTEGER                                                    :: i, j, k, it, c, reclsnaps
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
    WRITE(*,*) "######## Reading parameters file #########"
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

    120 format (A,F6.1)
    140 format (A,I4,X,I4)
    160 format (A,I8)

    write(*,160)  "Polynomial order                   -> ",SET%N
    write(*,160)  "Number of elements SEM             -> ",SEM%ne
    write(*,160)  "Number of elements DG              -> ",DG%ne
    write(*,120)  "Element size                       -> ",SET%h
    write(*,160)  "Number of time steps               -> ",SET%nt
    write(*,120)  "Time step                          -> ",SET%dt
    write(*,120)  "Ricker's peak frequency            -> ",SET%f0
    write(*,140)  "Source location [nel/ngll]         -> ",SET%esrc, SET%gsrc
    write(*,160)  "Snapshot interval                  -> ",SET%isnap


    SEM%ngll  = SET%N * SEM%ne + 1                                  ! Total GLL poINts
    DG%ngll   = (SET%N + 1) * DG%ne 
    SET%Jc    = SET%h / 2                                           ! Jacobian for structured 1D mesh
    SET%Jci   = 1 / SET%Jc                                          ! Jacobian INverse
    SET%ne    = SEM%ne + DG%ne



    ALLOCATE(globalsol_s(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll+DG%ngll))
    ALLOCATE(globalsol_v(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll+DG%ngll))

    ALLOCATE(SET%xi(SET%N+1))                                       ! GLL poINts
    ALLOCATE(SET%wi(SET%N+1))                                       ! GLL Quadrature weights
    ALLOCATE(SET%Cij(SET%N+1,SET%ne))                               ! Connectivity matrix
    ALLOCATE(SET%v1D(SET%ne))                                       ! 1D velocity model IN elements
    ALLOCATE(SET%rho1D(SET%ne))                                     ! Density velocity model IN elements
    ALLOCATE(SET%lprime(SET%N+1,SET%N+1))                           ! Dervatives of Lagrange polynomials

    !##########################################
    !######### SEM &  DG SOLVER INIT ##########
    !##########################################
    ! SEM var allocation
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
    ALLOCATE(SEM%uddot(SEM%ngll), SEM%udot(SEM%ngll))
    ALLOCATE(SEM%sigma(SEM%ngll), SEM%sigmaold(SEM%ngll), SEM%sigmanew(SEM%ngll))
    ALLOCATE(SEM%tauL(SEM%ngll),SEM%tauLold(SEM%ngll),SEM%tauLnew(SEM%ngll),SEM%T(SEM%ngll))
    ALLOCATE(SEM%uold(SEM%ngll))                                            ! displacement vector at time t + dt
    ALLOCATE(SET%src(SET%nt))                                               ! Source time FUNCTION
    ALLOCATE(SEM%F(SEM%ngll))                                               ! External force
    ALLOCATE(SEM%Uout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll))               ! Snapshots
    ALLOCATE(SEM%Udotout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll), & 
             SEM%sigmaout(NINT(REAL(SET%nt/SET%isnap)),SEM%ngll))
    ALLOCATE(SEM%mu1Dgll(SEM%ngll))                                         ! Shear modulus mapped
    ALLOCATE(SEM%xgll(SEM%ngll))                                            ! Array for global mappINg
    ALLOCATE(g(SEM%ngll))
    ALLOCATE(SEM%recsigma(SET%nt))
   
    ! DG var allocation
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

    CALL zwgljd(SET%xi,SET%wi,SET%N+1,0.,0.)                                ! GettINg GLL poINts and weights
    CALL lagrangeprime(SET%N,SET%lprime)                                    ! Lagrange polynomials derivatives

    CALL readmodelfiles1D(SET%v1D, SET%rho1D, SET%ne,modnameprefix)         ! READINg model files
    CALL shapefunc(SET%N,SET%h,SEM%ne,SET%Cij,SEM%xgll)                     ! Global domaIN mappINg
    CALL shapefuncDG(SET%N,SET%h,DG%ne, DG%xgll)                            ! Global domain mapping
    CALL mapmodel(SET%N,SET%ne,SEM%ne,SET%rho1D(1:SEM%ne),SET%v1D(1:SEM%ne),SEM%rho1Dgll,SEM%v1Dgll)   ! MappINg models
    CALL ricker(SET%nt,SET%f0,SET%dt,SET%src)                               ! Source time FUNCTION

    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"
    write(*,*)"##########################################"

    mindist = SEM%xgll(2) - SEM%xgll(1)
    CFL = (SET%dt/mindist) * maxval(SET%v1D(:))

    if (CFL > .4) then
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop

    else
        print"(a14,f6.3)","CFL value is ",CFL
        print*,"Simulation is stable"
    end if


    write(*,*)"##########################################"
    write(*,*)"########## Space Sampling check ##########"
    write(*,*)"##########################################"
    lambdamin = minval(SET%v1D)/(SET%f0*2.5)

    print"(a32,f6.2)", " Elements per minimum wavelength ->", lambdamin/SET%h

    if ((lambdamin/SET%h)<1) then
        print*,"Element size is too large"
        print*,"Numerical dispersion might be present"
        print*,"Do you wish to continue anyways? Yes/no"
        read*, filecheck

        if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
                filecheck=="Y") then
            print*,"Proceeding..."

        elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
            write(*,*) "The program had been terminated"
            write(*,*) "Reduce element size or frequency as such"
            write(*,*) "you have at least 1 element per min wavelength"
            write(*,*) "This holds for N=5 to N=10 as per"
            write(*,*) "Komatitsch and Tromp 1999"
            stop
        else
            write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
            write(*,*) "The program have been terminated, please star over"
            stop
        end if
    else
        print*, "Spatial sampling as OK!"
    end if


    ! DG & SEM Matrices setup
    CALL semDefine(SET,SEM)
    CALL dgDefine(SET, DG)

    k = 0
    DO it=1,SET%nt 
        
        CALL semSolve(SET,SEM,it)

        CALL dgSolve(SET, DG, SEM, it)



        SEM%T(SET%Cij(SET%N+1,SEM%ne))       = DG%u(1,1,1)
        SEM%udot(SET%Cij(SET%N+1,SEM%ne))    = DG%u(1,1,2)
        !SEM%udotnew(SET%Cij(SET%N+1,SEM%ne)) = DG%unew(1,1,2)


        if (mod(it,SET%isnap) == 0) then
            k = k + 1
            SEM%Uout(k,:)    = SEM%u
            SEM%Udotout(k,:) = SEM%udot
            SEM%sigmaout(k,:)= SEM%sigma
            c = 1
            do i=1,DG%ne
                do j=1,SET%N+1
                    DG%sigma(k,c) = DG%u(i,j,1)
                    DG%v(k,c)     = DG%u(i,j,2)
                    c = c + 1
                end do
            end do
        end if


        if (mod(it,NINT(SET%nt/10.))==0) then
            print*,"At time sample ->",it, "/",SET%nt
        end if

    END DO
    
    globalsol_s(:,1:SEM%ngll)          =  SEM%sigmaout
    globalsol_s(:,SEM%ngll+1:DG%ngll)  =  DG%sigma

    globalsol_v(:,1:SEM%ngll)          =  SEM%Udotout
    globalsol_v(:,SEM%ngll+1:DG%ngll)  =  DG%v


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

    inquire(iolength=reclsnaps) DG%sigma
    outname = "OUTPUT/DG_snapshots_Sigma.bin"

    open(19,file=outname,access="direct",recl=reclsnaps)
    write(19,rec=1) DG%sigma
    close(19)

    outname = "OUTPUT/DG_snapshots_V.bin"
    open(20,file=outname,access="direct",recl=reclsnaps)
    write(20,rec=1) DG%v
    close(20)


    inquire(iolength=reclsnaps) globalsol_s

    outname = "OUTPUT/global_S.bin"
    open(21,file=outname,access="direct",recl=reclsnaps)
    write(21,rec=1) globalsol_s
    close(21)

    outname = "OUTPUT/global_V.bin"
    open(22,file=outname,access="direct",recl=reclsnaps)
    write(22,rec=1) globalsol_v
    close(22)


    ! Ellapse time
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





    !DEALLOCATE(DG%Z)                                                       
    !DEALLOCATE(DG%xgll)                                                  
    !DEALLOCATE(DG%Minv)                                             
    !DEALLOCATE(DG%Ke)                                                   
    !DEALLOCATE(DG%Ar,DG%Al)                                          
    !DEALLOCATE(DG%u,DG%unew)                              
    !DEALLOCATE(DG%k1,DG%k2,DG%source)
    !DEALLOCATE(DG%flux)                                   
    !DEALLOCATE(DG%sigma, DG%v)



    !DEALLOCATE(globalsol_s,globalsol_v)

    !DEALLOCATE(SET%xi)                                       
    !DEALLOCATE(SET%wi)                                       
    !DEALLOCATE(SET%Cij)                              
    !DEALLOCATE(SET%v1D)                                
    !DEALLOCATE(SET%lprime)   
    !DEALLOCATE(SET%rho1D)                                   

  
                             
    !DEALLOCATE(SEM%unew, SEM%uddotnew, SEM%udotnew,SEM%uold,SET%src,SEM%F &                                               
    !          ,SEM%Uout, SEM%Udotout, SEM%sigmaout, SEM%mu1Dgll, SEM%recsigma &
    !          ,SEM%MINv, SEM%Me, SEM%Kg, SEM%Ke, SEM%sigma, SEM%sigmaold, SEM%sigmanew &
    !          ,SEM%tauL,SEM%tauLold,SEM%tauLnew,SEM%T, SEM%xgll, g, SEM%rho1Dgll, SEM%v1Dgll, & 
    !          SEM%M, SEM%u, SEM%uddot, SEM%udot)


END PROGRAM