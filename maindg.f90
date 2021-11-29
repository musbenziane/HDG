program DG1D
    !####################################################################################################################
    ! DG1D: Noda Discontinuous Galerkin Method in 1D for the Elastic Wave Equation with Regular Mesh
    ! NOTE: This program is yet to be tested against analytical solutions.
    ! 
    ! Language: Fortran 90, with parallel impelementation using OpenMP API
    ! 
    ! Sources:  Igel 2017, Leveque 2002, Hesthaven & Warburton 2008
    ! 
    ! The code used for arbitrary polynomial order to get GLL points and weights was created by M.I.T departement of engineering. 
    ! Link is hereafter
    ! https://geodynamics.org/cig/doxygen/release/specfem3d/    | file name: gll_library.f90
    ! 
    ! This is part of the Numerical Modelling Workshop.
    ! 
    ! Supervisor: Pr. Emmanual Chaljub
    ! Author    : Mus Benziane
    !
    ! Input file: [Example]
    ! testing              ! model name prefix (model names: prefix_vp, prefix_rho)
    ! 4000                 ! xmax
    ! 6                    ! Polynomial order
    ! 800                  ! Number of elements
    ! 5                    ! Element size
    ! 3.                   ! Wavelet's peak frequency
    ! 0.000040             ! Time step
    ! 50000                ! Number of time steps
    ! 3                    ! Source location - element number
    ! 1                    ! Source location - collocation point
    ! 100                  ! Snapshot interval
    ! 3                    ! [1/2/3] 1: Free surface, 2: Rigid wall, 3: Periodic
    ! 1                    ! Boundary condition only on the left side [not for periodic BC]
    ! 1                    ! Initial conditions instead of source injection [will produce analytical solution]
    ! 200                  ! Gaussian width for IC and for analytical solution
    ! 
    ! -> Model files in C-Style binary floats [doubles]: Vs, Rho files are needed.
    !                                                  : For simple models, use create1Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    !    Output files in OUTPUT directory:
    !
    !####################################################################################################################


!$ use omp_lib
implicit none
    real(kind=8)                                  :: xmax, h, f0, dt, CFL, mindist, Ja, Jai, time
    real(kind=8)                                  :: t_cpu_0, t_cpu_1, t_cpu, sd
    real(kind=8), dimension(:),     allocatable   :: xi, wi, v1D, rho1D, src, mu, Z, xglls
    real(kind=8), dimension(:,:),   allocatable   :: lprime, xgll, Minv, Ke, sigma, v, sgreens, vgreens, asgreens, & 
                                                     sas,sav
    real(kind=8), dimension(:,:,:), allocatable   :: Al, Ar, u, unew, k1, k2, flux,  as, av, source
    integer                                       :: N, ne, nt, esrc, gsrc, isnap, ngll, i, j, k, it, el,c, &
                                                     bc, sbc, IC, reclsnaps, ir, t0, t1, ircv, nrcv, is_rcv
    integer, dimension(:,:), allocatable          :: Cij
    character(len=40)                             :: filename, filecheck, outname_sigma, outname_v,  &
                                                     modnameprefix
    !$ integer                                    :: n_workers
    logical                                       :: OMPcheck = .false.


    write(*,*) "##########################################"
    write(*,*) "############### OpenMP     ###############"
    !$ OMPcheck = .true.
    if (OMPcheck) then
        !$OMP PARALLEL
        !$ n_workers = OMP_GET_NUM_THREADS()
        !$OMP END PARALLEL
        !$ print '(3X,"Number of workers ->  ",i2)',n_workers
    else
        write(*,*) "Program has been compiled without OpenMP; Time marching will run in serial"
    end if


    call cpu_time(t_cpu_0)
    call system_clock(count=t0, count_rate=ir)

    filename          = "parameters.in"
    outname_sigma     = "OUTPUT/snapshots_sigma.bin"
    outname_v         = "OUTPUT/snapshots_v.bin"

    write(*,*) "##########################################"
    write(*,*) "######## Reading parameters file #########"
    write(*,*) "##########################################"

    print*,"Is the parameters input file (parameters.in) [Yes/no]"
    read(*,*) filecheck

    if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") then
        write(*,*) "Reading simulation parameters..."

    elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
            filecheck=="N") then
        write(*,*) "Enter simulation parameters text file name with extension"
        write(*,*) "40 characters max"
        read(*,*) filename

    else
        write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        write(*,*) "The program have been terminated, please star over"
        stop
    end if

    open (2, file=filename, status = 'old')
    read(2,*) modnameprefix
    read(2,*) xmax
    read(2,*) N
    read(2,*) ne
    read(2,*) h
    read(2,*) f0
    read(2,*) dt
    read(2,*) nt
    read(2,*) esrc
    read(2,*) gsrc
    read(2,*) isnap
    read(2,*) is_rcv
    read(2,*) ircv
    read(2,*) bc
    read(2,*) sbc
    read(2,*) IC
    read(2,*) sd
    close(2)

    120 format (A,F6.1)
    140 format (A,I4,X,I4)
    160 format (A,I8)

    write(*,120)  "Maximum distance zmax, xmax        -> ",xmax
    write(*,160)  "Polynomial order                   -> ",N
    write(*,160)  "Number of elements                 -> ",ne
    write(*,120)  "Element size                       -> ",h
    write(*,160)  "Number of time steps               -> ",nt
    write(*,120)  "Time step                          -> ",dt
    write(*,120)  "Ricker's peak frequency            -> ",f0
    write(*,140)  "Source location [nel/ngll]         -> ",esrc, gsrc
    write(*,160)  "Snapshot interval                  -> ",isnap
    write(*,160)  "Free Surface BC [1/2/3]==[FS/R/P]  -> ",bc


    !##########################################
    !#####    Matrices allocation         #####
    !##########################################

    ngll = (N + 1) * ne                 ! Number of GLL points
    nrcv = NINT(REAL(ne/ircv))

    allocate(xi(N+1))                                                                  ! GLL points
    allocate(wi(N+1))                                                                  ! GLL Quadrature weights
    allocate(v1D(ne))                                                                  ! 1D velocity model in elements
    allocate(rho1D(ne))                                                                ! Density velocity model in elements
    allocate(mu(ne))                                                                   ! Shear modulus mapped
    allocate(Z(ne))                                                                    ! Impepdences   
    allocate(xgll(ne,N+1),xglls(ngll))                                                 ! Array for global mapping
    allocate(lprime(N+1,N+1))                                                          ! Dervatives of Lagrange polynomials
    allocate(src(nt))                                                                  ! Source time function
    allocate(Cij(N+1,ne))                                                              ! Connectivity matrix
    allocate(Minv(N+1,N+1))                                                            ! Elemental mass matrix
    allocate(Ke(N+1,N+1))                                                              ! Elemental stifness matrix
    allocate(Ar(ne,2,2),Al(ne,2,2))                                                    ! Wave PDE Coefficient Matrix
    allocate(u(ne,N+1,2),unew(ne,N+1,2))                                               ! Solution fields
    allocate(k1(ne,N+1,2),k2(ne,N+1,2),source(ne,N+1,2))                                  ! RK
    allocate(flux(ne,N+1,2))                                                           ! Flux matrix
    allocate(sigma(NINT(nt/REAL(isnap)),ngll),  v(NINT(nt/REAL(isnap)),ngll),&
                as(ne,N+1,NINT(nt/REAL(isnap))), av(ne,N+1,NINT(nt/REAL(isnap))))      ! Stretched solution fields
    allocate(sgreens(nt,nrcv),vgreens(nt,nrcv))
    allocate(asgreens(nt,nrcv))
    allocate(sas(ngll,NINT(nt/REAL(isnap))),sav(ngll,NINT(nt/REAL(isnap))))



    !##########################################
    !#####      Calling subroutines       #####
    !##########################################

    call readmodelfiles1D(v1D, rho1D, ne, modnameprefix)          ! Reading model files
    call shapefunc(N,h,ne, xgll)                                  ! Global domain mapping
    call lagrangeprime(N,lprime)                                  ! Lagrange polynomials derivatives
    call ricker(nt,f0,dt,src)                                     ! Source time function
    call connectivity_matrix(N,ne,Cij) 
    call zwgljd(xi,wi,N+1,0.,0.)                                  ! Getting GLL points and weights

    Ja   = h / 2.                ! Jacobian: Regular 1D mesh is used
    Jai  = 1 / ja                ! Jacobian's inverse

    mu = (v1D**2) * rho1D        ! Shear modulus

    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"

    mindist = xgll(1,2) - xgll(1,1)
    CFL = (dt / mindist) * maxval(v1D(:))
    if (CFL > .19) then
        print"( a14,f6.3)"," Courant number is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop
    else
        print"(a14,f6.3)","Courant number is ",CFL
        print*,"Simulation is stable"
    end if
    write(*,*)"##########################################"


    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,N+1
        Minv(i,i) = 1. / (wi(i) * Ja)
    end do

    !##########################################
    !##### Construct the stiffness matrix #####
    !##########################################

    do i=1,N+1
        do j=1,N+1
               Ke(i,j) = wi(j) * lprime(i,j)
        end do
    end do


    !##########################################
    !##### Construct the Flux matrices    #####
    !##########################################

    Z = rho1D * v1D

    !$OMP PARALLEL DO PRIVATE(i) SHARED(Ar,Al,v1D,Z) SCHEDULE(static)
    do i=1,ne-2
        Ar(i,1,1) =  .5 * v1D(i)
        Ar(i,1,2) = -.5 * Z(i) * v1D(i)
        Ar(i,2,1) = -.5 * v1D(i) / Z(i)
        Ar(i,2,2) =  .5 * v1D(i)

        Al(i,1,1) = -.5 * v1D(i)
        Al(i,1,2) = -.5 * Z(i) * v1D(i)
        Al(i,2,1) = -.5 * v1D(i) / Z(i)
        Al(i,2,2) = -.5 * v1D(i);
    end do
    !$OMP END PARALLEL DO



    write(*,*) "##########################################"
    write(*,*) "########### Begin time  loop  ############"

    sgreens(:,:) = 0
    vgreens(:,:) = 0
    source(:,:,:)   = 0
    k = 0

    if (IC==1) then
        do el=1,ne
            do j=1,N+1
                u(el,j,2) = exp(-(1/sd)**2*((xgll(el,j)-xgll(esrc,gsrc)))**2)
            end do
        end do
    end if

    do it=1,nt
        if (IC .ne. 1) then
            !u(esrc,gsrc,2)  = src(it)                    ! source injection - stress component
            source(esrc,gsrc,2) =  src(it) * wi(gsrc) * Ja
        end if

        call compute_flux(ne,u,N,Al,Ar,flux)

        !$OMP PARALLEL DO PRIVATE(el) SHARED(k1,Minv,mu,rho1D,ke,u,flux) SCHEDULE(static)
        do el=2,ne-1
            k1(el,:,1)   = MATMUL(Minv,(-mu(el) * MATMUL(Ke, u(el,:,2))) - flux(el,:,1) + source(el,:,1))
            k1(el,:,2)   = MATMUL(Minv,(-1. / rho1D(el)) * MATMUL(Ke,u(el,:,1)) - flux(el,:,2) + source(el,:,2))
        end do
        !$OMP END PARALLEL DO 

        !$OMP PARALLEL DO PRIVATE(el) SHARED(unew,dt,Minv,mu,rho1D,ke,u,flux) SCHEDULE(static)
        do el=2,ne-2
            unew(el,:,1) = dt * MATMUL(Minv,(-mu(el) * MATMUL(Ke,u(el,:,2))) - &
                                       flux(el,:,1) + source(el,:,1)) + u(el,:,1)
            unew(el,:,2) = dt * MATMUL(Minv,(-1. / rho1D(el)) * MATMUL(Ke,u(el,:,1)) - &
                                       flux(el,:,2) + source(el,:,2)) + u(el,:,2)
        end do
        !$OMP END PARALLEL DO

        call compute_flux(ne,unew,N,Al,Ar,flux)

        !$OMP PARALLEL DO PRIVATE(el) SHARED(k2,dt,Minv,mu,rho1D,ke,u,flux) SCHEDULE(static)
        do el=2,ne-2
            k2(el,:,1)   = MATMUL(Minv,(-mu(el) * MATMUL(Ke, unew(el,:,2))) - flux(el,:,1) + source(el,:,1))
            k2(el,:,2)   = MATMUL(Minv,(-1. / rho1D(el)) * MATMUL(Ke,unew(el,:,1)) - flux(el,:,2) + source(el,:,2))
        end do
        !$OMP END PARALLEL DO

        unew = u + .5 * dt * (k1 + k2)
        u    = unew;

        !##########################################
        !##### Boundary Conditions            #####
        !##### 1: Free surface                #####
        !##### 2: Rigid wall                  #####
        !##### 3: Periodic                    #####
        !##########################################

        if (bc .eq. 1) then ! Free surface BC
            u(1,1:N+1,1)    = -u(4,1:N+1,1)
            u(2,1:N+1,1)    = -u(3,1:N+1,1)
            u(1,1:N+1,2)    =  u(4,1:N+1,2)
            u(2,1:N+1,2)    =  u(3,1:N+1,2)

            if (sbc .ne. 1) then
                u(ne-1,1:N+1,1) = -u(ne-2,1:N+1,1)
                u(ne,1:N+1,1)   = -u(ne-3,1:N+1,1)
                u(ne-1,1:N+1,2) =  u(ne-2,1:N+1,2)
                u(ne,1:N+1,2)   =  u(ne-3,1:N+1,2)
            end if
        end if

        if (bc .eq. 2) then ! Rigid BC
            u(1,1:N+1,1)    =   u(4,1:N+1,1)
            u(2,1:N+1,1)    =   u(3,1:N+1,1)
            u(1,1:N+1,2)    =  -u(4,1:N+1,2)
            u(2,1:N+1,2)    =  -u(3,1:N+1,2)

            if (sbc .ne. 1) then
                u(ne-1,1:N+1,1) =   u(ne-2,1:N+1,1)
                u(ne,1:N+1,1)   =   u(ne-3,1:N+1,1)
                u(ne-1,1:N+1,2) =  -u(ne-2,1:N+1,2)
                u(ne,1:N+1,2)   =  -u(ne-3,1:N+1,2)
            end if
        end if

        if (bc .eq. 3) then ! Periodic
            u(1,1:N+1,1)    =   u(ne-3,1:N+1,1)
            u(2,1:N+1,1)    =   u(ne-2,1:N+1,1)
            u(1,1:N+1,2)    =   u(ne-3,1:N+1,2)
            u(2,1:N+1,2)    =   u(ne-2,1:N+1,2)

            u(ne-1,1:N+1,1) =   u(3,1:N+1,1)
            u(ne,1:N+1,1)   =   u(4,1:N+1,1)
            u(ne-1,1:N+1,2) =   u(3,1:N+1,2)
            u(ne,1:N+1,2)   =   u(4,1:N+1,2)
        end if
        
        !##########################################
        !### Extract solution at rcv position   ###
        !##########################################

        if (is_rcv .eq. 1) then
            do i=1,nrcv
                sgreens(it,i) = u(i*ircv,CEILING(REAL((N+1)/2)),1) 
                vgreens(it,i) = u(i*ircv,CEILING(REAL((N+1)/2)),2) 
            end do
        end if

        !##########################################
        !##### Stretch solution matrices      #####
        !##########################################

        if ( mod(it,isnap) == 0) then
            k = k + 1
            c = 1
            do i=1,ne
                do j=1,N+1
                    sigma(k,Cij(j,i)) = u(i,j,1)
                    v(k,Cij(j,i))     = u(i,j,2)
                    c = c + 1
                end do
            end do
        end if

        if (mod(it,NINT(nt/10.)) == 0) then
            print*, "########### At time sample ->",it, "/",nt
        end if

    end do



    !##########################################
    !##### Analytical Solution            #####
    !##########################################
    if (IC==1) then
        k = 1
        do it=1,nt,isnap
            !$OMP PARALLEL DO PRIVATE(i,j)  SCHEDULE(static)
            do i=1,ne
                do j=1,N+1
                    av(i,j,k) = 1./2.*(exp(-1./sd**2 * (xgll(i,j)-xgll(esrc,gsrc) + MAXVAL(v1D)*it*dt)**2) + &
                                       exp(-1./sd**2 * (xgll(i,j)-xgll(esrc,gsrc) - MAXVAL(v1D)*it*dt)**2))

                    as(i,j,k) = 1./(2.*MAXVAL(Z))*(exp(-1./sd**2 * (xgll(i,j)-xgll(esrc,gsrc) + MAXVAL(v1D)*it*dt)**2) - &
                                                   exp(-1./sd**2 * (xgll(i,j)-xgll(esrc,gsrc) - MAXVAL(v1D)*it*dt)**2))
                end do
            enddo
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(i,j)  SCHEDULE(static)
            do i=1,ne
                do j=1,N+1
                    sas(Cij(j,i),k) = as(i,j,k)
                    sav(Cij(j,i),k) = av(i,j,k)
                end do
            end do
            !$OMP END PARALLEL DO

            k = k + 1
        end do
    end if



    write(*,*) "##########################################"
    write(*,*)  "Number of snapshots recorded          -> ",k

    write(*,*) "##########################################"
    write(*,*) "######### Write solution binaries ########"
    write(*,*) "##########################################"

    inquire(iolength=reclsnaps) sigma

    open(3,file=outname_sigma,access="direct",recl=reclsnaps)
    write(3,rec=1) sigma
    close(3)

    open(4,file=outname_v,access="direct",recl=reclsnaps)
    write(4,rec=1) v
    close(4)

    open(15,file="OUTPUT/source.bin",access="direct",recl=nt*8)
    write(15,rec=1) src
    close(15)

    if (is_rcv .eq. 1) then
        reclsnaps = 0 
        inquire(iolength=reclsnaps) sgreens
        open(16,file="OUTPUT/seism_sigma.bin",access="direct",recl=reclsnaps)
        write(16,rec=1) sgreens
        close(16)
        
        open(17,file="OUTPUT/seism_v.bin",access="direct",recl=reclsnaps)
        write(17,rec=1) vgreens
        close(17)
    end if

    if (IC .eq. 1) then
        reclsnaps = 0 
        inquire(iolength=reclsnaps) sas

        open(18,file="OUTPUT/snapshots_sigma_analytical.bin",access="direct",recl=reclsnaps)
        write(18,rec=1) sas
        close(18)

        open(19,file="OUTPUT/snapshots_v_analytical.bin",access="direct",recl=reclsnaps)
        write(19,rec=1) sav
        close(19)
    end if

    reclsnaps = 0 
    inquire(iolength=reclsnaps) xglls
    
    open(20,file="OUTPUT/xcoord_global_mapped.bin",access="direct",recl=reclsnaps)
    write(20,rec=1) xgll
    close(20)
    
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


    deallocate(xi)
    deallocate(wi)
    deallocate(v1D)
    deallocate(rho1D)
    deallocate(mu)
    deallocate(Z)
    deallocate(xgll)
    deallocate(lprime)
    deallocate(src)
    deallocate(Minv)
    deallocate(Ke)
    deallocate(Ar,Al)
    deallocate(u,unew)
    deallocate(k1,k2)
    deallocate(flux)
    deallocate(sigma,v,as,av,sgreens,vgreens)
end program