subroutine semSolve(SET, SEM)

    !$ USE OMP_LIB
    implicit none
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
        INTEGER                                         :: ne,ngll
    END TYPE    

    TYPE(semProblem)                                           :: SEM
    TYPE(SETup)                                                :: SET
    REAL(kind=8)                                               :: sum, temp1(SEM%ngll), temp2(SEM%ngll), temp3(SEM%ngll),&
                                                                  tmp
    INTEGER                                                    :: i, j, k, l, el, it

    


    !##########################################
    !####### Construct the mass matrix ########
    !##########################################

    do i=1,SEM%ne
        do j=1,SET%N+1
            SEM%Me(j)           = SEM%rho1Dgll(SET%Cij(j,i)) * SET%wi(j) * SET%Jc   ! Elemental mass matrix construction
            SEM%M(SET%Cij(j,i)) =  SEM%M(SET%Cij(j,i)) + SEM%Me(j)                  ! Filling the global mass matrix
        end do
    end do

    SEM%Minv(:,:) = 0
    ! Invert mass matrix
    do i=1,SEM%ngll
        SEM%Minv(i,i) = 1 / SEM%M(i)
    end do


    !###############################################
    !####### Construct the Stiffness matrix ########
    !###############################################
    SEM%mu1Dgll(:) =  SEM%rho1Dgll(:) * SEM%v1Dgll(:)**2.          ! Shear modulus
    SEM%Kg(:,:) = 0
    !$omp parallel do private(el,i,j,k,sum) shared(SEM,SET) schedule(static)
    do el=1,SEM%ne
        do i=1,SET%N+1                                    ! Elemental stifness matrix
            do j=1,SET%N+1
                sum = 0.
                do k=1,SET%N+1
                    sum = sum + SEM%mu1Dgll(SET%Cij(k,el)) * SET%wi(k) * SET%lprime(j,k) *SET%lprime(i,k) & 
                                                           * SET%Jc * SET%Jci**2
                end do
                SEM%Ke(i,j) = sum
            end do
        end do

        do i=1,SET%N+1                                      ! Global assembly
            do j=1,SET%N+1
                SEM%Kg(SET%Cij(i,el),SET%Cij(j,el)) = SEM%Kg(SET%Cij(i,el),SET%Cij(j,el)) + SEM%Ke(i,j)
            end do
        end do
    end do
    !$omp end parallel do



    do it=1,SET%nt

        !SEM%tauL(SET%Cij(SET%N+1,200)) = SEM%recsigma(it)

        SEM%F(SET%Cij(SET%gsrc,SET%esrc)) =  SET%src(it) * SET%wi(SET%gsrc) * SET%Jc
    
        SEM%unew(:)     = SEM%u(:) + SET%dt * SEM%udot(:) + (SET%dt**2)/2 * SEM%uddot(:)

        temp1(:) = 0
        !$OMP PARALLEL DO PRIVATE(k) SHARED(SEM,SET) SCHEDULE(static) 
        do l=1, SEM%ngll
            temp1(l) =  DOT_PRODUCT(SEM%Kg( l , : ), SEM%unew(:))
        end do
        !$OMP END PARALLEL DO

        !$OMP PARALLEL WORKSHARE
        temp2 = SEM%F(:) - temp1(:) + SEM%tauL(:)
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
                    tmp = tmp + SEM%mu1Dgll(SET%Cij(i,el)) * SEM%u(SET%Cij(j,el)) * SET%lprime(j,i) * SET%Jci    
                end do
                SEM%sigma(SET%Cij(i,el)) = tmp
            end do 
        end do
        
        !recsigma(t)  = sigma(Cij(N+1,200))
       
        SEM%u     = SEM%unew
        SEM%udot  = SEM%udotnew
        SEM%uddot = SEM%uddotnew

        if (mod(it,SET%isnap) == 0) then
            k = k + 1
            SEM%Uout(k,:)    = SEM%u
            SEM%Udotout(k,:) = SEM%udot
            SEM%sigmaout(k,:)= SEM%sigma

            if (mod(it,NINT(SET%nt/10.))==0) then
                print*,"At time sample ->",it, "/",SET%nt

            end if
        end if
    end do


end subroutine