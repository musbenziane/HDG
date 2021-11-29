        !$omp parallel do private(i,k,tmp) shared(Kg,u,temp1,ngll) schedule(static)
        do i=1,ngll
            tmp = 0.0
            do k=1,ngll
                tmp = tmp + kg(i,k) * u(k)
            enddo
            temp1(i) = tmp
        enddo
        !$omp end parallel do

        !$OMP PARALLEL WORKSHARE
        temp2 = F(:) - temp1(:)
        !$OMP END PARALLEL WORKSHARE

        !$omp parallel do private(i) shared(Minv,temp2,temp3,ngll) schedule(static)
        do i=1,ngll
            temp3(i) = Minv(i,i) * temp2(i)
        enddo
        !$omp end parallel do




        do i=1,ne
            fe(:) = 0
            do j=1,N+1
                tmpf = 0
                do k=1,N+1
                    tmpf = tmpf + u(Cij(k,i)) * lprime(k,j)
                end do
                
                sigma    = mu1Dgll(Cij(j,i)) * tmpf * Jci
                tmpf1(j) = sigma * Jc * Jci 
            end do

            do j=1,N+1
                tmpf = 0
                do k=1,N+1
                    tmpf = tmpf + tmpf1(k) * lprime(j,k) * wi(k)
                end do

                fe(j) = -tmpf 
                if (i==esrc .and. j==gsrc) then
                    !$omp parallel do private(i,k,tmp) shared(Kg,u,temp1,ngll) schedule(static)
        do i=1,ngll
            tmp = 0.0
            do k=1,ngll
                tmp = tmp + kg(i,k) * unew(k)
            enddo
            temp1(i) = tmp
        enddo
        !$omp end parallel do

        !$OMP PARALLEL WORKSHARE
        temp2 = F(:) - temp1(:)
        !$OMP END PARALLEL WORKSHARE


        !$omp parallel do private(i) shared(Minv,temp2,temp3,ngll) schedule(static)
        do i=1,ngll
            temp3(i) = Minv(i,i) * temp2(i)
        enddo
        !$omp end parallel do
                end if    
            end do

            do j=1,N+1
                F(Cij(j,i)) = F(Cij(j,i)) + fe(j)
            end do
        end do


        do i=1,ngll
            tmp = 0.0
            do k=1,ngll
                tmp = tmp + kg(i,k) * unew(k)
            enddo
            temp1(i) = tmp
        enddo

            temp1(:) = 0
            do k=1, ngll
                do i=1, ngll
                    temp1(:) = temp1(:) + kg( i , k ) * unew( k  )
                end do
            end do
