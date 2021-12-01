subroutine compute_flux(ne,u,N,Al,Ar,flux,SEM,SET)
    !$ use omp_lib
    implicit none
    INCLUDE "sem.h"
    INCLUDE "set.h"

    TYPE(semProblem), INTENT(INOUT)                            :: SEM
    TYPE(SETup),      INTENT(INOUT)                            :: SET
    integer, intent(in)                                        :: ne, N
    real(kind=8), intent(in)                                   :: u(ne,N+1,2), Al(ne,2,2), Ar(ne,2,2)
    real(kind=8), intent(out)                                  :: flux(ne,N+1,2)

    integer                                                    :: i
    real(kind=8)                                               :: zeros(2),boundaryflux(2)


    flux(:,:,:) = 0.
    zeros(:)    = 0.

    !##########################################
    !#### Inside the computational domain  ####
    !##########################################
    
    !$OMP PARALLEL DO PRIVATE(i) SHARED(flux,Ar,AL,u) SCHEDULE(static)
    do i=2,ne-2
        flux(i,1,:)     = MATMUL(RESHAPE(Ar(i,:,:),(/2,2/)),  RESHAPE((-u(i-1,N+1,:)),(/2/))) + &
                          MATMUL(RESHAPE(Al(i,:,:),(/2,2/)),  RESHAPE((-u(i  ,1,:)),(/2/)))
        flux(i,N+1,:)   = MATMUL(RESHAPE(Ar(i,:,:),(/2,2/)),  RESHAPE(  u(i,N+1,:),(/2/)))    + &
                          MATMUL(RESHAPE(Al(i,:,:),(/2,2/)),  RESHAPE(u(i+1,1,:),(/2/)))
    end do
    !$OMP END PARALLEL DO

    !##########################################
    !#####  At the  domain's boundaries  ######
    !##########################################

    boundaryflux(1)                                    = SEM%sigma(SET%Cij(1,SEM%ne))
    boundaryflux(2)                                    = SEM%udot(SET%Cij(1,SEM%ne))
    SEM%T(SET%Cij(N+1,SEM%ne))       = u(2,1,1)
    SEM%udot(SET%Cij(N+1,SEM%ne))    = u(2,1,2)
    SEM%udotnew(SET%Cij(N+1,SEM%ne)) = u(2,1,2)
    

    flux(1,1,:)         = MATMUL(RESHAPE(Ar(1,:,:),(/2,2/)), RESHAPE(-boundaryflux,(/2/)))    + &
                          MATMUL(RESHAPE(Al(1,:,:),(/2,2/)), RESHAPE(-u(1  ,1,:),(/2/)))
    flux(1,N+1,:)       = MATMUL(RESHAPE(Ar(1,:,:),(/2,2/)), RESHAPE( u(1,N+1,:),(/2/))) +   &
                          MATMUL(RESHAPE(Al(1,:,:),(/2,2/)), RESHAPE( u(2,1,:),(/2/)))

    flux(ne,1, :)       = MATMUL(RESHAPE(Ar(ne,:,:), (/2,2/)),RESHAPE((-u(ne-1,N+1,:)),(/2/))) +   &
                          MATMUL(RESHAPE(Al(ne,:,:), (/2,2/)),RESHAPE((-u(ne,1,:)),(/2/)))
    flux(ne,N+1,:)      = MATMUL(RESHAPE(Ar(ne,:,:), (/2,2/)) ,RESHAPE(u(ne, N+1,:),(/2/)))  + &
                          MATMUL(RESHAPE(Al(ne,:,:), (/2,2/)), zeros)



end subroutine compute_flux

