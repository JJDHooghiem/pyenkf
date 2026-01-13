! Ensemble Kalman Square Root Filter 
! Copyright (C) 2022 J.J.D. Hooghiem 

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

module enkf_core
  
  implicit none (type,external)

  ! Reference and definition of the BLAS double dot function
  real(kind=8),external :: ddot

  ! Reference and definition of the BLAS routines required 
  external :: dgemv,dger,daxpy

  private 
  
  ! enksrf is the ensemble square root filter
  ! dtrkmm is a fused matmul( kron(A1,A2), B ) with A1 and A2 lower triangular matrices

  public :: enksrf, dtrkmm

  ! cache blocking parameters for kron(a1,a2) @ B
  ! to use all ymm registers on avx2 machines
  integer, parameter :: mr = 8
  integer, parameter :: nr = 6

contains

  subroutine enksrf(nobs,nmembers,nparams,obs,mrej,rej_thr,hx,hxp,xp,r, hphr,x,rejected,assimilate)
    !
    ! This subroutine implements the Ensemble Kalman Square Root Filter 
    ! algorithm as described in 
    !    "Ensemble Data Assimilation without Perturbed Observations"
    ! by Whitaker and Hamill 2002. DOI: 10.1175/1520-0493(2002)130<1913:EDAWPO>2.0.CO;2
    ! and  
    !    "An ensemble data assimilation system to estimate CO2 surface 
    !    "fluxes from atmospheric trace gas observations" 
    ! by Peters et al. 2005. DOI: 10.1029/2005JD006157
    ! 

    ! Input variables, not changed during routine
    integer*8,intent(in)                              :: nobs     ! the amount of observations
    integer*8,intent(in)                              :: nmembers ! the amount of ensemble members used to represent the covariance
                                                                  ! strucure (sometimes called particle number) 
    integer*8,intent(in)                              :: nparams  ! number of statevector elements to be estimated
    logical,dimension(nobs), intent(in)               :: mrej     ! may reject wether an observations may be rejected
    logical,dimension(nobs), intent(in)               :: assimilate ! wether we want to assimilate this observation 
    real*8, dimension(nobs), intent(in)               :: obs      ! observation values
    real*8, dimension(nobs), intent(in)               :: R        ! observation error  
    real*8, dimension(nobs), intent(in)               :: rej_thr  ! rejection threshold for observation that may be rejected
! local variables
    integer*8                                         :: i,j       !
    real*8                                            :: n_fac     ! for the repeated factor of  1/(N-1), n=nmembers
    real*8                                            :: res,alpha ! residual, alpha factor (see Whitaker and Hamill 2002)
    real*8,dimension(nparams)                         :: PHt       ! Estimated product of PH-transpose
    real*8,dimension(nparams)                         :: KG        ! Kalman Gain 
    real*8,dimension(nobs)                            :: fac       ! factor used for updates
    real*8,dimension(nmembers)                        :: HXp_n     ! temporary row of deviations  

    !
    ! Output (updated in the loop)
    ! 
    real*8, dimension(nobs,nmembers),intent(inout)    :: HXp     ! observed deviations
    real*8, dimension(nobs),intent(inout)             :: Hx,HPHR ! simulationed mean value of the observations, HPH+R
    real*8, dimension(nparams,nmembers),intent(inout) :: Xp      ! statevecter deviations
    real*8, dimension(nparams),intent(inout)          :: x       ! statevector
    integer*8, dimension(nobs),intent(inout)          :: rejected ! returns True if observation is rejected

    ! initialize variables
    PHt=0.0 
    KG=0.0 
    fac=0.0 
    n_fac=1.0/(float(nmembers)-1.0)
       
    ! start processing observations one at a time
    do i=1,nobs
        if (assimilate(i).eqv..False.) then
                cycle 
        endif
        ! compute difference between forecast and observed:
        res           = obs(i) - Hx(i)

        ! We should be checking if we reject observations here by comparing to res(idual)
        if (mrej(i).eqv..True.) then 
                if ( abs(res) > rej_thr(i) * sqrt(R(i) )) then 
                        rejected(i)=1
                        cycle
                endif
        endif 

        ! Local copy of the deviations row HXp(i,:) will be updated later
        HXp_n=HXp(i,:)

        ! calculate PHt using blas/lapac DGEMV 
        ! Computes PHt = 1/(N-1) HXp_n . Xp 
        call DGEMV('N',nparams,nmembers,n_fac,Xp,nparams,HXp_n,1,0.0,PHt,1)

        ! estimates HPH+R=1/(N-1) HXp_n . HXp_n + R  , using DDOT. 
        HPHR(i) = n_fac * DDOT(nmembers,HXp_n,1,HXp_n,1) + R(i)

        ! Compute the Kalman Gain
        KG = PHt/HPHR(i)

        ! alpha factor
        alpha=1.0 / (1.0 + sqrt( R(i) / HPHR(i) ))

        ! We now have the info to compute the updated version of
        ! x, Xp, HXp, and Hx  
        ! some of this stuf is not final yet 

        ! could it be more effiecient to do this and the part of fac in order to update Hx?
        ! call DGEMV('N',nobs,nmembers,res*n_fac*(1/HPHR(i)),HXp,nmembers,HXp(i,:),1,1.0,Hx,1)

        ! Approximate HK, and store in fac
        call DGEMV( 'N' , nobs , nmembers , n_fac * (1.0 / HPHR(i)),HXp,nobs,HXp_n,1,0.0,fac,1)

        ! update Hx = Hx + HK * residual 
        call daxpy(nobs,res,fac,1,Hx,1)

        ! update x = x  + res * KG 
        call daxpy(nparams,res,KG,1,x,1)

        ! update deviations HX = HX - HK HXp(i) using outer product
        call DGER(nobs,nmembers,-1.0*alpha,fac,1,HXp_n,1,HXp,nobs)

        ! update Xp 
        call DGER(nparams,nmembers,-1.0*alpha,KG,1,HXp_n,1,Xp,nparams)

        ! update deviations HX = HX - HK HXp(i) using outer product in loop
        ! do j=1,nmembers
        !         ! update Xp
        !         call daxpy(nparams,-1.0*alpha*HXp_n(j),KG,1,Xp(:,j),1)
        !         ! update HXp for column j
        !         call daxpy(nobs,-1.0*alpha*HXp_n(j) ,fac,1,HXp(:,j),1)
        ! enddo ! loop over nmembers 

    ! end loop over observations
    enddo 

  end subroutine  enksrf
  ! reference implementation that is much slower than the block version
   subroutine dtrkmm_ref(N, M, K, NM, A1, A2, B, C)
      !
      ! Computes the
      !       matmul( kron(A1,A2), B )
      !       where A(N,N) and B(M,M) are cholesky decompositions
      !       B is a matrix with dimension (NM,K)
      !       and stores the result in C(NM,K)
      !
      integer*8, intent(in) :: N, M, K, NM
      real*8, dimension(N, N), intent(in) :: A1
      real*8, dimension(M, M), intent(in) :: A2
      real*8, dimension(NM, K), intent(in) :: B
      real*8, dimension(NM, K), intent(inout) :: C
      integer*8            :: i, j, h, l, g

      do i = 1, N
         do j = 1, i
            do h = 1, M
               do l = 1, h
                  do g = 1, K
                     C(h + (i - 1)*M, g) = C(h + (i - 1)*M, g) + A1(i, j)*A2(h, l)*B(l + (j - 1)*M, g)
                  end do
               end do
            end do
         end do
      end do

   end subroutine dtrkmm_ref

   subroutine dtrkmm(N, M, K, NM, A1, A2, B, C)
      !
      ! Computes the
      !       matmul( kron(A1,A2), B )
      !       where A1(N,N) and A2(M,M) are cholesky decompositions
      !       B is a matrix with dimension (NM,K)
      !       and stores the result in C(NM,K)
      !
      integer, intent(in) :: N, M, K, NM
      real*8, dimension(N, N), intent(in) :: A1
      real*8, dimension(M, M), intent(in) :: A2
      real*8, dimension(NM, K), intent(in) :: B
      real*8, dimension(NM, K), intent(inout) :: C
      integer            :: i, j, h, l, g, gg, hh, jnrmax, imrmax, imr, jnr, iblnr, jblnr
      real*8 :: aa
      real*8, dimension(mr, M*(M/mr + 1)) :: ablock
      real*8, dimension(nr, M*(K/nr + 1)) :: bblock


      !! repack A2 so that it can be accessed contiguously later on...
      !! for large A2 we should consider implementing a blocking along M
      !! or a pivoting scheme and reversing A2 and A1  
      do h = 1, M, mr 
         imrmax = min(mr, M - h + 1)
         iblnr = h/mr
         do l = 1, M
            do imr = 1, imrmax
               ablock(imr, l + iblnr*M) = A2(h - 1 + imr, l)
            end do
            do imr = imrmax + 1, mr
               ablock(imr, l + iblnr*M) = 0.0
            end do
         end do
      end do

      do i = 1, N
         do j = 1, i
            ! skip A1 values that are 0 as it is 
            ! factor for all that follows
            if (A1(i, j) .eq. 0.0) cycle
            ! pack a panel of B and transpose for contiguous access 

            do g = 1, K, nr
               jnrmax = min(nr, K - g + 1)
               jblnr = g/nr
               do l = 1, M
                  do jnr = 1, jnrmax
                     bblock(jnr, l + jblnr*M) = B(l + (j - 1)*M, jnr + g - 1)
                  end do
                  do jnr = jnrmax + 1, nr
                     bblock(jnr, l + jblnr*M) = 0.0
                  end do
               end do
            end do
            ! original mloop
            !do l=1,M ! this loop over M signals columns of A2 and rows of B and is effectively a

            !$omp parallel do private(jnrmax,jblnr,g,h,imrmax,iblnr)
            do g = 1, K, nr
               jnrmax = min(nr, K - g + 1)
               jblnr = g/nr
               do h = 1, M, mr !needed to create ablocks
                  imrmax = min(mr, M - h + 1)
                  iblnr = h/mr
                  ! within the kernel well loop over ablock(mr,M) values A2
                  ! and bblock(nr,M)
                  call mm_kernel_8x6(C((i - 1)*M + h:, g:), M,  &
                  ablock(:, iblnr*M + 1:(iblnr + 1)*M - 1),     &
                  bblock(:, jblnr*M + 1:(jblnr + 1)*M - 1),     & 
                  imrmax, jnrmax, A1(i, j))
                  ! reference for what happens in the kernel in terms of original A1 A2 and B
                  !do gg=1,jnrmax
                  !do hh=1,imrmax
                  !C(hh+(i-1)*M+h-1,gg+g-1) =C(hh+(i-1)*M+h-1,gg+g-1)+  A1(i,j)* A2(hh+h-1,l) *B(l+(j-1)*M,gg+g-1)
                  !enddo
                  !enddo
                  !enddo
               end do
            end do
            !$omp end parallel do
            !end do M loop in original version, i.e. without kernel and ablock and bblock
         end do
      end do

   end subroutine dtrkmm

   subroutine mm_kernel_8x6(Cout, pmax, ablock, bblock, imrmax, jnrmax, alpha)
      ! compute kernel that targets avx2 ymm registers
      ! imrmax/jmrmax, max copyout indices for edge cases
      ! alpha to ultiply with a constant factor
      integer, intent(in) :: pmax, imrmax, jnrmax
      real*8, intent(in) :: alpha
      real*8, contiguous, dimension(:, :) :: ablock
      real*8, contiguous, dimension(:, :) :: bblock
      real*8, intent(inout), dimension(:, :) ::  Cout
      real*8, dimension(mr, nr) :: C ! 12 registers ymm
      integer ::jjr, iir, pr

      ! initialize
      C = 0.0

      ! unroll all mr nr loops
      do pr = 1, pmax
         C(:, 1) = C(:, 1) + bblock(1, pr)*ablock(:, pr)
         C(:, 2) = C(:, 2) + bblock(2, pr)*ablock(:, pr)
         C(:, 3) = C(:, 3) + bblock(3, pr)*ablock(:, pr)
         C(:, 4) = C(:, 4) + bblock(4, pr)*ablock(:, pr)
         C(:, 5) = C(:, 5) + bblock(5, pr)*ablock(:, pr)
         C(:,  6) = C(:, 6) + bblock( 6, pr)*ablock(:, pr)
         !! on avx512 we might want to use the zmm registers (512 bits/ 32 registers)
         ! C(:,  7) = C(:, 7) + bblock( 7, pr)*ablock(:, pr)
         ! C(:,  8) = C(:, 8) + bblock( 8, pr)*ablock(:, pr)
         ! C(:,  9) = C(:, 9) + bblock( 9, pr)*ablock(:, pr)
         ! C(:,  9) = C(:, 9) + bblock( 9, pr)*ablock(:, pr)
         ! C(:, 10) = C(:,10) + bblock(10, pr)*ablock(:, pr)
         ! C(:, 11) = C(:,11) + bblock(11, pr)*ablock(:, pr)
         ! C(:, 12) = C(:,12) + bblock(12, pr)*ablock(:, pr)
      end do

      ! multiply 
      ! C = alpha*C
      ! copyout only those result required and multiply with constant (FMA)
      ! in the future consider smaller kernels to handle edge cases
      do jjr = 1, jnrmax
      do iir = 1, imrmax
         Cout(iir, jjr) = Cout(iir, jjr) + alpha*C(iir, jjr)
      end do
      end do

   end subroutine mm_kernel_8x6

end module enkf_core
