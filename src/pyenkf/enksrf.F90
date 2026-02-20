! Ensembln Kalman Square Root Filter
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

   !! timer function
   use omp_lib

   implicit none(type, external)

   ! Reference and definition of the BLAS double dot function
   real(kind=8), external :: ddot

   !
   ! Reference and definition of the BLAS routines required
   external :: dgemv, dger, daxpy, dgemm,openblas_set_num_threads,openblas_get_num_threads,dsymv,dpotri,dpotrf,dsyrk

   private

   ! enksrf is the ensemble square root filter
   ! dtrkmm is a fused matmul( kron(A1,A2), B ) with A1 and A2 lower triangular matrices

   public :: enksrf, enksrf_ref, dtrkmm

   ! cache blocking parameters for kron(a1,a2) @ B
   ! to use all ymm registers on avx2 machines
   integer, parameter :: mr = 8
   integer, parameter :: nr = 6

   ! Blocksizes used in openmp parallelization of the ensemble kalman filter
   ! system tuning required! depends on cache sizes
   integer*8, parameter :: blsize_l1d = 64*10
   integer*8, parameter :: nrblocksl2 = 1
   integer*8, parameter :: blsize_l2d = nrblocksl2*blsize_l1d

contains

   subroutine enksrf_ref(nobs, nmembers, nparams, obs, mrej, rej_thr, hx, hxp, xp, r, hphr, x, rejected, assimilate)
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
      integer*8, intent(in)                              :: nobs     ! the amount of observations
      integer*8, intent(in)                              :: nmembers ! the amount of ensemble members used to represent the covariance
      ! strucure (sometimes called particle number)
      integer*8, intent(in)                              :: nparams  ! number of statevector elements to be estimated
      logical, dimension(nobs), intent(in)               :: mrej     ! may reject wether an observations may be rejected
      logical, dimension(nobs), intent(in)               :: assimilate ! wether we want to assimilate this observation
      real*8, dimension(nobs), intent(in)               :: obs      ! observation values
      real*8, dimension(nobs), intent(in)               :: R        ! observation error
      real*8, dimension(nobs), intent(in)               :: rej_thr  ! rejection threshold for observation that may be rejected
! local variables
      integer*8                                         :: i, j       !
      real*8                                            :: n_fac     ! for the repeated factor of  1/(N-1), n=nmembers
      real*8                                            :: res, alpha ! residual, alpha factor (see Whitaker and Hamill 2002)
      real*8, dimension(nparams)                         :: PHt       ! Estimated product of PH-transpose
      real*8, dimension(nparams)                         :: KG        ! Kalman Gain
      real*8, dimension(nobs)                            :: fac       ! factor used for updates
      real*8, dimension(nmembers)                        :: HXp_n     ! temporary row of deviations

      !
      ! Output (updated in the loop)
      !
      real*8, dimension(nobs, nmembers), intent(inout)    :: HXp     ! observed deviations
      real*8, dimension(nobs), intent(inout)             :: Hx, HPHR ! simulationed mean value of the observations, HPH+R
      real*8, dimension(nparams, nmembers), intent(inout) :: Xp      ! statevecter deviations
      real*8, dimension(nparams), intent(inout)          :: x       ! statevector
      integer*8, dimension(nobs), intent(inout)          :: rejected ! returns True if observation is rejected

      ! initialize variables
      PHt = 0.0
      KG = 0.0
      fac = 0.0
      n_fac = 1.0/(float(nmembers) - 1.0)

      ! start processing observations one at a time
      do i = 1, nobs
         if (assimilate(i) .eqv. .False.) then
            cycle
         end if
         ! compute difference between forecast and observed:
         res = obs(i) - Hx(i)

         ! We should be checking if we reject observations here by comparing to res(idual)
         if (mrej(i) .eqv. .True.) then
            if (abs(res) > rej_thr(i)*sqrt(R(i))) then
               rejected(i) = 1
               cycle
            end if
         end if

         ! Local copy of the deviations row HXp(i,:) will be updated later
         HXp_n = HXp(i, :)

         ! calculate PHt using blas/lapac DGEMV
         ! Computes PHt = 1/(N-1) HXp_n . Xp
         call DGEMV('N', nparams, nmembers, n_fac, Xp, nparams, HXp_n, 1, 0.0, PHt, 1)

         ! estimates HPH+R=1/(N-1) HXp_n . HXp_n + R  , using DDOT.
         HPHR(i) = n_fac*DDOT(nmembers, HXp_n, 1, HXp_n, 1) + R(i)

         ! Compute the Kalman Gain
         KG = PHt/HPHR(i)

         ! alpha factor
         alpha = 1.0/(1.0 + sqrt(R(i)/HPHR(i)))

         ! We now have the info to compute the updated version of
         ! x, Xp, HXp, and Hx
         ! some of this stuf is not final yet

         ! could it be more effiecient to do this and the part of fac in order to update Hx?
         ! call DGEMV('N',nobs,nmembers,res*n_fac*(1/HPHR(i)),HXp,nmembers,HXp(i,:),1,1.0,Hx,1)

         ! Approximate HK, and store in fac
         call DGEMV('N', nobs, nmembers, n_fac*(1.0/HPHR(i)), HXp, nobs, HXp_n, 1, 0.0, fac, 1)

         ! update Hx = Hx + HK * residual
         call daxpy(nobs, res, fac, 1, Hx, 1)

         ! update x = x  + res * KG
         call daxpy(nparams, res, KG, 1, x, 1)

         ! update deviations HX = HX - HK HXp(i) using outer product
         ! call DGER(nobs, nmembers, -1.0*alpha, fac, 1, HXp_n, 1, HXp, nobs)

         ! update Xp
         ! call DGER(nparams, nmembers, -1.0*alpha, KG, 1, HXp_n, 1, Xp, nparams)

         ! update deviations HX = HX - HK HXp(i) using outer product in loop
         do j=1,nmembers
                 ! update Xp
                 call daxpy(nparams,-1.0*alpha*HXp_n(j),KG,1,Xp(:,j),1)
                 ! update HXp for column j
                 call daxpy(nobs,-1.0*alpha*HXp_n(j) ,fac,1,HXp(:,j),1)
         enddo ! loop over nmembers

         ! end loop over observations
      end do

   end subroutine enksrf_ref
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
                  call mm_kernel_8x6(C((i - 1)*M + h:, g:), M, &
                                     ablock(:, iblnr*M + 1:(iblnr + 1)*M - 1), &
                                     bblock(:, jblnr*M + 1:(jblnr + 1)*M - 1), &
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
         C(:, 6) = C(:, 6) + bblock(6, pr)*ablock(:, pr)
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

   subroutine enksrf_kernel(nmembers, nparams, meanvec, devvec, devmat, res, n_fac, R)
      !
      ! Kernel for updating a vector and deviation matrix in enkf routines
      ! either called as
      !         Hx = Hx + HK * res
      !         HX_prime = HX_prime + alpha*HKHm
      ! or as:
      !         x = x + HK * res
      !         X_prime = X_prime + alpha*KHm
      !
      !
      !

      ! Input variables, not changed during routine
      integer*8, intent(in)                                  :: nmembers         ! the amount of ensemble members used to represent the covariance
      integer*8, intent(in)                                  :: nparams          ! number of statevector elements to be estimated
      real*8, intent(in)                                      :: res, R,n_fac ! residual,R,and (1/N-1) 

      ! Output (updated in the loop)
      real*8, dimension(blsize_l1d, nmembers), intent(inout) :: devmat           ! deviation matrix (state|observations)
      real*8, dimension(blsize_l1d), intent(inout)           :: meanvec          ! mean vector (state or observations)

      ! local
      real*8, dimension(nmembers)                            :: devvec           ! devs to analyse (obs)
      real*8, dimension(blsize_l1d)                          :: K                ! Either HK or K
      integer                                                :: j
      real*8                                                 :: tstart,n_fac2, alpha,HPHR !alpha factor (see Whitaker and Hamill 2002)


      ! estimates HPH+R=1/(N-1) HXp_n . HXp_n + R  , using DDOT.
      HPHR = n_fac*DDOT(nmembers, devvec, 1, devvec, 1) + R

      ! alpha factor
      alpha = -1.0/(1.0 + sqrt(R/HPHR))

      ! recurring factor
      n_fac2 = n_fac/HPHR

      ! compute regression strenth and store in K
      call DGEMV('N', blsize_l1d, nmembers, n_fac2, devmat, blsize_l1d, devvec, 1, 0.0, K, 1)

      ! Update deviations using alpha K
      call DGER(blsize_l1d, nmembers, alpha, K, 1, devvec, 1, devmat, blsize_l1d)

      ! update the mean x = x  + res * K  or Hx = Hx + HK
      call daxpy(blsize_l1d, res, K, 1, meanvec, 1)

   end subroutine enksrf_kernel

   subroutine enksrf(nobs, nmembers, nparams, obs, mrej, rej_thr, hx, hxp, xp, r, x, rejected, assimilate)
      !
      ! This subroutine implements the Ensemble Kalman Square Root Filter algorithm as described in
      !    "Ensemble Data Assimilation without Perturbed Observations"
      ! by Whitaker and Hamill 2002. DOI: 10.1175/1520-0493(2002)130<1913:EDAWPO>2.0.CO;2
      ! and
      !    "An ensemble data assimilation system to estimate CO2 surface
      !    "fluxes from atmospheric trace gas observations"
      ! by Peters et al. 2005. DOI: 10.1029/2005JD006157
      !
      ! It achieves parallelziation trough cache blocking using the global parameter blsize_l1d and nrblocksl2 Tuning is required
      ! and depends on the amount of ensemble members

      ! Input variables, not changed during routine
      integer*8, intent(in)                                            :: nobs                                     ! the amount of observations
      integer*8, intent(in)                                            :: nmembers                                 ! the amount of ensemble members used to represent the covariance
      integer*8, intent(in)                                            :: nparams                                  ! number of statevector elements to be estimated
      logical, dimension(nobs), intent(in)                             :: mrej                                     ! may reject wether an observations may be rejected
      logical, dimension(nobs), intent(in)                             :: assimilate                               ! wether we want to assimilate this observation
      real*8, dimension(nobs), intent(in)                              :: obs                                      ! observation values
      real*8, dimension(nobs), intent(in)                              :: R                                        ! observation error
      real*8, dimension(nobs), intent(in)                              :: rej_thr                                  ! rejection threshold for observation that may be rejected
      ! local variables

      ! Output (updated in the loop)
      real*8, dimension(nobs, nmembers), intent(inout)                 :: HXp                                      ! observed deviations
      real*8, dimension(nobs), intent(inout)                           :: Hx                                       ! simulationed mean value of the observations
      real*8, dimension(nparams, nmembers), intent(inout)              :: Xp                                       ! statevecter deviations
      real*8, dimension(nparams), intent(inout)                        :: x                                        ! statevector
      integer*8, dimension(nobs), intent(inout)                        :: rejected                                 ! returns True if observation is rejected

      ! Local/temprorary
      integer*8                                                        :: i, j, jstart, jblock, jj, ii, jjstart, jjj, n !
      real*8                                                           :: n_fac                          ! for the repeated factor of  1/(N-1), n=nmembers
      real*8                                                           :: res                                ! residual, alpha factor (see Whitaker and Hamill 2002)
      real*8, dimension(nmembers)                                      :: HXp_n                                    ! temporary row of deviations
      integer*8, dimension(nobs)                                       :: rejectedp                                ! returns True if observation is rejected

      real*8, dimension(blsize_l1d, nmembers*(nobs/blsize_l1d + 1))    :: HXblock ! Transposed copyin of HX for contiguous access and first touch
      real*8, dimension(blsize_l1d, nmembers*(nparams/blsize_l1d + 1)) :: Xblock  ! Transposed copyin of X for contiguous access and first touch
      real*8, dimension(nmembers, nobs)                                :: HXblock_regres !Local variable recording of the state of HX_prime at loop index of assimilation so obs/state regression can be seperated
      real*8, dimension(nobs)                                          :: res_regres
      real*8, dimension(nmembers, blsize_l1d*nrblocksl2)               :: HXblock_regres_p ! threadprivate Local variable recording of the state of HX_prime at loop index of assimilation so obs/state regression can be seperated
      real*8, dimension(blsize_l2d)                                    :: res_regres_p ! threadprivate Records the residual at time of assimilation
      real*8, dimension(blsize_l1d, nmembers*nrblocksl2)               :: sblock ! temp cache sized block (threadprivate)
      real*8, dimension(blsize_l1d*nrblocksl2)                         :: svec, Rvec ! temp cache sized vectors (threadprivate)
      real*8, dimension((int(nobs/blsize_l1d) + 1)*blsize_l1d)         :: Hxpad! Padded vector for regularized kernel computations
      real*8, dimension((int(nparams/blsize_l1d) + 1)*blsize_l1d)      :: xpad ! Padded vector for regularized kernel computations
      real*8                                                           :: tstart ! Timer

      !! start parallel region
      !$omp parallel default(shared) private(res,HXp_n,jstart,j,n_fac,svec,sblock,ii,jjj,jj,jjstart,n,HXblock_regres_p,res_regres_p,jblock,rejectedp) firstprivate(mrej,assimilate,R,obs,rej_thr)
      !! copy in and reshape HXprime and Hx
      !! in the future this could potentially become a simple transpose, which makes the code cleaner and loop
      !! access much easier.
      !! in the enksrf_kernel then the dgemv can take in the matrix with the flag 'T' set. 
      !$omp do schedule(static)
      do jj = 1, nobs, blsize_l2d
      do ii = 1, nmembers
      do jjj = jj, blsize_l2d + jj - 1, blsize_l1d
         jjstart = 1 + ((jjj - 1)/(blsize_l1d))*nmembers
         HXblock(1:min(jjj + blsize_l1d - 1, nobs) - jjj + 1, jjstart + ii - 1) = HXp(jjj:min(jjj + blsize_l1d - 1, nobs), ii)
      end do
      end do
      end do
      !$omp end do

      !$omp do schedule(static)
      do j = 1, nobs, blsize_l1d
         Hxpad(j:min(j + blsize_l1d - 1, nobs)) = Hx(j:min(j + blsize_l1d - 1, nobs))
      end do
      !$omp end do

      ! init blocks and constants
      sblock = 0.0
      svec = 0.0

      call openblas_set_num_threads(1)

      n_fac = 1.0/(float(nmembers) - 1.0)

      !$omp master
      tstart = omp_get_wtime()
      !$omp end master
      !! start processing observations one at a time
      !! in blsize_l2d blocks
      do i = 1, nobs, blsize_l2d
         !! copy into sblock (each thread)
         jjstart = 1 + ((i - 1)/blsize_l2d)*nmembers*nrblocksl2
         sblock(:, 1:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1)) + 1 - jjstart) = HXblock(:, jjstart:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1)))
         svec(1:min(i + blsize_l2d - 1, nobs) - i + 1) = Hxpad(i:min(i + blsize_l2d - 1, nobs))

         !! process observations in this block
         !! storing the intermedite iith ensemble and residual 
         !! done by all cores
         do ii = i, min(i + blsize_l2d - 1, nobs)
            if (assimilate(ii) .eqv. .False.) then
               cycle
            end if
            ! compute difference between forecast and observed:
            res_regres_p(ii - i + 1) = obs(ii) - svec(ii - i + 1)

            ! We should be checking if we reject observations here by comparing to res(idual)
            if (mrej(ii) .eqv. .True.) then
               if (abs(res_regres_p(ii - i + 1)) > rej_thr(ii)*sqrt(R(ii))) then
                  rejectedp(ii) = 1
                  cycle
               end if
            end if
            !! Local copy of the deviations row HXp(i,:) will be updated later
            jblock = 1 + ((ii - i + 1 - 1)/(blsize_l1d))*nmembers
            HXblock_regres_p(:, ii - i + 1) = sblock(MOD(ii - 1, blsize_l1d) + 1, jblock:jblock + nmembers - 1)

            !! to copyout
            !! we don't use this so often, so it could be removed altogether
            !! also, the computation of HPRI/alhpa,n_fac2 can be moved inside the kernel to save some code lines
            !!$omp master 
            !HPHR(ii) = HPHRi
            !!$omp end master 

            ! Store local regression coefficient and residual to update HX_prime and later X_prime
            ! HXblock_regres_p(:, ii - i + 1) = HXp_n
            

            ! Update current block to proceed 
            do j = 1, nrblocksl2
               jstart = 1 + (j - 1)*nmembers
               call enksrf_kernel(nmembers, blsize_l1d, svec((j - 1)*blsize_l1d + 1:j*blsize_l1d),HXblock_regres_p(:, ii - i + 1), sblock(:, jstart:jstart + nmembers - 1),res_regres_p(ii - i + 1), n_fac, R(ii))

            end do
         end do

         !! potentially not required if HXblock_regres and res_regres become threadprivate
         !! store in the global variable
         !$omp master 
         HXblock_regres(:, i:min(blsize_l2d + i, nobs)) = HXblock_regres_p(:, 1:min(blsize_l2d + i, nobs) - i + 1)
         res_regres(i:min(blsize_l2d + i, nobs)) = res_regres_p(1:min(blsize_l2d + i, nobs) - i + 1)
         !$omp end master 
         !$omp barrier

         ! update all blocks in parallel
         ! per block do all obs in above block
         !$omp do schedule(static)
         do jj = 1, nobs, blsize_l2d
            ! copy in block
            jjstart = 1 + ((jj - 1)/blsize_l2d)*nmembers*nrblocksl2
            sblock(:, 1:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1)) + 1 - jjstart) = HXblock(:, jjstart:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1)))
            svec(1:min(jj + blsize_l2d - 1, nobs) - jj + 1) = Hxpad(jj:min(jj + blsize_l2d - 1, nobs))

            do j = 1, nrblocksl2
               jstart = 1 + (j - 1)*nmembers
               do ii = i, min(i + blsize_l2d - 1, nobs)
                  if (assimilate(ii) .eqv. .False.) then
                     cycle
                  end if
                  if (rejectedp(ii) == 1) then
                     cycle
                  end if

                  call enksrf_kernel(nmembers, blsize_l1d, svec((j - 1)*blsize_l1d + 1:j*blsize_l1d),HXblock_regres_p(:, ii - i + 1), sblock(:, jstart:jstart + nmembers - 1),res_regres_p(ii - i + 1), n_fac,R(ii))

               end do
            end do

            ! store block results until fetched later
            HXblock(:, jjstart:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1))) = sblock(:, 1:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nobs/blsize_l1d + 1)) + 1 - jjstart)
            Hxpad(jj:min(jj + blsize_l2d - 1, nobs)) = svec(1:min(jj + blsize_l2d - 1, nobs) - jj + 1)

         end do
         !$omp end do
         end do 
      !! end of  loop over observations
      !! all states of HX_prime and residuals have been deterimend and can opperate on Xprime/x
      !! first we can copyout HXp and Hx as we are done with those
      !$omp do schedule(static)
      do jj = 1, nobs, blsize_l2d
      do ii = 1, nmembers
      do jjj = jj, blsize_l2d + jj - 1, blsize_l1d
         jjstart = 1 + ((jjj - 1)/(blsize_l1d))*nmembers
         HXp(jjj:min(jjj + blsize_l1d - 1, nobs), ii) = HXblock(1:min(jjj + blsize_l1d - 1, nobs) - jjj + 1, jjstart + ii - 1)
      end do
      end do
      end do
      !$omp end do
      !$omp do schedule(static)
      do j = 1, nobs, blsize_l1d
         Hx(j:min(j + blsize_l1d - 1, nobs)) = Hxpad(j:min(j + blsize_l1d - 1, nobs))
      end do
      !$omp end do

      !$omp master
      rejected(:)=rejectedp(:)
      write (*, *) "after obs-obs loop computation ", (omp_get_wtime() - tstart)
      !$omp end master

      !! copy in now Xprime and x (also "first touch for NUMA")
      !$omp do schedule(static)
      do jj = 1, nparams, blsize_l2d
         do ii = 1, nmembers
            do jjj = jj, blsize_l2d + jj - 1, blsize_l1d
               jjstart = 1 + ((jjj - 1)/(blsize_l1d))*nmembers
               Xblock(1:min(jjj + blsize_l1d - 1, nparams) - jjj + 1, jjstart + ii - 1) = Xp(jjj:min(jjj + blsize_l1d - 1, nparams), ii)
            end do
         end do
      end do
      !$omp end do
      !$omp do schedule(static)
      do j = 1, nparams, blsize_l1d
         xpad(j:min(j + blsize_l1d - 1, nparams)) = x(j:min(j + blsize_l1d - 1, nparams))
      end do
      !$omp end do

      !! in parallel fetch block of x/HX_prime
      !! apply regression to block and be done with them
      !$omp do schedule(static)
      do jj = 1, nparams, blsize_l2d
         !! copy block Xprime and x
         jjstart = 1 + ((jj - 1)/blsize_l2d)*nmembers*nrblocksl2
         sblock(:, 1:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nparams/blsize_l1d + 1)) + 1 - jjstart) = Xblock(:, jjstart:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nparams/blsize_l1d + 1)))
         svec(1:min(jj + blsize_l2d - 1, nparams) - jj + 1) = xpad(jj:min(jj + blsize_l2d - 1, nparams))

         !! apply regression loop
         do n = 1, nobs
            if (assimilate(n) .eqv. .False.) then
               cycle
            end if
            if (rejectedp(n) == 1) then
               cycle
            end if

            !!  Apply
            do j = 1, nrblocksl2
               jstart = 1 + (j - 1)*nmembers
               call enksrf_kernel(nmembers, blsize_l1d, svec((j - 1)*blsize_l1d + 1:j*blsize_l1d), HXblock_regres(:, n), sblock(:, jstart:jstart + nmembers - 1), res_regres(n), n_fac, R(n))
            end do

            ! end loop over observations
         end do

         !! copy back into block
         Xblock(:, jjstart:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nparams/blsize_l1d + 1))) = sblock(:, 1:min(jjstart + nrblocksl2*nmembers - 1, nmembers*(nparams/blsize_l1d + 1)) + 1 - jjstart)
         xpad(jj:min(jj + blsize_l2d - 1, nparams)) = svec(1:min(jj + blsize_l2d - 1, nparams) - jj + 1)
      end do
      !$omp end do

      !$omp master
      write (*, *) "Time after regression computation ", (omp_get_wtime() - tstart)
      !$omp end master

      !! copy out results for x and Xprime
      !$omp do schedule(static)
      do j = 1, nparams, blsize_l1d
         x(j:min(j + blsize_l1d - 1, nparams)) = xpad(j:min(j + blsize_l1d - 1, nparams))
      end do
      !$omp end do
      !$omp do schedule(static)
      do jj = 1, nparams, blsize_l2d
      do ii = 1, nmembers
      do jjj = jj, blsize_l2d + jj - 1, blsize_l1d
         jjstart = 1 + ((jjj - 1)/(blsize_l1d))*nmembers
         Xp(jjj:min(jjj + blsize_l1d - 1, nparams), ii) = Xblock(1:min(jjj + blsize_l1d - 1, nparams) - jjj + 1, jjstart + ii - 1)
      end do
      end do
      end do
      !$omp end do

      !! end parallel region
      !$omp end parallel

   end subroutine enksrf

end module enkf_core
