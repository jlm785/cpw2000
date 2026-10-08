!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  Calculates, for the generalized Kohn-Sham meta-GGA, the product
!>  of vtau and of the derivative with respect to k of the vtau term
!>  of the hamiltonian times neig wavevectors.
!>
!>    < k+G | H_tau | k+G' > = 1/2 q . bdot . q' vtau(G-G'),   q = k+G
!>
!>    < k+G | d H_tau / d k_l | k+G' > =
!>                    1/2 [ (bdot q)_l + (bdot q')_l ] vtau(G-G')
!>
!>  with k and G in reciprocal lattice coordinates, the same
!>  convention as berry_dhdk_psi.  The second derivative is
!>  bdot_lm vtau(G-G'), so vtaupsi is also returned.
!>
!>  Generic version (only the generic FFT is used).  Complex version.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         8 October 2026.
!>  \copyright    GNU Public License v2

subroutine dvtaudk_psi_c16(mtxd, neig, psi, vtaupsi, dvtaupsi,           &
    rkpt, adot, isort, kgv, vtaumsh, kmscr,                              &
    mxddim, mxdbnd, mxdgve, mxdscr)

! Written 8 October 2026 from hk_psi_mgga_c16. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vtaumsh

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevectors

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  real(REAL64), intent(in)           ::  vtaumsh(mxdscr)                 !<  d (rho eps_xc) / d tau in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  complex(REAL64), intent(out)       ::  vtaupsi(mxddim,mxdbnd)          !<  vtau |psi> (convolution in G-space)
  complex(REAL64), intent(out)       ::  dvtaupsi(mxddim,mxdbnd,3)       !<  (d H_tau / d k_l) |psi>, k in lattice coordinates

! local allocatable arrays

  integer,allocatable          ::  ipoint(:)
  real(REAL64),allocatable     ::  bq(:,:)                               !  bdot (k+G), contravariant components
  complex(REAL64),allocatable  ::  chd(:)
  real(REAL64),allocatable     ::  wrkfft(:)

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  integer         ::  kd1, kd2, kd3
  integer         ::  n1, n2, n3, id
  integer         ::  ntot, it
  integer         ::  nsfft(3)

  real(REAL64)    ::  vcell, bdot(3,3)
  real(REAL64)    ::  q(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer   ::   i, j, l, n, k1, k2, k3


  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in dvtaudk_psi_c16.  mxdwrk = ",i8,            &
           &  " is greater than mxdscr = ",i8)') mxdwrk, mxdscr

    stop

  endif

  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)
  id = kmscr(7)
  ntot = id * n2 * n3

  if(ntot > mxdfft) mxdfft = ntot

  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))

  call adot_to_bdot(adot, vcell, bdot)

! fills array ipoint and the contravariant components of k+G

  allocate(ipoint(mtxd))
  allocate(bq(3,mtxd))

  kd1 = 0
  kd2 = 0
  kd3 = 0
  do i = 1,mtxd
    it = isort(i)
    do j = 1,3
      q(j) = rkpt(j) + kgv(j,it)
    enddo
    do j = 1,3
      bq(j,i) = bdot(j,1)*q(1) + bdot(j,2)*q(2) + bdot(j,3)*q(3)
    enddo
    k1 = kgv(1,it)
    if (iabs(k1) > kd1) kd1 = iabs(k1)
    if (k1 < 0) k1 = n1 + k1
    k2 = kgv(2,it)
    if (iabs(k2) > kd2) kd2 = iabs(k2)
    if (k2 < 0) k2 = n2 + k2
    k3 = kgv(3,it)
    if (iabs(k3) > kd3) kd3 = iabs(k3)
    if (k3 < 0) k3 = n3 + k3
    ipoint(i) = (k3*n2 + k2)*id + k1 + 1
  enddo

! no wraparound for eigenvectors

  if(kd1 > (n1-1)/2 .or. kd2 > (n2-1)/2 .or. kd3 > (n3-1)/2) then

    write(6,'("   STOPPED in dvtaudk_psi_c16:   dimension of k ",        &
         &  "index= ",3i5," exceeds ",3i5)') kd1,kd2,kd3,                &
                     (n1-1)/2,(n2-1)/2,(n3-1)/2

    stop

  endif

  do n = 1,neig

!   vtau |psi>

!$omp parallel do default(shared) private(i)
    do i = 1,ntot
      chd(i) = C_ZERO
    enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
    do i = 1,mtxd
      chd(ipoint(i)) = psi(i,n)
    enddo
!$omp end parallel do

    call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

!$omp parallel do default(shared) private(i)
    do i = 1,ntot
      chd(i) = vtaumsh(i)*chd(i)
    enddo
!$omp end parallel do

    call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)

!$omp parallel do default(shared) private(i)
    do i = 1,mtxd
      vtaupsi(i,n) = chd(ipoint(i))
    enddo
!$omp end parallel do

!   1/2 [ (bdot q)_l vtau |psi> + vtau (bdot q)_l |psi> ]

    do l = 1,3

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd(i) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        chd(ipoint(i)) = bq(l,i)*psi(i,n)
      enddo
!$omp end parallel do

      call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd(i) = vtaumsh(i)*chd(i)
      enddo
!$omp end parallel do

      call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        dvtaupsi(i,n,l) = (bq(l,i)*vtaupsi(i,n) + chd(ipoint(i))) / 2
      enddo
!$omp end parallel do

    enddo

  enddo

  deallocate(chd)
  deallocate(wrkfft)
  deallocate(ipoint)
  deallocate(bq)

  return

end subroutine dvtaudk_psi_c16
