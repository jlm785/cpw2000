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

!>  Calculates the product of the hamiltonian times
!>  neig wavevectors. The non-local pseudopotential
!>  is separable. The local potential is dealt with
!>  fast Fourier transforms. Complex version

!>  This is the generic version. There are special versions for
!>  combinations of libraries and CPU/GPU

subroutine hk_psi_c16(mtxd, neig, psi, hpsi, lnewanl,                    &
    ng, kgv,                                                             &
    ekpg, isort, vscr, kmscr,                                            &
    anlga, xnlkb, nanl,                                                  &
    mxddim, mxdbnd, mxdanl, mxdgve, mxdscr)

! written february 2, 1990. jlm
! modified june 3,1999. jlm
! modified (chd) 25 july 2002. jlm
! modified openmp july 2013.   jlm
! modified complex February 6 2014. jlm
! Modified kmscr, 28 October 2015. JLM
! Modified, documentation, January 2020. JLM
! Replaced qmod with ekpg.  13 February 2021. JLM         WARNING  API modification
! copyright INESC-MN/Jose Luis Martins

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdanl                          !<  array dimension of number of projectors
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nanl                            !<  number of projectors
  complex(REAL64), intent(in)        ::  anlga(mxddim,mxdanl)            !<  Kleinman-Bylander projectors
  real(REAL64), intent(in)           ::  xnlkb(mxdanl)                   !<  Kleinman-Bylander normalization
  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

! input and output

  logical, intent(inout)             ::  lnewanl                         !<  indicates that anlga has been recalculated (not used in default implementation)

! output

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  V_NL |psi>

! local allocatable arrays

  integer,allocatable          ::  ipoint(:)
  complex(REAL64),allocatable  ::  chd(:)
  real(REAL64),allocatable     ::  wrkfft(:)

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  integer         ::  kd1, kd2, kd3
  integer         ::  n1, n2, n3, id
  integer         ::  ntot,it
  integer         ::  nsfft(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer   ::   i, n, k1, k2, k3


! find n for fast fourier transform
! ni is the number of points used in direction i.

  if(ng > mxdgve) stop

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_c16.  mxdwrk = ",i8,                 &
           " is greater than mxdscr = ",2i8)') mxdwrk, mxdscr, ng

    stop

  endif

  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)
!  id = nsfft(1) + 1
  id = kmscr(7)
  ntot = id * n2 * n3

  if(ntot > mxdfft) then
    mxdfft = ntot
    write(6,*)
    write(6,*)  '   WARNING   in hk_psi_c16:   mxdwrk may be',           &
      '  incorrectly calculated ', lnewanl
    write(6,*)
  endif

  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))

! fills array ipoint

  allocate(ipoint(mtxd))

  kd1 = 0
  kd2 = 0
  kd3 = 0
  do i=1,mtxd
    it = isort(i)
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

    write(6,'("   STOPPED in hk_psi_c16:      dimension of k ",          &
          "index= ",3i5," exceeds ",3i5)') kd1,kd2,kd3,                  &
                     (n1-1)/2,(n2-1)/2,(n3-1)/2
    write(6,*) ' You are probably using the dual approximation'
    write(6,*) ' with a k-point far away from the Brillouin zone'

    stop

  endif

! local potential and kinetic energy

  do n=1,neig

!   loop over eigenvectors

!$omp parallel do default(shared) private(i)
    do i=1,ntot
       chd(i) = cmplx(ZERO,ZERO,REAL64)
    enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
    do i=1,mtxd
      chd(ipoint(i)) = psi(i,n)
    enddo
!$omp end parallel do

!   fourier transform to real space

    call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

!   calculates the product of vscr and chd

!$omp parallel do default(shared) private(i)
    do i=1,ntot
      chd(i) = vscr(i)*chd(i)
    enddo
!$omp end parallel do

!   fourier transforms to g-space

    call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)

!   fills hpsi

!$omp parallel do default(shared) private(i)
    do i=1,mtxd
      hpsi(i,n) = chd(ipoint(i))
    enddo
!$omp end parallel do

!   kinetic energy

!$omp parallel do default(shared) private(i)
    do i=1,mtxd
      hpsi(i,n) = hpsi(i,n) + ekpg(i)*psi(i,n)
    enddo
!$omp end parallel do


!   end of loop over eigenvectors

  enddo

! non local potential

  call hk_psi_nl_c16(mtxd, neig, psi, hpsi, anlga, xnlkb, nanl, .TRUE.,  &
      mxddim,mxdbnd,mxdanl)


  deallocate(chd)
  deallocate(wrkfft)
  deallocate(ipoint)

  return
end subroutine hk_psi_c16
