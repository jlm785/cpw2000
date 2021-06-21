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

!>  Calculates the product of the local potential times
!>  neig wavevectors. The local potential is dealt with
!>  fast fourier transforms. complex version

  subroutine hk_psi_loc_c16(mtxd, neig, psi, hpsi, ladd,                 &
    ng, kgv,                                                             &
    isort, vscr, kmscr,                                                  &
    mxddim, mxdbnd, mxdgve, mxdscr)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified 8 November 2015. Compatibility new libpw. JLM
! Modified, documentation, 26 January 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  logical, intent(in)                ::  ladd                            !<  true: adds to existing hpsi, false: input hpsi is zeroed

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  complex(REAL64), intent(in)        ::  psi(mxddim,mxdbnd)              !<  wavevector

! input and output

  complex(REAL64), intent(inout)     ::  hpsi(mxddim,mxdbnd)             !<  |hpsi> =  V_loc |psi>

! local allocatable arrays

  integer,allocatable                ::  ipoint(:)
  complex(REAL64),allocatable        ::  chd(:)
  real(REAL64),allocatable           ::  wrkfft(:)

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  integer                ::  kd1, kd2, kd3
  integer                ::  n1, n2, n3, id
  integer                ::  ntot,it
  integer                ::  nsfft(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer   ::   i, n, k1, k2, k3


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_loc_c16.  mxdwrk = ",i8,             &
           &  " is greater than mxdscr = ",2i8)') mxdwrk, mxdscr, ng

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
    write(6,*)  '   WARNING   in hk_psi_loc_c16:   mxdwrk may be',       &
      '  incorrectly calculated '
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

    write(6,'("   STOPPED in hk_psi_loc_c16:    dimension of k ",        &
          "index= ",3i5," exceeds ",3i5)') kd1,kd2,kd3,                  &
                     (n1-1)/2,(n2-1)/2,(n3-1)/2

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

    if(ladd) then

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        hpsi(i,n) = hpsi(i,n) + chd(ipoint(i))
      enddo
!$omp end parallel do

    else

!$omp parallel do default(shared) private(i)
      do i=1,mtxd
        hpsi(i,n) = chd(ipoint(i))
      enddo
!$omp end parallel do

    endif

!   end of loop over eigenvectors

  enddo



  deallocate(chd)
  deallocate(wrkfft)
  deallocate(ipoint)

  return
end subroutine hk_psi_loc_c16
