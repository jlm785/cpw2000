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
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         February 18, 2014. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine hk_psi_loc_c16(mtxd, neig, psi, hpsi, ladd, nspin, nsp,       &
    ng, kgv,                                                             &
    isort, vscr, kmscr,                                                  &
    mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

! Written February 18, 2014, from hk_psi_c16.   jlm
! See that file for historical record.
! Modified 8 November 2015. Compatibility new libpw. JLM
! Modified, documentation, 26 January 2020. JLM
! Added spin polarization option, 4 November 2025. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves (not counting spin)
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands (including spin)
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr
  integer, intent(in)                ::  mxdnsp                          !<  array dimension for number of spin components (1,2,4)

  integer, intent(in)                ::  nspin                           !<  spin components (1:no spin or 2:spin present)

  integer, intent(in)                ::  mtxd                            !<  wavefunction dimension (basis size, not counting spin)
  integer, intent(in)                ::  neig                            !<  number of wavefunctions (including spin)
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  logical, intent(in)                ::  ladd                            !<  true: adds to existing hpsi, false: input hpsi is zeroed

  real(REAL64), intent(in)           ::  vscr(mxdscr,mxdnsp)             !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh
  integer, intent(in)                ::  nsp                             !<  number of spin components ox xc-potential (1,2,4)

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  complex(REAL64), intent(in)        ::  psi(nspin*mxddim,mxdbnd)        !<  wavevector

! input and output

  complex(REAL64), intent(inout)     ::  hpsi(nspin*mxddim,mxdbnd)       !<  |hpsi> =  V_loc |psi>

! local allocatable arrays

  integer,allocatable                ::  ipoint(:)
  complex(REAL64),allocatable        ::  chd(:)
  complex(REAL64),allocatable        ::  chd_m(:)                        !  s = -1/2
  complex(REAL64),allocatable        ::  chd_p(:)                        !  s =  1/2
  real(REAL64),allocatable           ::  wrkfft(:)

! local variables

  integer                ::  mxdfft                                      !  array dimension for fft transform
  integer                ::  mxdwrk                                      !  array dimension for fft transform workspace

  complex(REAL64)        ::  tmp_m, tmp_p

  integer                ::  kd1, kd2, kd3
  integer                ::  n1, n2, n3, id
  integer                ::  ntot,it
  integer                ::  nsfft(3)

! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter :: C_I = cmplx(ZERO,UM,REAL64)

! counters

  integer   ::   i, n, k1, k2, k3


! paranoid check, stops complaint about unused ng...

  if(ng > mxdgve) stop

! checks input

  if(nspin /=1 .and. nspin /= 2) then
    write(6,*)
    write(6,*) '    STOPPED in hk_psi_loc_c16, nspin = ', nspin
    write(6,*)

    STOP

  endif

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdwrk > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_loc_c16.  mxdwrk = ",i8,             &
           &  " is greater than mxdscr = ",2i8)') mxdwrk, mxdscr, ng

    stop

  endif

  if(nsp /=1 .and. nsp /=2 .and. nsp /= 4) then
    write(6,*)
    write(6,'("   STOPPED in hk_psi_loc_c16.  nsp = ",i8,                &
           &  " is an incorrect value, should be 1,2,4")') nsp

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

  if(nspin == 1) then
    allocate(chd(mxdfft))
  else
    allocate(chd_m(mxdfft), chd_p(mxdfft))
  endif

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
    write(6,*) ' You are probably using the dual approximation'
    write(6,*) ' with a k-point far away from the Brillouin zone'

    stop

  endif

! local potential

  if(nspin == 1) then                                                    !  no spin

    do n=1,neig

!     loop over eigenvectors

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

!     fourier transform to real space

      call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

!     calculates the product of vscr and chd

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd(i) = vscr(i,1)*chd(i)
      enddo
!$omp end parallel do

!     fourier transforms to g-space

      call cfft_wf_c16(chd, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)

!     fills hpsi

      if(ladd) then

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(i,n) = hpsi(i,n) + chd(ipoint(i))
        enddo
!$omp end parallel do

      else

!$omp parallel do default(shared) private(i)
        do i = 1,mtxd
          hpsi(i,n) = chd(ipoint(i))
        enddo
!$omp end parallel do

      endif

    enddo                                                                !   end of loop over eigenvectors

  else



    do n = 1,neig

!     loop over eigenvectors

!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd_p(i) = C_ZERO
      enddo
!$omp end parallel do
!$omp parallel do default(shared) private(i)
      do i = 1,ntot
        chd_m(i) = C_ZERO
      enddo
!$omp end parallel do

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        chd_p(ipoint(i)) = psi(2*i-1,n)
      enddo
!$omp end parallel do
!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        chd_m(ipoint(i)) = psi(2*i  ,n)
      enddo
!$omp end parallel do

!     fourier transform to real space

      call cfft_wf_c16(chd_p, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)
      call cfft_wf_c16(chd_m, id, n1,n2,n3, kd1,kd2,kd3, -1, wrkfft, mxdwrk)

!     calculates the product of vscr_sp and chd

      if(nsp == 1) then

!$omp parallel do default(shared) private(i)
        do i = 1,ntot
          chd_p(i) = vscr(i,1)*chd_p(i)
        enddo
!$omp end parallel do
!$omp parallel do default(shared) private(i)
        do i = 1,ntot
          chd_m(i) = vscr(i,1)*chd_m(i)
        enddo
!$omp end parallel do

      elseif(nsp ==2) then

!$omp parallel do default(shared) private(i)
        do i = 1,ntot
          chd_p(i) = (vscr(i,1)+vscr(i,2))*chd_p(i)
        enddo
!$omp end parallel do
!$omp parallel do default(shared) private(i)
        do i = 1,ntot
          chd_m(i) = (vscr(i,1)-vscr(i,2))*chd_m(i)
        enddo
!$omp end parallel do

      elseif(nsp ==4) then

!$omp parallel do default(shared) private(i, tmp_m, tmp_p)
        do i = 1,ntot
          tmp_p = ( vscr(i,1) + vscr(i,2) )*chd_p(i) +                   &
                ( vscr(i,3) - C_I*vscr(i,4) )*chd_m(i)
          tmp_m = ( vscr(i,1) - vscr(i,2) )*chd_m(i) +                   &
                ( vscr(i,3) + C_I*vscr(i,4) )*chd_p(i)
          chd_p(i) = tmp_p
          chd_m(i) = tmp_m
        enddo
!$omp end parallel do

      endif

!     fourier transforms to g-space

      call cfft_wf_c16(chd_p, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)
      call cfft_wf_c16(chd_m, id, n1,n2,n3, kd1,kd2,kd3,  1, wrkfft, mxdwrk)

!     fills hpsi

!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        hpsi(2*i-1,n) = chd_p(ipoint(i))
      enddo
!$omp end parallel do
!$omp parallel do default(shared) private(i)
      do i = 1,mtxd
        hpsi(2*i  ,n) = chd_m(ipoint(i))
      enddo
!$omp end parallel do

    enddo

  endif

  deallocate(wrkfft)
  if(nspin == 1) then
    deallocate(chd)
  else
    deallocate(chd_m, chd_p)
  endif
  deallocate(ipoint)

  return

end subroutine hk_psi_loc_c16
