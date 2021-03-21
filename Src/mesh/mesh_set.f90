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

!>  Distributes a quantity on the representative G-vector
!>  in an FFT mesh

  subroutine mesh_set(ipr, purpose, adot, den, rhomsh, ncheck,           &
    ng, kgv, phase, conj, inds, kmax,                                    &
    mxdgve, mxdnst, mxdscr)

! Written September 30, 2015 from setinmesh(CLR)
! Modified 12 September 2019, documentation test of mxdfft.  JLM
! Modified ipr, icheck, 13 February 2021. JLM                      WARNING   changed API
! copyright INESC-MN/Jose Luis Martins/Carlos Loia Reis

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdscr                          !<  array dimension for rhomsh

  integer, intent(in)                ::  ipr                             !<  contrlos printing
  character(len=*),  intent(in)      ::  purpose                         !<  characterization of den
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  ng                              !<  size of g-space
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  G-vectors in reciprocal lattice coordinates
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|

  complex(REAL64), intent(in)        ::  den(mxdnst)                     !<  density or other quantity in prototype G-vector

! output

  real(REAL64), intent(out)          ::  rhomsh(mxdscr)                  !<  density or other quantity on regular mesh in real space
  integer, intent(out)               ::  ncheck(4)                       !<  n1,n2,n3,id for consistency checking
! local allocatable arrays

  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:)

! local variables

  integer         ::  mxdfft, mxdwrk
  real(REAL64)    ::  vcell, bdot(3,3)
  integer         ::  nsfft(3), n1, n2, n3, id
  integer         ::  nn1, nn2, nn3
  real(REAL64)    ::  sum_rho, dmax, dmin, cmax, abschd
  logical         ::  lwrap
  integer         ::  ierr, iadd

! counters

  integer         ::  k1, k2, k3

! parameters

  real(REAL64), parameter :: EPS = 1.0E-09_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64


  call adot_to_bdot(adot, vcell, bdot)

  if(ipr > 2) write(6,*) ' mesh_set  ', purpose

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  if(mxdfft > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in mesh_set.  mxdfft = ",i8,                   &
           " is greater than mxdscr = ",i8)') mxdfft, mxdscr
    write(6,*) purpose

    stop

  endif

  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
  id = nsfft(1) + 1

  ncheck(1) = n1
  ncheck(2) = n2
  ncheck(3) = n3
  ncheck(4) = id

  nn1 = (n1-1) / 2
  nn2 = (n2-1) / 2
  nn3 = (n3-1) / 2

  if(ipr > 2) then
    write(6,*)
    write(6,'("  mesh_set  n = ",3i5)') n1,n2,n3
    write(6,*)
  endif

! initialize charge density array and enter symmetrized

  lwrap = .TRUE.
  if(kmax(1) < nn1 .and. kmax(2) < nn2 .and.                             &
      kmax(3) < nn3) then
!   this should be the normal case
    lwrap = .FALSE.
  endif

  call mesh_unfold(den, chd, id, n1,n2,n3, lwrap,                        &
      ng, kgv, phase, conj, inds,                                        &
      mxdgve, mxdnst, mxdfft)

! fourier transform to real space

  call cfft_c16(chd, id, n1,n2,n3, -1, wrkfft, mxdwrk)

! checks the correctness of the charge density
! and converts to electrons(whatever)/bohr^3

  sum_rho = ZERO
  dmax = real(chd(1))
  dmin = dmax
  cmax = ZERO
  ierr = 0
  do k3 = 1,n3
  do k2 = 1,n2
  do k1 = 1,n1
    iadd = ((k3-1)*n2 + (k2-1))*id + (k1-1) + 1
    if (real(chd(iadd)) > dmax) dmax = real(chd(iadd))
    if (real(chd(iadd)) < dmin) dmin = real(chd(iadd))
    abschd = abs(aimag(chd(iadd)))
    if (abschd > cmax) cmax = abschd
    if (abschd > eps) ierr = ierr+1
    rhomsh(iadd) = real(chd(iadd),REAL64) / vcell
    sum_rho = sum_rho + rhomsh(iadd)
  enddo
  enddo
  enddo
  sum_rho = sum_rho*vcell/(n1*n2*n3)

  if(ipr > 1) then
    write(6,*) purpose,' sum_rho =', sum_rho
    write(6,'("  max and min values:",3f14.6)') dmax,dmin,cmax
    write(6,*)
  endif

  deallocate(chd)
  deallocate(wrkfft)

  return
end subroutine mesh_set
