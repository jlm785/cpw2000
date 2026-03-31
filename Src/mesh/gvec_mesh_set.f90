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
!>
!>  \author       Carlos Loia Reis, José Luís Martins
!>  \version      5.13
!>  \date         September 30 2015, 10 March 2026.
!>  \copyright    GNU Public License v2

subroutine gvec_mesh_set(ipr, purpose, adot, den,                        &
    rhomsh, id,n1,n2,n3, lvol,                                           &
    ng, kgv, phase, conj, inds, kmax,                                    &
    mxdgve, mxdnst, mxdscr)

! Written September 30, 2015 from setinmesh(CLR)
! Modified 12 September 2019, documentation test of mxdfft.  JLM
! Modified ipr, icheck, 13 February 2021. JLM
! Name. n1,n2,n3, lvol. 11 March 2026. JLM                       WARNING NEW API

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

  integer, intent(in)                ::  id,n1,n2,n3                     !<  packing of rhomsh(id,n2,n3), id >= n1

  logical, intent(in)                ::  lvol                            !<  if true scales result by volume

! output

  real(REAL64), intent(out)          ::  rhomsh(mxdscr)                  !<  density or other quantity on regular mesh in real space

! local allocatable arrays

  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:)
  complex(REAL64), allocatable       ::  deng(:)

! local variables

  integer         ::  mxdfft, mxdwrk
  real(REAL64)    ::  vcell, bdot(3,3)
  integer         ::  nsfft(3)
  integer         ::  knn(3)
  real(REAL64)    ::  sum_rho, dmax, dmin, cmax, abschd
  logical         ::  lwrap
  integer         ::  ierr, iadd

  real(REAL64)    ::  fac

! counters

  integer         ::  i, k1, k2, k3

! parameters

  real(REAL64), parameter :: EPS = 1.0E-09_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64



  if(lvol) then
    call adot_to_bdot(adot, vcell, bdot)
    fac = UM / vcell
  else
    fac = UM
  endif

  if(ipr > 2) write(6,*) '  gvec_mesh_set  ', purpose

  knn(1) = (n1-1) / 2
  knn(2) = (n2-1) / 2
  knn(3) = (n3-1) / 2

! some paranoid checks
! find n for fast fourier transform
! ni is the number of points used in direction i.
! note that mxdfft mxdwrk may be dependent on fft package
! so it is safer to call size_fft again

  call size_fft(knn, nsfft, mxdfft, mxdwrk)

  if(nsfft(1) /= n1 .or. nsfft(2) /= n2 .or. nsfft(3) /= n3              &
      .or. id < n1) then
    write(6,*)
    write(6,*)  "  STOPPED in gvec_mesh_set applied to ", purpose
    write(6,'("  mesh in calling sub: ",4i6," in gvec_mesh_set: ",3i6)') &
               n1,n2,n3,id, (nsfft(i),i=1,3)
    write(6,*)

    stop

  endif

  if(mxdfft > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in gvec_mesh_set.  mxdfft = ",i8,              &
          & " is greater than mxdscr = ",i8)') mxdfft, mxdscr
    write(6,*) purpose

    stop

  endif

  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))
  allocate(deng(mxdgve))

  if(ipr > 2) then
    write(6,*)
    write(6,'("  gvec_mesh_set  n = ",3i5)') n1,n2,n3
    write(6,*)
  endif

! initialize charge density array and enter symmetrized

  lwrap = .TRUE.
  if(kmax(1) < knn(1) .and. kmax(2) < knn(2) .and.                       &
      kmax(3) < knn(3)) then
!   this should be the normal case
    lwrap = .FALSE.
  endif

  call gvec_star_of_g_unfold(deng, den, .FALSE.,                         &
      ng, phase, conj, inds,                                             &
      mxdgve, mxdnst)

  call gvec_mesh_unfold(deng, chd, id, n1,n2,n3, lwrap,                  &
      ng, kgv,                                                           &
      mxdgve, mxdfft)

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
    if (abschd > EPS) ierr = ierr+1
    rhomsh(iadd) = real(chd(iadd),REAL64) * fac
    sum_rho = sum_rho + rhomsh(iadd)
  enddo
  enddo
  enddo
  sum_rho = sum_rho / (n1*n2*n3 * fac)

  if (ierr /= 0) then
    write(6,*)
    write(6,*) '    WARNING in gvec_mesh_set for ', purpose
    write(6,'("   complex function in", i12," points, cmax = ",e12.4)')  &
                 ierr,cmax
  endif

  if (cmax > 1000.0*EPS) then
    write(6,*)
    write(6,*) '   STOPPED in gvec_mesh_set for ', purpose
    write(6,'("   cmax = ",e12.4)') cmax

    stop

  endif

  if(ipr > 1) then
    write(6,*) purpose,' sum_rho =', sum_rho
    write(6,'("  max and min values:",3f14.6)') dmax,dmin,cmax
    write(6,*)
  endif

  deallocate(chd)
  deallocate(wrkfft)

  return

end subroutine gvec_mesh_set
