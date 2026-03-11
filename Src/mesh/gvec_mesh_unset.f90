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

!>  Gathers from an FFT mesh a quantity on the representative G-vector
!>  Inverse of gvec_mesh_set
!>
!>  \author       José Luís Martins
!>  \version      5.13
!>  \date         11 March 2026.
!>  \copyright    GNU Public License v2

subroutine gvec_mesh_unset(ipr, purpose, adot, den,                      &
    rhomsh, id,n1,n2,n3, lvol,                                           &
    ng, kgv, phase, conj, ns, inds, kmax, mstar,                         &
    mxdgve, mxdnst, mxdscr)

! Written inverting mesh_set. 11 march 2026. JLM

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
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  mstar(mxdnst)                   !<  number of g-vectors in the j-th star

  real(REAL64), intent(in)           ::  rhomsh(mxdscr)                  !<  density or other quantity on regular mesh in real space

  integer, intent(in)                ::  id,n1,n2,n3                     !<  packing of rhomsh(id,n2,n3), id >= n1

  logical, intent(in)                ::  lvol                            !<  if true scales result by volume

! output

  complex(REAL64), intent(out)       ::  den(mxdnst)                     !<  density or other quantity in prototype G-vector

! local allocatable arrays

  real(REAL64), allocatable          ::  wrkfft(:)
  complex(REAL64), allocatable       ::  chd(:)
  complex(REAL64), allocatable       ::  deng(:)

! local variables

  integer         ::  mxdfft, mxdwrk
  real(REAL64)    ::  vcell, bdot(3,3)
  integer         ::  nsfft(3)
  integer         ::  nn1, nn2, nn3
  integer         ::  iadd

  real(REAL64)    ::  fac

! counters

  integer         ::  i, k1, k2, k3

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64


  if(lvol) then
    call adot_to_bdot(adot, vcell, bdot)
    fac = vcell
  else
    fac = UM
  endif

  if(ipr > 2) write(6,*) ' gvec_mesh_unset  ', purpose

! some paranoid checks
! find n for fast fourier transform
! ni is the number of points used in direction i.
! note that mxdfft mxdwrk may be dependent on fft package
! so it is safer to call size_fft again

  call size_fft(kmax, nsfft, mxdfft, mxdwrk)

  if(nsfft(1) /= n1 .or. nsfft(2) /= n2 .or. nsfft(3) /= n3              &
      .or. id < n1) then
    write(6,*)
    write(6,*)  "  STOPPED in gvec_mesh_unset applied to ", purpose
    write(6,'("  in v_hxc ",4i6," in gvec_mesh_set ",3i6)') n1,n2,n3,id, &
               (nsfft(i),i=1,3)
    write(6,*)

    stop

  endif

  if(mxdfft > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in gvec_mesh_unset.  mxdfft = ",i8,            &
          & " is greater than mxdscr = ",i8)') mxdfft, mxdscr
    write(6,*) purpose

    stop

  endif

  allocate(chd(mxdfft))
  allocate(wrkfft(mxdwrk))
  allocate(deng(mxdgve))

  nn1 = (n1-1) / 2
  nn2 = (n2-1) / 2
  nn3 = (n3-1) / 2

  if(ipr > 2) then
    write(6,*)
    write(6,'("  gvec_mesh_set  n = ",3i5)') n1,n2,n3
    write(6,*)
  endif

  do k3 = 1,n3
  do k2 = 1,n2
  do k1 = 1,n1
    iadd = ((k3-1)*n2 + (k2-1))*id + (k1-1) + 1
    chd(iadd) = rhomsh(iadd) * fac
  enddo
  enddo
  enddo

! fourier transform to reciprocal space

  call cfft_c16(chd, id, n1,n2,n3, 1, wrkfft, mxdwrk)

! initialize charge density array and enter symmetrized

  call gvec_mesh_fold(deng, chd, id, n1,n2,n3,                           &
      ng, kgv,                                                           &
      mxdgve, mxdfft)

  call gvec_star_of_g_fold(den, deng, .FALSE.,                           &
      ng, phase, conj, ns, inds, mstar,                                  &
      mxdgve, mxdnst)

  deallocate(chd)
  deallocate(wrkfft)

  return

end subroutine gvec_mesh_unset
