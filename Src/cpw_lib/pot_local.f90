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

!>  Computes the local potential on a grid
!>  using fast fourier transforms.
!>
!>  The description of the fft mesh (kmscr) is calculated
!>  outside, by size_kmscr, so that it can be shared with other
!>  quantities in real space.
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         20 February 2018. 24 September 2026.
!>  \copyright    GNU Public License v2

subroutine pot_local(ipr, vscr, vmax, vmin, veff, kmscr, kmax,           &
    ng, kgv, phase, conj, ns, inds,                                      &
    mxdscr,mxdgve, mxdnst)

! written august 6 1987. jlm
! modified august 31 1987. jlm
! modified 4.0. 17 october 93. jlm
! modified (conj) 19 march 99. jlm
! modified (chd) 26 july 2002. jlm
! modified (f90) 13 January 2014. jlm
! modified mesh_unfold. 11 October 2015. JLM
! Modified kmscr, 28 October 2015. JLM
! Modified, documentation, January 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, calls gvec_mesh_set. 11 March 2026. JLM
! Modified, kmscr(4:7) and idshift moved to size_kmscr,
!           kmax is input.  24 September 2026. JLM            WARNING NEW API

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars

  integer, intent(in)                ::  ipr                             !<  print option, ipr=0 no printing, ipr=1 some printing, ipr=2 a lot of printing
  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  ionic potential (local+Hartree+XC) for the prototype g-vector in star j
  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax

  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size (from size_kmscr)
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|

! output

  real(REAL64),  intent(out)         ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  real(REAL64),  intent(out)         ::  vmax, vmin                      !<  maximum and minimum values of vscr

! local variables

  integer        ::  id,n1,n2,n3

  real(REAL64)   ::  adot(3,3)                                           !  unused in gvec_mesh_set

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer    ::  i, j, k, ijk


! compatibility with gvec_mesh_set, adot not used

  adot(:,:) = ZERO
  do i = 1,3
    adot(i,i) = UM
  enddo

! printout local potential

  if (ipr > 2) then
    write(6,*)
    write(6,*)  '  local potential in g-space'
    write(6,*)
    write(6,'(20f10.3)') (veff(i),i=1,ns)
    write(6,*)
  endif

! fft mesh calculated by size_kmscr

  n1 = kmscr(4)
  n2 = kmscr(5)
  n3 = kmscr(6)
  id = kmscr(7)

  if (ipr /= 0) then
    write(6,*)
    write(6,'("  in fft for local potential n =",3i5)') n1,n2,n3
  endif

! initialize charge density array and enter symmetrized charge.
! gvec_mesh_set checks that the mesh fits in vscr

  call gvec_mesh_set(ipr, 'potential', adot, veff,                       &
    vscr, id,n1,n2,n3, .FALSE.,                                          &
    ng, kgv, phase, conj, inds, kmax,                                    &
    mxdgve, mxdnst, mxdscr)

  vmax = vscr(1)
  vmin = vmax
  do k = 1,n3
  do j = 1,n2
  do i = 1,n1
    ijk = ((k-1)*n2 + j-1)*id + i
    if (vscr(ijk) > vmax) vmax = vscr(ijk)
    if (vscr(ijk) < vmin) vmin = vscr(ijk)
  enddo
  enddo
  enddo

  return

end subroutine pot_local
