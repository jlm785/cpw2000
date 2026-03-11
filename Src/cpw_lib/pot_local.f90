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
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         20 February 2018. 11 March 2026.
!>  \copyright    GNU Public License v2

subroutine pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,        &
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

  integer, intent(in)                ::  idshift                         !<  shift of the fft mesh, used /= 0 only in highly banked memory.

! input and output

  integer, intent(inout)             ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

! output

  real(REAL64),  intent(out)         ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  real(REAL64),  intent(out)         ::  vmax, vmin                      !<  maximum and minimum values of vscr

! local variables

  integer        ::  mxdfft                                              !  array dimension for fft transform
  integer        ::  mxdwrk                                              !  array dimension for fft transform workspace

  integer        ::  id,n1,n2,n3
  integer        ::  nsfft(3)

  real(REAL64)   ::  adot(3,3)                                           !  unused in gvec_mesh_set
  integer        ::  kmax(3)                                             !  max value of |kgv(i,n)|



! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter  :: SMALL = 1.0E-9_REAL64

! counters

  integer    ::  i, j, k, ijk


! compatibility with new gvec_mesh_set, adot not used, kmax recalculated to avoid new API

  adot(:,:) = ZERO
  do i = 1,3
    adot(i,i) = UM
  enddo

  do i = 1,3
    kmax(i) = 0
    do j = 1,ng
      if(iabs(kgv(i,j)) > kmax(i)) kmax(i) = iabs(kgv(i,j))
    enddo
  enddo

! printout local potential

  if (ipr > 2) then
    write(6,*)
    write(6,*)  '  local potential in g-space'
    write(6,*)
    write(6,'(20f10.3)') (veff(i),i=1,ns)
    write(6,*)
  endif

! find n for fast fourier transform
! ni is the number of points used in direction i.

  call size_fft(kmscr, nsfft, mxdfft, mxdwrk)

  if(mxdfft > mxdscr) then
    write(6,*)
    write(6,'("   STOPPED in pot_local.  mxdfft = ",i8,                  &
          & " is greater than mxdscr = ",i8)') mxdwrk, mxdscr

    stop

  endif

  n1 = nsfft(1)
  n2 = nsfft(2)
  n3 = nsfft(3)
!  id = nsfft(1) + 1
  id = nsfft(1) + idshift

  kmscr(4) = n1
  kmscr(5) = n2
  kmscr(6) = n3
  kmscr(7) = id

  if (ipr /= 0) then
    write(6,*)
    write(6,'("  in fft for local potential n =",3i5)') n1,n2,n3
  endif

! initialize charge density array and enter symmetrized charge.

  call gvec_mesh_set(ipr, 'potential', adot, veff,                       &
    vscr, id,n1,n2,n3, .FALSE.,                                          &
    ng, kgv, phase, conj, inds, kmax,                                    &
    mxdgve, mxdnst, mxdfft)

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
