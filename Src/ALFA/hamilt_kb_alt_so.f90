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

!>  Sets up the full hamiltonian matrix including spin-orbit
!>  for a given k-point and for a Kleinman-Bylander type pseudopotential
!>  Alternative method using proj_nl_kb
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         22 March 2024,
!>  \copyright    GNU Public License v2

subroutine hamilt_kb_alt_so(rkpt, mtxd, isort, qmod, ekpg,               &
    hamk,                                                                &
    ng, kgv, phase, conj, inds, kmax, indv,                              &
    veff, nqnl, delqnl, vkb, nkb,                                        &
    ntype, natom, rat, adot,                                             &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim, mxdsml)

! Adapted from hamilt_kb_alt.  22 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer , intent(in)               ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdsml                          !<  array dimension of hamiltonian

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  j-th component in lattice coordinates of the k-point
  integer, intent(in)                ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(in)                ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian
  real(REAL64), intent(in)           ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(in)           ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  real part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  real part of the ionic potential (hartree) for the prototype g-vector in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms

  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

! output

  complex(REAL64), intent(out)       ::  hamk(mxdsml,mxdsml)             !<  hamiltonian

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)

! counters

  integer   ::  i, j


! initialization

  hamk(:,:) = C_ZERO

! local potential without spin

  call hamilt_local(mtxd, isort, qmod, hamk,                             &
      ng, kgv, phase, conj, inds, kmax, indv,                            &
      veff,                                                              &
      mxdgve, mxdnst, mxdcub, mxddim, mxdsml)

! folds into large matrix

  do i = mtxd,1,-1
  do j = mtxd,1,-1
    hamk(2*i  ,2*j  ) = hamk(i,j)
    hamk(2*i-1,2*j-1) = hamk(i,j)
    hamk(2*i  ,2*j-1) = C_ZERO
    hamk(2*i-1,2*j  ) = C_ZERO
  enddo
  enddo

! kinetic energy

  do i = 1,mtxd
    hamk(2*i  ,2*i  ) = hamk(2*i  ,2*i  ) + cmplx(ekpg(i),ZERO,REAL64)
    hamk(2*i-1,2*i-1) = hamk(2*i-1,2*i-1) + cmplx(ekpg(i),ZERO,REAL64)
  enddo

! non-local part

  call hamilt_kb_non_local_so(rkpt, mtxd, isort, hamk,                   &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdsml)

  return

end subroutine hamilt_kb_alt_so
