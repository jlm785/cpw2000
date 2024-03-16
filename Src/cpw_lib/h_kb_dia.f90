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

!>  Calculates the hamiltonian for one k-point
!>  and diagonalizes it
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         18 october 1993. 30 November 2023.
!>  \copyright    GNU Public License v2

subroutine h_kb_dia(emax, rkpt, neig, flgpsd,                            &
    ipr, ifail, icmax, iguess, epspsi,                                   &
    ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                      &
    sfact, veff, icmplx,                                                 &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    mtxd, hdiag, isort, qmod, ekpg, lkpg,                                &
    psi, hpsi, ei,                                                       &
    vscr, kmscr,                                                         &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim, mxdbnd, mxdscr)

! version 4.0. 18 october 1993. jlm
! modified 25 march 1999, 19 april 1999.
! modified (chd) 24 july 2002.
! modified 19 november 2013. veffr(1). jlm
! modified 19 december 2013. kmscr. jlm
! modified December 19, 2013 (f90). jlm
! modified, dimensions vkb, March 31, 2014. jlm
! modified, f90, complex, September-October 2015. JLM
! Modified kmscr, 28 October 2015. JLM
! Modified lkpg, October 2018. JLM
! Modified documentation, August 2019. JLM
! Modified, dimension of ei fixed. icmax,ifail. 18 February 2020. JLM
! Modified, hamilt_pw, 2 June 2020. JLM
! Modified, qmod-->ekpg in ditsp_c16. 13 February 2021. JLM
! Modified, nanlspin, 30 November 2023. JLM
! Added the commented out alternative call to hamilt_kb_alt. 16 Mrch 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdscr                          !<  array dimension for vscr

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential

  integer, intent(in)                ::  ipr                             !<  print control. 0 no printing, 2 lots of stuff
  integer, intent(in)                ::  icmax                           !<  maximum value of outer iteration
  integer, intent(in)                ::  iguess                          !<  tells if guess eigenvectors are available
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  real part of the phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of kgv(i,n)
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)            !<  real part of the structure factor
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex
  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  real part of the ionic potential (hartree) for the prototype g-vector in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh and fft mesh size
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  logical, intent(in)                ::  lkpg                            !<  If true use the previous G-vectors (same mtxd and isort)

! input and output

  integer, intent(inout)             ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)
  complex(REAL64), intent(inout)     ::  psi(mxddim,mxdbnd)              !<  component j of eigenvector i (guess on input)

  integer, intent(inout)             ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(inout)             ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

! output

  integer, intent(out)               ::  ifail                           !<  if ifail=0 ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  real(REAL64), intent(out)          ::  ei(mxdbnd)                      !<  eigenvalues (Hartree)
  real(REAL64), intent(out)          ::  hdiag(mxddim)                   !<  hamiltonian diagonal
  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  component j of eigenvector i (guess on input)


! local allocatable arrays

  real(REAL64), allocatable           ::  xnlkb(:)
  complex(REAL64), allocatable        ::  anlga(:,:)                     !  KB projectors
  complex(REAL64), allocatable        ::  hamsm(:,:)                     !  small hamiltonian

! local varaibles

  real(REAL64)    ::  veffr1
  integer         ::  mxdanl,mxdsml
  integer         ::  nanl, nanlso, nanlspin
  integer         ::  mtxds
  INTEGER         ::  NDUM

  logical         ::  lnewanl                                            !  indicates that anlga has been recalculated (not used in default implementation)


  if(flgpsd /= 'PSEUKB') then
    write(6,*)
    write(6,'("   stopped in hkbdia: wrong pseudo ",a6)') flgpsd

    stop

  endif

  veffr1 = real(veff(1),REAL64)

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
      mxdtyp)

  mxdanl = nanl

  allocate(xnlkb(mxdanl))
  allocate(anlga(mxddim,mxdanl))


  call hamilt_pw(emax,  rkpt,  lkpg,  veffr1,  nanl,                     &
      mtxd, isort, qmod, ekpg, hdiag,                                    &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  call size_mtxds(ipr, hdiag, neig, mtxd, mtxds)

  mxdsml = mtxds

  allocate(hamsm(mxdsml,mxdsml))

  call hamilt_kb(mtxds, hdiag, isort, qmod,                              &
      ng, kgv, phase, conj, inds, kmax, indv, ek,                        &
      sfact, veff, nqnl, delqnl, vkb, nkb,                               &
      ntype, adot, hamsm,                                                &
      mxdtyp, mxdgve, mxdnst, mxdcub, mxdlqp, mxdsml)

! the call to "old" hamilt_kb can be replaced to the call
! to the new alternative hamilt_kb_alt

!   call hamilt_kb_alt(rkpt, mtxds, isort, qmod, ekpg,                     &
!       hamsm,                                                             &
!       ng, kgv, phase, conj, inds, kmax, indv,                            &
!       veff, nqnl, delqnl, vkb, nkb,                                      &
!       ntype, natom, rat, adot,                                           &
!       mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim, mxdsml)

  NDUM = NEIG

  lnewanl = .TRUE.

  call ditsp_c16(ipr, ifail, icmax, iguess, epspsi, lnewanl,             &
      NDUM, mtxd, mtxds,                                                 &
      psi, hpsi, ei,                                                     &
      ekpg, isort, vscr, kmscr,                                          &
      ng, kgv,                                                           &
      anlga, xnlkb, nanl,                                                &
      hamsm, hdiag,                                                      &
      mxddim, mxdsml, mxdbnd, mxdgve, mxdscr, mxdanl)

  deallocate(hamsm)
  deallocate(xnlkb)
  deallocate(anlga)

  return
end subroutine h_kb_dia
