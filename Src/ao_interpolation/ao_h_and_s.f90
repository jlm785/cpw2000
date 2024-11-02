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

!>  Calculates the hamiltonian and overlap for one k-point
!>  with non-orthogonal atomic orbitals.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         before 2015. 1 November 2024.
!>  \copyright    GNU Public License v2

subroutine ao_h_and_s(emax, rkpt, nbaslcao, flgpsd,                      &
    ng, kgv,                                                             &
    veff, icmplx,                                                        &
    nqnl, delqnl, vkb, nkb,                                              &
    ntype, natom, rat, adot,                                             &
    norbat, nqwf, delqwf, wvfao, lorb,                                   &
    psi, hpsi,                                                           &
    Hao, S, dh0drk,                                                      &
    vscr, kmscr,                                                         &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxddim, mxdorb,              &
    mxdscr, mxdlao)

! Written by Carlos Loia Reis adapting previous code.
! Modified 25 November 2015. JLM
! Modified    November 2018.  CLR (dh0drk computaion in ao's)
! Modified, documentation, 25 May 2020, JLM
! Modified, hamilt_pw, 6 June 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Modified, nanlspin, 30 November 2023. JLM
! Modified, allocation of d2h0drk2, 14 January 2024. JLM
! Modified ao_atomic_orbital, 6 October 2024. JLM
! Modified, removed mxdbnd, hdiag, etc... from API. 1 November 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdscr                          !<  array dimension of vscr
  integer, intent(in)                ::  mxdorb                          !<  array dimension for fft transform workspace
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
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

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)

! output

  integer, intent(inout)             ::  nbaslcao                        !<  number of eigenvectors
  complex(REAL64), intent(out)       ::  psi(mxddim,mxdorb)              !<  component j of eigenvector i
  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdorb)             !<  component j of eigenvector i

  complex(REAL64), intent(out)       ::  Hao(mxdorb,mxdorb)              !<  hamiltonian in non-orthogonal atomic orbitals
  complex(REAL64), intent(out)       ::  S(mxdorb,mxdorb)                !<  overlap matrix of non-orthogonal atomic orbitals

  complex(REAL64), intent(out)       ::  dh0drk(mxdorb,mxdorb,3)         !<  d <Psi|H|Psi> d k (computed with atomic orbitals)

! local allocatable arrays

  real(REAL64), allocatable          ::  xnlkb(:)
  complex(REAL64), allocatable       ::  anlga(:,:)

  complex(REAL64), allocatable       ::  h0(:,:)                         !  <Psi|H|Psi> without spin-orbit
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               !  d^2 <Psi|H|Psi> d k^2
  real(REAL64), allocatable          ::  ev_fake(:)
  integer, allocatable               ::  infolcao(:,:)                   !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)

  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i

! local varaibles

  integer       ::  mtxd                                                 !  dimension of the hamiltonian
  integer       ::  nder
  logical       ::  lkplusg                                              !  If true use the previous G-vectors (same mtxd and isort)

  logical       ::  lnewanl                                              !  indicates that anlga has been recalculated (not used in default implementation)       integer    ::  mxdanl, mxdsml
  integer       ::  mxdanl
  integer       ::  nanl, nanlso, nanlspin
  real(REAL64)  ::  bdot(3,3),vcell,veffr1

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)


  allocate(h0(mxdorb,mxdorb))
  allocate(ev_fake(mxdorb))


  if(flgpsd /= 'PSEUKB') then
    write(6,*)
    write(6,'("   stopped in ao_h_and_s: wrong pseudo ",a6)') flgpsd

    stop

  endif

  call adot_to_bdot(adot, vcell, bdot)

  veffr1 = real(veff(1),REAL64)


  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
      mxdtyp)

  mxdanl = nanl
  allocate(xnlkb(mxdanl))
  allocate(anlga(mxddim,mxdanl))

  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))

  lkplusg = .FALSE.

  call hamilt_pw(emax, rkpt, lkplusg, veffr1, nanl,                      &
      mtxd, isort, qmod, ekpg, hdiag,                                    &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  lnewanl = .TRUE.

  allocate(infolcao(5,mxdorb))

  call atomic_orbital_c16(rkpt, mtxd, isort, icmplx,                     &
      nbaslcao, psi, infolcao,                                           &
      ng, kgv,                                                           &
      norbat, nqwf, delqwf, wvfao, lorb,                                 &
      ntype, natom, rat, adot,                                           &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve, mxdlao)

  deallocate(infolcao)

  call zgemm('c','n', nbaslcao, nbaslcao, mtxd, C_UM, psi, mxddim,       &
      psi, mxddim, C_ZERO, S, mxdorb)

  call hk_psi_c16(mtxd, nbaslcao, psi, hpsi, lnewanl,                    &
      ng, kgv,                                                           &
      ekpg, isort, vscr, kmscr,                                          &
      anlga, xnlkb, nanl,                                                &
      mxddim, mxdorb, mxdanl, mxdgve, mxdscr)

  nder = 1

  allocate(d2h0drk2(1,1,3,3))

  call kdotp_matrix(mtxd, nbaslcao, psi, ev_fake, rkpt, isort, nder,     &
      h0, dh0drk, d2h0drk2,                                              &
      ng, kgv,                                                           &
      ntype, natom, rat, adot,                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      mxdtyp,mxdatm,mxdlqp,mxddim,mxdorb,mxdgve)

  deallocate(d2h0drk2)

  call zgemm('c','n', nbaslcao, nbaslcao, mtxd, C_UM, hpsi, mxddim,      &
      psi, mxddim, C_ZERO, Hao, mxdorb)

  deallocate(xnlkb)
  deallocate(anlga)

  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)

  deallocate(h0)
  deallocate(ev_fake)

  return

end subroutine ao_h_and_s
