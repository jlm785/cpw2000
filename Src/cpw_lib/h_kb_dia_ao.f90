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

!>  Calculates the hamiltonian for one k-point in the atomic basis
!>  and diagonalizes it in an LCAO basis.  It can do also a Jacobian update improvement.
!>
!>  \author       José Luís Martins
!>  \version      5.09
!>  \date         18 october 1993. 30 November 2023.
!>  \copyright    GNU Public License v2



subroutine h_kb_dia_ao(emax, rkpt, neig, flgpsd, flgscf,                 &
    veffr1, icmplx, lhpsi,                                               &
    ng, kgv,                                                             &
    nqnl, delqnl, vkb, nkb,                                              &
    norbat, nqwf, delqwf, wvfao, lorb,                                   &
    ntype, natom, rat, adot,                                             &
    mtxd, hdiag, isort, qmod, ekpg, lkpg,                                &
    psi, hpsi, ei,                                                       &
    vscr, kmscr,                                                         &
    mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdbnd, mxdscr, mxdlao)

! version 4.0. 18 october 1993. jlm
! modified 25 march 1999, 19 april 1999.
! modified (chd) 24 july 2002.
! modified 19 november 2013. veffr(1). jlm
! modified 19 december 2013. kmscr. jlm
! modified December 19, 2013 (f90). jlm
! modified, dimensions vkb, March 31, 2014. jlm
! modified 5 May 2014. rq_jacr. EI. JLM
! modified 14 May 2014. xveci = 0. JLM
! modified complex version 19 September 2015. JLM
! Modified kmscr, 28 October 2015. JLM
! Modified lkpg, August 2019. JLM
! Modified documentation, August 2019. JLM
! Modified, corrected dimension ei. 18 February 2020. JLM
! Modified, hamilt_pw, diag_c16_gen, hk_from_bas_diag, 2 June 2020. JLM
! Modified, qmod-->ekpg in hk_psi. 13 February 2021. JLM
! Modified, nanlspin, 30 November 2023. JLM
! prefix diag_rq_jac. 17 March 2024. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands
  integer, intent(in)                ::  mxdscr                          !<  array dimension for vscr
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  character(len=6), intent(in)       ::  flgscf                          !<  type of scf calculation

  real(REAL64), intent(in)           ::  veffr1                          !<  Average value (veff(1)) of the effective potential.
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex
  logical, intent(in)                ::  lhpsi                           !<  indicates that H | psi> should be also calculated

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * kb nonlocal pseudo. for atom k, ang. mom. l. NOT normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   kb pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  real(REAL64), intent(in)           ::  vscr(mxdscr)                    !<  screened potential in the fft real space mesh
  integer, intent(in)                ::  kmscr(7)                        !<  max value of kgv(i,n) used for the potential fft mesh and fft mesh size

  logical, intent(in)                ::  lkpg                            !<  If true use the previous G-vectors (same mtxd and isort)

! input and output

  integer, intent(inout)             ::  neig                            !<  number of eigenvectors (requested on input, modified by degeneracies on output)

  integer, intent(inout)             ::  mtxd                            !<  dimension of the hamiltonian
  integer, intent(inout)             ::  isort(mxddim)                   !<  g-vector associated with row/column i of hamiltonian

! output

  real(REAL64), intent(out)          ::  ei(mxdbnd)                      !<  eigenvalues (Hartree)
  complex(REAL64), intent(out)       ::  psi(mxddim,mxdbnd)              !<  | psi > component j of eigenvector i
  complex(REAL64), intent(out)       ::  hpsi(mxddim,mxdbnd)             !<  H | psi >
  real(REAL64), intent(out)          ::  hdiag(mxddim)                   !<  hamiltonian diagonal
  real(REAL64), intent(out)          ::  qmod(mxddim)                    !<  length of k+g-vector of row/column i
  real(REAL64), intent(out)          ::  ekpg(mxddim)                    !<  kinetic energy (hartree) of k+g-vector of row/column i

! local allocatable arrays

  real(REAL64), allocatable           ::  xnlkb(:)                       !  KB projector normalization
  complex(REAL64), allocatable        ::  anlga(:,:)                     !  KB projectors
  complex(REAL64), allocatable        ::  bas(:,:)                       !  | bas >  atomic-orbital wave-functions
  complex(REAL64), allocatable        ::  hbas(:,:)                      !  H | bas >
  integer, allocatable                ::  infolcao(:,:)                  !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  complex(REAL64), allocatable        ::  psibas(:,:)                    !  wavefunctions in atomic basis

! local variables

  integer       ::  mxdorb, mxdanl

  integer       ::  nbasorb
  integer       ::  nanl, nanlso, nanlspin, ndeg

  logical       ::  lnewanl                                              !  indicates that anlga has been recalculated (not used in default implementation)
  logical       ::  lhpsiloc                                             !  local value of lhpsi

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64

! counters

  integer       ::  n


  if(flgpsd /= 'PSEUKB') then
    write(6,*)
    write(6,'("   STOPPED in h_kb_dia_ao   wrong pseudo ",a6)') flgpsd

    stop

  endif

  if(flgscf /= 'AO    ' .and. flgscf /= 'AOJC  ' .and.                   &
     flgscf /= 'AOJCPW') then
    write(6,*)
    write(6,'("   STOPPED in h_kb_dia_ao: not an AO calculation ", a6)') &
          flgscf

    stop

  endif

  call size_proj_nl_kb(ntype, natom, nkb, nanl, nanlso, nanlspin,        &
       mxdtyp)

  mxdanl = nanl

  allocate(xnlkb(mxdanl))
  allocate(anlga(mxddim,mxdanl))

  call hamilt_pw(emax, rkpt, lkpg, veffr1, nanl,                         &
      mtxd, isort, qmod, ekpg, hdiag,                                    &
      ng, kgv,                                                           &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      anlga, xnlkb,                                                      &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdanl, mxdgve)

  lnewanl = .TRUE.

  call size_nbaslcao(ntype, natom, norbat, lorb, mxdorb, mxdtyp, mxdlao)

  allocate(bas(mxddim,mxdorb))
  allocate(hbas(mxddim,mxdorb))

  allocate(infolcao(5,mxdorb))

  call atomic_orbital_c16(rkpt, mtxd, isort, icmplx,                     &
      nbasorb, bas, infolcao,                                            &
      ng, kgv,                                                           &
      norbat, nqwf, delqwf, wvfao, lorb,                                 &
      ntype, natom, rat, adot,                                           &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve, mxdlao)

  deallocate(infolcao)

  if(neig > nbasorb) then
    write(6,'("   STOPPED in h_kb_dia_ao   number of orbitals ",         &
      & "is smaller than number of eigenvectors",2i8)') nbasorb,neig

    stop

  endif

  call hk_psi_c16(mtxd, nbasorb, bas, hbas,  lnewanl,                    &
      ng, kgv,                                                           &
      ekpg, isort, vscr, kmscr,                                          &
      anlga, xnlkb, nanl,                                                &
      mxddim, mxdorb, mxdanl, mxdgve, mxdscr)

  lhpsiloc = lhpsi
  if(flgscf == 'AOJC  ' .or. flgscf == 'AOJCPW') lhpsiloc = .TRUE.

! do not use it

  allocate(psibas(1,1))

  call hk_from_bas_diag(mtxd, nbasorb, neig, ndeg,                       &
      .FALSE., .FALSE., .TRUE., lhpsiloc,                                &
      bas, hbas, ei, psibas, psi, hpsi,                                  &
      mxddim, mxdorb, mxdbnd)

! improves the vectors from the previous calculation

  if(flgscf == 'AOJC  ' .or. flgscf == 'AOJCPW') then

    deallocate(bas)
    deallocate(hbas)
    allocate(bas(mxddim,2*ndeg))
    allocate(hbas(mxddim,2*ndeg))

    do n = 1,ndeg
      call zcopy(mtxd, psi(:,n), 1, bas(:,ndeg+n), 1)
    enddo

    do n = 1,ndeg
      call zcopy(mtxd, hpsi(:,n), 1, hbas(:,ndeg+n), 1)
    enddo

    do n = 1,ndeg
      call zaxpy(mtxd, cmplx(-ei(n), ZERO, REAL64), psi(:,n), 1, hpsi(:,n), 1)
    enddo

    call diag_rq_jac_c16(hpsi, ei, hdiag, mtxd, ndeg, mxddim, mxdbnd)

    do n = 1,ndeg
      call zcopy(mtxd, hpsi(:,n), 1, bas(:,n), 1)
    enddo

    call hk_psi_c16(mtxd, ndeg, bas, hbas, lnewanl,                      &
        ng, kgv,                                                         &
        ekpg, isort, vscr, kmscr,                                        &
        anlga, xnlkb, nanl,                                              &
        mxddim, 2*ndeg, mxdanl, mxdgve, mxdscr)

    call hk_from_bas_diag(mtxd, 2*ndeg, neig, ndeg,                      &
        .FALSE., .FALSE., .TRUE., lhpsi,                                 &
        bas, hbas, ei, psibas, psi, hpsi,                                &
        mxddim, 2*ndeg, mxdbnd)

  endif

  deallocate(bas)
  deallocate(hbas)

  deallocate(xnlkb)
  deallocate(anlga)

  deallocate(psibas)

  return

end subroutine h_kb_dia_ao
