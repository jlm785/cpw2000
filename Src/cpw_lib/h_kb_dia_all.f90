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

!>  calculates the hamiltonian for one k-point and diagonalizes
!>  either in a plane wave basis basis
!>  or in an LCAO basis.
!>
!>
!>  \author       Carlos Loia reis, Jose Luis Martins
!>  \version      5.09
!>  \date         May 2020. 11 November 2023.
!>  \copyright    GNU Public License v2

subroutine h_kb_dia_all(diag_type, emax, rkpt, neig, nocc,               &
  flgpsd, ipr, ifail, icmax, iguess, epspsi,                             &
  ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                        &
  sfact, veff, icmplx,                                                   &
  nqnl, delqnl, vkb, nkb,                                                &
  ntype, natom, rat, adot,                                               &
  mtxd, hdiag, isort, qmod, ekpg, lkpg,                                  &
  psi, hpsi, ei,                                                         &
  vscr, kmscr,                                                           &
  latorb, norbat, nqwf, delqwf, wvfao, lorb,                             &
  mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,                &
  mxdbnd, mxdscr, mxdlao)

! Adapted May 2020. CLR
! Modified, latorb, 7 June 2020. JLM
! Added nocc for future use. 12 June 2020. JLM
! Modified, norbtot bug, 30 November 2020. JLM
! Modified, k-point far from 1st BZ, 30 January 2021. JLM
! default value of ifail. 11 November 2023. JLM

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
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  real(REAL64), intent(in)           ::  rkpt(3)                         !<  component in lattice coordinates of the k-point
  integer, intent(in)                ::  nocc                            !<  number of ocupied eigenvectors for strict convergence (in reserve for later use).
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

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k


  real(REAL64)                       ::  veffr1


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

! allocatable local arrays

  integer, allocatable               ::  isort_tr(:)

! local variables

  integer        ::  iguess_local
  logical        ::  lhpsi
  real(REAL64)   ::  xmax
  integer        ::  norbtot

  logical        ::  lkshift                                             !  shift is needed
  real(REAL64)   ::  qualpsi                                             !  estimate of the quality of psi after applying the shift.  1 good, 0 terrible.

  real(REAL64)   ::  rk_tmp(3)
  integer        ::  kgshift(3)

  integer        ::  mtxd_tr

! parameters

  real(REAL64), parameter    ::  UM = 1.0_REAL64

! counter

  integer        ::  i

! AVOIDS WARNINGS UNTIL NOCC IS IMPLEMENTED
  IFAIL = NOCC

  ifail = 0


! deals with k-points far away from the 1st Brillouin zone

  allocate(isort_tr(mxddim))

  call h_kb_dia_check_k(emax, rkpt,                                      &
      lkshift, qualpsi, mtxd_tr, isort_tr,                               &
      ng, kgv, adot, kmscr,                                              &
      mxdgve, mxddim)

  rk_tmp(:) = rkpt(:)
  if(lkshift) then
    write(6,*)
    write(6,'("   WARNING    in h_kb_dia_all.  k-point ", 3f10.3,        &
       &  "   was shifted")') (rkpt(i),i=1,3)
    write(6,'("   Energies will be accurate, quality of wave-functions", &
       &  "   is (0--1): ",f10.3)') qualpsi
    write(6,*)
    do i = 1,3
      kgshift(i) = nint(rkpt(i))
      rk_tmp(i) = rkpt(i) - UM*kgshift(i)
    enddo
  endif


  lhpsi = .FALSE.

  iguess_local = 0

  veffr1 = real(veff(1),REAL64)

  if(diag_type == "pw  ") then

    norbtot = 0
    if(latorb) then
      call size_nbaslcao(ntype, natom, norbat, lorb, norbtot,             &
           mxdtyp, mxdlao)
    endif
    iguess_local = iguess

    if(iguess == 0 .and. latorb .and. neig <= norbtot) then

      lhpsi = .TRUE.

      call h_kb_dia_ao(emax, rk_tmp, neig, flgpsd, 'AOJC  ',             &
          veffr1, icmplx, lhpsi,                                         &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          norbat, nqwf, delqwf, wvfao, lorb,                             &
          ntype, natom, rat, adot,                                       &
          mtxd, hdiag, isort, qmod, ekpg, lkpg,                          &
          psi, hpsi, ei,                                                 &
          vscr, kmscr,                                                   &
          mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdbnd, mxdscr, mxdlao)

      iguess_local = 1


    endif

    call h_kb_dia(emax, rk_tmp, neig, flgpsd,                            &
        ipr, ifail, icmax, iguess_local, epspsi,                         &
        ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                  &
        sfact, veff, icmplx,                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, lkpg,                            &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim, mxdbnd,  &
        mxdscr)

    if(ifail /= 0) then
      call ditsp_error(xmax, neig, mtxd, psi, hpsi,                      &
          mxddim, mxdbnd)

      if(xmax > 4.0*epspsi*epspsi) then
        write(6,*)
        write(6,*) '  WARNING       WARNING:   After h_kb_dia '
        write(6,*) '  The estimated error in energy has an accuracy'
        write(6,'("  of ",f8.1," digits")') -log10(xmax)
        write(6,*)

        if(ifail < 4) then
          write(6,*)
          write(6,*) '  STOPPED   in h_kb_dia_all.'
          write(6,*) '  sqrt(xmax), epspsi = ', sqrt(xmax), epspsi
          write(6,*)

          stop

        endif

      endif

    endif

  elseif(diag_type == "aojc") then

    call h_kb_dia_ao(emax, rk_tmp, neig, flgpsd, 'AOJC  ',               &
        veffr1, icmplx, lhpsi,                                           &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        norbat, nqwf, delqwf, wvfao, lorb,                               &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, lkpg,                            &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdbnd, mxdscr, mxdlao)


  elseif(diag_type == "ao  ") then

    call h_kb_dia_ao(emax, rk_tmp, neig, flgpsd, 'AO    ',               &
        veffr1, icmplx, lhpsi,                                           &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        norbat, nqwf, delqwf, wvfao, lorb,                               &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, lkpg,                            &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        mxdtyp, mxdatm, mxdgve, mxdlqp, mxddim, mxdbnd, mxdscr, mxdlao)

  endif

  if(lkshift) then
    call psi_translate(neig, kgshift,                                    &
        mtxd, isort, psi, mtxd_tr, isort_tr,                             &
        kgv,                                                             &
        mxddim, mxdbnd, mxdgve)
  endif

  deallocate(isort_tr)

  return
end subroutine h_kb_dia_all
