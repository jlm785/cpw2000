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

!>  Performs the interpolation for the calculation the effective mass
!>  for a given k-vector and direction with finite differences
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      5.11
!>  \date         7 November 2023 .5 March 2024
!>  \copyright    GNU Public License v2

subroutine out_mass_fd_xk(rkpt, xk, neig, npt, delta, lso,               &
    deidk, d2eidk2,                                                      &
    ei_so, deidk_so, d2eidk2_so,                                         &
    emax, flgdal, flgpsd, epspsi, icmax,                                 &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao,              &
    mxddim, mxdbnd)

! January 2023. JLM
! Modified, spin_perturb to spin_improve. 5 March 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type
  integer, intent(in)                ::  mxddim                          !<  array dimension of plane-waves
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for number of bands

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  xk(3)                           !<  k-direction in reciprocal lattice coordinates

  integer, intent(in)                ::  neig                            !<  wavefunction dimension

  integer, intent(in)                ::  npt                             !<  2*npt+1 is the total number of interpolation points
  real(REAL64), intent(in)           ::  delta                           !<  spacing between the poins used in the interpolation
  logical, intent(in)                ::  lso                             !<  calculates with spin-orbit in perturbation


  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization

   real(REAL64), intent(in)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,mxdgve)                   !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  phase(mxdgve)                   !<  phase factor of G-vector n
  real(REAL64), intent(in)           ::  conj(mxdgve)                    !<  is -1 if one must take the complex conjugate of x*phase

  integer, intent(in)                ::  ns                              !<  number os stars with length less than gmax
  integer, intent(in)                ::  inds(mxdgve)                    !<  star to which g-vector n belongs
  integer, intent(in)                ::  kmax(3)                         !<  max value of |kgv(i,n)|
  integer, intent(in)                ::  indv(mxdcub)                    !<  kgv(i,indv(jadd)) is the g-vector associated with jadd. jadd is defined by the g-vector components and kmax
  real(REAL64), intent(in)           ::  ek(mxdnst)                      !<  kinetic energy (hartree) of g-vectors in star j

  complex(REAL64), intent(in)        ::  sfact(mxdtyp,mxdnst)            !<  structure factor
  integer, intent(in)                ::  icmplx                          !<  indicates if the structure factor is complex

  complex(REAL64), intent(in)        ::  veff(mxdnst)                    !<  ionic potential (local+Hartree+XC) for the prototype g-vector in star j

  integer, intent(in)                ::  nqnl(mxdtyp)                    !<  number of points for the non-local pseudopotential interpolation
  real(REAL64), intent(in)           ::  delqnl(mxdtyp)                  !<  step used in the interpolation
  real(REAL64), intent(in)           ::  vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)  !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<  KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k


! output

  real(REAL64), intent(out)          ::  deidk(mxdbnd)                   !<  d E / d xk  (lattice coordinates)
  real(REAL64), intent(out)          ::  d2eidk2(mxdbnd)                 !<  d^2 E / d xk^2  (lattice coordinates)

  real(REAL64), intent(out)          ::  ei_so(2*mxdbnd)                 !<  eigenvalue no. i. (hartree)with perturbation spin-orbit
  real(REAL64), intent(out)          ::  deidk_so(2*mxdbnd)              !<  d E / d xk  with perturbation spin-orbit (lattice coordinates)
  real(REAL64), intent(out)          ::  d2eidk2_so(2*mxdbnd)            !<  d^2 E / d xk^2  with perturbation spin-orbit (lattice coordinates)

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi>
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)

  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  real(REAL64), allocatable          ::  ei_l(:,:)                       !  eigenvalue no. i. in the line (hartree)
  real(REAL64), allocatable          ::  rk_l(:,:)                       !  k-point on the line
  real(REAL64), allocatable          ::  ei_l_so(:,:)                    !  eigenvalue no. i. in the line (hartree) with so

  real(REAL64), allocatable          ::  ei_pert(:)                      !  eigenvalue in spin-peerturbation

  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i with so
  complex(REAL64), allocatable       ::  hpsi_so(:,:)                    !  H | psi>

  real(REAL64), allocatable          ::  vscr_sp(:,:)

  real(REAL64), allocatable          ::  xin(:), yin(:)

  integer, allocatable               ::  ipl(:)                          !  most similar levels.

! local variables

  integer           ::  mxdscr                                           !  array dimension for screening potential

  integer           ::  ifail                                            !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  integer           ::  mxdwrk                                           !  array dimension for fft transform workspace

  integer           ::  mtxd                                             !  dimension of the hamiltonian
  integer           ::  kmscr(7)                                         !  max value of kgv(i,n) used for the potential fft mesh
  integer           ::  idshift                                          !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)      ::  vmax, vmin                                       !  maximum and minimum values of vscr

  integer           ::  nrka
  character(len=5)  ::  labelk
  integer           ::  ipr
  integer           ::  nsfft(3)

  integer           ::  nocc
  integer           ::  iguess

  integer           ::  nsp, mxdnsp
  integer           ::  njac, nritz

  real(REAL64)      ::  yk(3), xknorm
  real(REAL64)      ::  vcell, bdot(3,3)

  real(REAL64)      ::  y(0:2), dy(0:2)

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-8_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-12_REAL64

! counters

  integer    ::  i, j, k, n


  if(npt < 1) then
    write(6,*)
    write(6,*) '    stopped in out_mass_fd_xk    npt = ', npt
    write(6,*)

    stop

  endif

  if(delta < TOL) then
    write(6,*)
    write(6,*) '    stopped in out_mass_fd_xk    delta = ', delta
    write(6,*)

    stop

  endif

! normalizes xk

  call adot_to_bdot(adot, vcell, bdot)

  xknorm = ZERO
  do i = 1,3
  do j = 1,3
    xknorm = xknorm + xk(i)*bdot(i,j)*xk(j)
  enddo
  enddo

  if(xknorm < EPS) then
    write(6,*)
    write(6,*) '   stopped in out_mass_fd_xk, xknorm = ', xknorm
    write(6,*)

    stop

  endif

  do i = 1,3
    yk(i) = xk(i) / sqrt(xknorm)
  enddo


! compatibility dummy argument

  iguess = 0

! calculates local potential in fft mesh


  if(flgdal == 'DUAL') then
    kmscr(1) = kmax(1)/2 + 2
    kmscr(2) = kmax(2)/2 + 2
    kmscr(3) = kmax(3)/2 + 2
  else
    kmscr(1) = kmax(1)
    kmscr(2) = kmax(2)
    kmscr(3) = kmax(3)
  endif

  call size_fft(kmscr, nsfft, mxdscr, mxdwrk)

  allocate(vscr(mxdscr))

  ipr = 1

  idshift = 0
  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
      ng, kgv, phase, conj, ns, inds,                                    &
      mxdscr, mxdgve, mxdnst)

! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))

  allocate(ei_l(mxdbnd,-npt:npt))
  allocate(rk_l(3,-npt:npt))

  if(lso) then
    allocate(ei_l_so(2*mxdbnd,-npt:npt))
    allocate(psi_so(2*mxddim,2*mxdbnd))
    allocate(hpsi_so(2*mxddim,2*mxdbnd))
    allocate(ei_pert(2*mxdbnd))
    mxdnsp = 1
    allocate(vscr_sp(mxdscr,mxdnsp))
  endif

  do n = -npt,npt

    do k = 1,3
      rk_l(k,n) = rkpt(k) - n*delta*yk(k)
    enddo

    nocc = neig
    iguess = 0

    ipr = 1

    call h_kb_dia_all('pw  ', emax, rk_l(:,n), neig, nocc,               &
        flgpsd, ipr, ifail, icmax, iguess, epspsi,                       &
        ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                  &
        sfact, veff, icmplx,                                             &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                         &
        psi, hpsi, ei_l(:,n),                                            &
        vscr, kmscr,                                                     &
        latorb, norbat, nqwf, delqwf, wvfao, lorb,                       &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,          &
        mxdbnd, mxdscr, mxdlao)

    ipr = 1

    nrka = -1
    call print_eig(ipr, 1, labelk, nrka, rk_l(:,n),                      &
        mtxd, icmplx, neig, psi,                                         &
        adot, ei_l(:,n), ekpsi, isort, kgv,                              &
        mxddim, mxdbnd, mxdgve)

    if(lso) then

      njac = 3
      nritz = 5

      nsp = 1

      vscr_sp(:,1) = vscr(:)

      call diag_improve_psi_spin(rk_l(:,n), mtxd, neig,                  &
          njac, nritz, TOL,                                              &
          ei_l(:,n), psi,                                                &
          ei_pert, ei_l_so(:,n), psi_so, hpsi_so,                        &
          ng, kgv,                                                       &
          ekpg, isort, vscr_sp, kmscr, nsp,                              &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

      ipr = 1

      nrka = -1
      call print_eig_so(ipr, 1, labelk, nrka, rk_l(:,n),                 &
          mtxd, neig, psi_so,                                            &
          adot, ei_l_so(:,n), ekpsi, isort, kgv,                         &
          mxddim, mxdbnd, mxdgve)

    endif

  enddo

! tries to match left and right sides

  allocate(ipl(neig))

  call out_mass_match(neig, npt, ei_l, ipl,                              &
       mxdbnd)

! does the fit

  allocate(xin(-npt:npt))
  allocate(yin(-npt:npt))

  do n = -npt,npt
    xin(n) = n*delta
  enddo

  do n = 1,neig

    do j = -npt,-1
      yin(j) = ei_l(ipl(n),j) - ei_l(n,0)
    enddo

    do j = 0,npt
      yin(j) = ei_l(n,j) - ei_l(n,0)
    enddo

    call poly_interp(y, dy, xin, yin, 2*npt, 2)

    d2eidk2(n) = y(2)
    deidk(n) = y(1)

  enddo

! repeats for spin-orbit

  if(lso) then

    deallocate(ipl)
    allocate(ipl(2*neig))

    call out_mass_match(2*neig, npt, ei_l_so, ipl,                       &
       2*mxdbnd)

  do n = 1,2*neig

    do j = -npt,-1
      yin(j) = ei_l_so(ipl(n),j) - ei_l_so(n,0)
    enddo

    do j = 0,npt
      yin(j) = ei_l_so(n,j) - ei_l_so(n,0)
    enddo

    call poly_interp(y, dy, xin, yin, 2*npt, 2)

    d2eidk2_so(n) = y(2)
    deidk_so(n) = y(1)
    ei_so(n) = y(0) + ei_l_so(n,0)

  enddo



  endif

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(hpsi)
  deallocate(ekpsi)

  deallocate(ei_l)
  deallocate(rk_l)
  if(lso) then
    deallocate(ei_l_so)
    deallocate(psi_so)
    deallocate(hpsi_so)
    deallocate(ei_pert)
    deallocate(vscr_sp)
  endif

  deallocate(xin,yin)
  deallocate(ipl)

  return

end subroutine out_mass_fd_xk





