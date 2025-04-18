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

!>  Calculates the quantum geometric properties for a given k-vector
!>
!>  \author       Jose Luis Martins
!>  \version      5.11.
!>  \date         4 April 2024.
!>  \copyright    GNU Public License v2

subroutine out_qgeom(ioreplay,                                           &
    emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot, efermi,           &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

! Should be merged with out_mass_berry someday...

! Adapted 4 April 2024 from out_mass_berry. JLM
! Orientation. 11 April 2024. JLM


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

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  integer, intent(in)                ::  iguess                          !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64)                       ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
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

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi>

  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  real(REAL64), allocatable          ::  deidxk(:)                       !  derivative of energy
  real(REAL64), allocatable          ::  d2eidxk2(:)                     !  second derivative of energy

  complex(REAL64), allocatable       ::  dhdkpsi(:,:,:)                  !  d H d k | psi>
  complex(REAL64), allocatable       ::  dpsidk(:,:,:)                   !  d |psi> / d k
  complex(REAL64), allocatable       ::  psidhdkpsi(:,:,:,:)             !  <psi_n| d H / d k |psi_m> for each energy level (lattice coordinates)

  complex(REAL64), allocatable       ::  tfqg(:,:,:,:,:)                 !<  Tensor F quantum geomtric, curvature and metric, (lattice coordinates)
  complex(REAL64), allocatable       ::  tgammamf(:,:,:,:,:)             !<  Tensor Gamma-F, orbital magnetization and contribution to effective mass (lattice coordinates)
  complex(REAL64), allocatable       ::  td2hdk2(:,:,:,:,:)              !<  <psi| d^2 H / d k^2 |psi> (lattice coordinates)

  real(REAL64), allocatable          ::  ei_pert(:)                      !  eigenvalue no. i. (hartree) spin-orbit perturbation
  real(REAL64), allocatable          ::  ei_sp(:)                        !  eigenvalue no. i. (hartree)

  complex(REAL64), allocatable       ::  psi_sp(:,:)                     !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hpsi_sp(:,:)                    !  H | psi>

  real(REAL64), allocatable          ::  vscr_sp(:,:)                    !  screened potential (spin)

  complex(REAL64), allocatable       ::  hamk_so(:,:)                    !  hamiltonian matrix with spin-orbit
  real(REAL64), allocatable          ::  eiall_so(:)                     !  energies
  complex(REAL64), allocatable       ::  vecall_so(:,:)                  !  wave-functions

  complex(REAL64), allocatable       ::  pauli(:,:)                      !  <psi|sigma_xyz|psi>

! local variables

  integer           ::  mxdscr                                           !  array dimension for screening potential
  integer           ::  mxdlev                                           !  array dimension for number of levels
  integer           ::  mxddeg                                           !  array dimension for number of levels

  integer           ::  ifail                                            !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  integer           ::  mxddim                                           !  array dimension for the hamiltonian
  integer           ::  mxdbnd                                           !  array dimension for the number of bands
  integer           ::  mxdwrk                                           !  array dimension for fft transform workspace

  integer           ::  mtxd                                             !  dimension of the hamiltonian
  integer           ::  neig                                             !  number of eigenvectors
  integer           ::  kmscr(7)                                         !  max value of kgv(i,n) used for the potential fft mesh
  integer           ::  idshift                                          !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)      ::  vmax, vmin                                       !  maximum and minimum values of vscr

  integer           ::  ipr
  integer           ::  nsfft(3)
  real(REAL64)      ::  rkcar(3)

  real(REAL64)      ::  rkpt(3)
  character(len=20) ::  typeofk

  character(len=1)  ::  yesno_so

  real(REAL64)      ::  bdot(3,3), vcell

  integer           ::  nlevel              !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg              !  maximum number of degeneracies

  integer           ::  neigin              !  initial value of neig
  integer           ::  nocc

  integer           ::  nsp, mxdnsp
  integer           ::  njac, nritz

  integer           ::  imethod, imethod2
  integer           ::  info

! constants

  real(REAL64), parameter     ::  TOL = 1.0E-8_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer    ::  j, n, k


! hard coded values

  njac = 1
  nritz = 120

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

  write(6,*)
  write(6,'(" enter number of bands (greater than ~",i4,")")') nint(ztot/2)
  read(5,*) neig
  write(ioreplay,*) neig,'   initial desired number of bands'

! allow for degeneracies

  neigin = neig+5
  mxdbnd = neigin+1

! gets the k-point

  call adot_to_bdot(adot,vcell,bdot)

  typeofk = 'reference k-point'

  call cpw_pp_get_k_vector(rkpt, rkcar, adot, typeofk, ioreplay)

  write(6,*)
  write(6,*) '  coordinates of the chosen k-point:'
  write(6,*) '      lattice coord.                  cartesian coord.'
  write(6,*)
  write(6,'(4x,3f9.4,5x,3f9.4)') (rkpt(j),j=1,3), (rkcar(j),j=1,3)
  write(6,*)

! Tries to get good suggestions for the k.p model


  call size_mtxd(emax, rkpt, adot, ng, kgv, mxddim)

  mxddim = mxddim + 30

! finds mxddim, mxdbnd

! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))

  nocc = neigin

  ipr = 1

  call h_kb_dia_all('pw  ', emax, rkpt, neigin, nocc,                    &
      flgpsd, ipr, ifail, icmax, iguess, epspsi,                         &
      ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                    &
      sfact, veff, icmplx,                                               &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                           &
      psi, hpsi, ei,                                                     &
      vscr, kmscr,                                                       &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,            &
      mxdbnd, mxdscr, mxdlao)

! tries to "align" eigenvectors

  call psi_orient_xyz(rkpt, adot, mtxd, neig, isort,                     &
      psi, ei,                                                           &
      ng, kgv,                                                           &
      mxddim, mxdbnd, mxdgve)

! energy levels.  first finds dimensions. recalculates neig according to degeneracies.

  write(6,*)
  write(6,*) '   Do you wan the results with spin-orbit (y/n)'
  write(6,*)

  read(5,*) yesno_so

  write(ioreplay,*) yesno_so,'   with spin-orbit'

  if(yesno_so == 'y' .or. yesno_so == 'Y') then

!   choose method for spin-orbitals

    write(6,*)
    write(6,*) '   How do you wanto to treat the reference wave-functions?'
    write(6,*)
    write(6,*) '   (1) Full diagonalization of the spin-orbit hamiltonian'
    write(6,*) '       Recommended if you can afford to diagonalize a matrix'
    write(6,'("       of size ",i6 " by ",i6)') 2*mtxd,2*mtxd
    write(6,*) '   (2) Iterative diagonalization of the spin-orbit hamiltonian'
    write(6,*) '   (3) First order perturbation in spin-orbit'
    write(6,*)
    write(6,*) '   Enter your choice 1-3'
    write(6,*)

    read(5,*) imethod

    write(ioreplay,*) imethod,'   spin-orbit method'

    if(imethod < 1 .or. imethod > 3) then
      write(6,*)
      write(6,*) '   Wrong value, enter method again 1,2,3'
      write(6,*)

      read(5,*) imethod

      write(ioreplay,*) imethod2,'   spin-orbit method again'

      if(imethod < 1 .or. imethod > 3) then
        if(mtxd < 1000) then
          write(6,*) '   Wrong value again, using full diagonalization'
          imethod = 1
        else
          write(6,*) '   Wrong value again, using iterative diagonalization'
          imethod = 2
        endif
      endif

    endif

    allocate(ei_pert(2*mxdbnd))
    allocate(ei_sp(2*mxdbnd))
    allocate(psi_sp(2*mxddim,2*mxdbnd))
    allocate(hpsi_sp(2*mxddim,2*mxdbnd))

    nsp = 1
    mxdnsp = 1

    allocate(vscr_sp(mxdscr,mxdnsp))

    vscr_sp(:,1) = vscr(:)

    if(imethod == 1) then

      allocate(hamk_so(2*mxddim,2*mxddim))
      allocate(eiall_so(2*mxddim))
      allocate(vecall_so(2*mxddim,2*mxddim))

      call hamilt_kb_alt_so(rkpt, mtxd, isort, qmod, ekpg,               &
           hamk_so,                                                      &
           ng, kgv, phase, conj, inds, kmax, indv,                       &
           veff, nqnl, delqnl, vkb, nkb,                                 &
           ntype, natom, rat, adot,                                      &
           mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim, 2*mxddim)

      call diag_c16(2*mtxd, hamk_so, eiall_so, vecall_so, 2*mxddim, info)

      if(info /= 0) then
        write(6,*)
        write(6,*) '   Stopped in out_mass_berry.  Failed diagonalization'
        write(6,*) '   info = ',info

        stop

      endif

      do n = 1,2*neig
        ei_sp(n) = eiall_so(n)
        psi_sp(:,n) = vecall_so(:,n)
      enddo

      deallocate(hamk_so)
      deallocate(eiall_so)
      deallocate(vecall_so)

    elseif(imethod == 2) then

      call diag_improve_psi_spin(rkpt, mtxd, neig, njac, nritz, TOL,     &
          ei, psi,                                                       &
          ei_pert, ei_sp, psi_sp, hpsi_sp,                               &
          ng, kgv,                                                       &
          ekpg, isort, vscr_sp, kmscr, nsp,                              &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr, mxdnsp)

    elseif(imethod == 3) then

      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
          neig, psi, ei, ei_pert, psi_sp, .TRUE.,                        &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      ei_sp(:) = ei_pert(:)

    else

      write(6,*) '   Stopped in out_qgeom.  paranoid check!'

      stop

    endif

! tries to "align" eigenvectors with spin

    call psi_orient_spin(mtxd, 2*neig, psi_sp, ei_sp,                    &
        mxddim, mxdbnd)

!   calculates the spin polarization

    allocate(pauli(3, 2*mxdbnd))

    call psi_pauli_psi(mtxd, 2*neig, psi_sp, pauli,                      &
        mxddim, mxdbnd)

!   now finds degeneracies

    allocate(levdeg(1))
    allocate(leveigs(1,1))

    neig = 2*neig
    neigin = neig

    call berry_degeneracy(.TRUE., neigin, neig, ei_sp, TOL,              &
           nlevel, maxdeg, levdeg, leveigs,                              &
           2*mxdbnd, 1, 1)

    mxdlev = nlevel
    mxddeg = maxdeg

    deallocate(levdeg)
    deallocate(leveigs)

    allocate(levdeg(mxdlev))
    allocate(leveigs(mxdlev,mxddeg))

!   fills the information

    call berry_degeneracy(.FALSE., neigin, neig, ei_sp, TOL,             &
           nlevel, maxdeg, levdeg, leveigs,                              &
           2*mxdbnd, mxdlev, mxddeg)

    write(6,*)
    write(6,'("   The system has ",i5," levels")') nlevel
    write(6,*)
    write(6,'(25x," -<sigma_x>/2        -<sigma_y>/2        -<sigma_z>/2")')
    write(6,*)
    do n = 1,nlevel
      write(6,*)
      write(6,'("  Level ",i5," has degeneracy ",i3)') n, levdeg(n)
      write(6,*)
      do j = 1,levdeg(n)
        write(6,'(f12.6,5x,3(5x,2f8.3))') ei_sp(leveigs(n,j))*HARTREE,    &
              (-pauli(k,leveigs(n,j))/2,k=1,3)
      enddo
    enddo
    write(6,*)

    deallocate(pauli)

!   get berry quantities.  starts with allocations.

    allocate(dhdkpsi(2*mxddim,2*mxdbnd,3))
    allocate(dpsidk(2*mxddim,2*mxdbnd,3))

    allocate(psidhdkpsi(mxddeg,mxddeg,3,mxdlev))

    allocate(tfqg(3,3,mxddeg,mxddeg,mxdlev))
    allocate(tgammamf(3,3,mxddeg,mxddeg,mxdlev))
    allocate(td2hdk2(3,3,mxddeg,mxddeg,mxdlev))

    call berry_derivative_spin(rkpt, mtxd, neig, isort, ekpg, .TRUE.,    &
      nlevel, levdeg, leveigs,                                           &
      psi_sp, ei_sp,                                                     &
      dhdkpsi, dpsidk, psidhdkpsi, tfqg, tgammamf, td2hdk2,              &
      ng, kgv,                                                           &
      vscr_sp, kmscr, nsp,                                               &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mxdtyp, mxdatm, mxdlqp, mxddim, 2*mxdbnd, mxdgve, mxdscr,          &
      mxdlev, mxddeg, mxdnsp)

    deallocate(ei_pert)
    deallocate(psi_sp)
    deallocate(hpsi_sp)
    deallocate(vscr_sp)

    allocate(deidxk(2*mxdbnd))
    allocate(d2eidxk2(2*mxdbnd))

  else

    allocate(levdeg(1))
    allocate(leveigs(1,1))

    call berry_degeneracy(.TRUE., neigin, neig, ei, TOL,                 &
           nlevel, maxdeg, levdeg, leveigs,                              &
           mxdbnd, 1, 1)

    mxdlev = nlevel
    mxddeg = maxdeg

    deallocate(levdeg)
    deallocate(leveigs)

    allocate(levdeg(mxdlev))
    allocate(leveigs(mxdlev,mxddeg))

!   fills the information

    call berry_degeneracy(.FALSE., neigin, neig, ei, TOL,                &
           nlevel, maxdeg, levdeg, leveigs,                              &
           mxdbnd, mxdlev, mxddeg)

    write(6,*)
    write(6,'("   The system has ",i5," levels")') nlevel
    do n = 1,nlevel
      write(6,*)
      write(6,'("  Level ",i5," has degeneracy ",i3)') n, levdeg(n)
      write(6,*)
      do j = 1,levdeg(n)
        write(6,'(f12.6)') ei(leveigs(n,j))*HARTREE
      enddo
    enddo
    write(6,*)

!   get berry quantities.  starts with allocations.

    allocate(dhdkpsi(mxddim,mxdbnd,3))
    allocate(dpsidk(mxddim,mxdbnd,3))

    allocate(psidhdkpsi(mxddeg,mxddeg,3,mxdlev))

    allocate(tfqg(3,3,mxddeg,mxddeg,mxdlev))
    allocate(tgammamf(3,3,mxddeg,mxddeg,mxdlev))
    allocate(td2hdk2(3,3,mxddeg,mxddeg,mxdlev))

    call berry_derivative(rkpt, mtxd, neig, isort, ekpg, .TRUE.,         &
      nlevel, levdeg, leveigs,                                           &
      psi, ei,                                                           &
      dhdkpsi, dpsidk, psidhdkpsi, tfqg, tgammamf, td2hdk2,              &
      ng, kgv,                                                           &
      vscr, kmscr,                                                       &
      nqnl, delqnl, vkb, nkb,                                            &
      ntype, natom, rat, adot,                                           &
      mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve, mxdscr,            &
      mxdlev, mxddeg)

    allocate(deidxk(mxdbnd))
    allocate(d2eidxk2(mxdbnd))

  endif

  if(yesno_so == 'y' .or. yesno_so == 'Y') then

    call out_qgeom_print(ioreplay, nlevel, levdeg, leveigs,              &
        ei_sp, adot, efermi,                                             &
        tfqg, tgammamf, td2hdk2,                                         &
        2*mxdbnd, mxdlev, mxddeg)

  else

    call out_qgeom_print(ioreplay, nlevel, levdeg, leveigs,              &
        ei, adot, efermi,                                                &
        tfqg, tgammamf, td2hdk2,                                         &
        mxdbnd, mxdlev, mxddeg)

  endif

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)

  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)

  deallocate(psi)
  deallocate(hpsi)

  deallocate(levdeg)
  deallocate(leveigs)

  deallocate(dhdkpsi)
  deallocate(dpsidk)
  deallocate(psidhdkpsi)

  deallocate(tfqg)
  deallocate(tgammamf)
  deallocate(td2hdk2)

  deallocate(deidxk)
  deallocate(d2eidxk2)

  if(yesno_so == 'y' .or. yesno_so == 'Y') deallocate(ei_sp)

  return

end subroutine out_qgeom
