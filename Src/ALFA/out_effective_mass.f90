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

!>  Calculates the effective mass for a given k-vector
!>  and direction (effective mass tensor for non-degenerate levels)
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         18 january 2022. 1 April 2023.
!>  \copyright    GNU Public License v2

subroutine out_effective_mass(ioreplay,                                  &
    emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,                   &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)



! Adapted from out_band_onek plus old "Silvaco" subroutines. 18 january 2022. JLM
! Better user interface, 1 April 2023. JLM

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

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  integer, intent(in)                ::  iguess                          !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64)                       ::  ztot                            !<  total charge density (electrons/cell)

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
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)

  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i
  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)


  complex(REAL64), allocatable       ::  h0(:,:)                         !  <Psi|H|Psi> without spin-orbit
  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               !  d^2 <Psi|H|Psi> d k^2

  complex(REAL64), allocatable       ::  hso0(:,:)                       !  <Psi|H|Psi> with spin-orbit
  complex(REAL64), allocatable       ::  dhso0drk(:,:,:)                 !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2hso0drk2(:,:,:,:)             !  d^2 <Psi|H|Psi> d k^2

  real(REAL64), allocatable          ::  ei_p(:), ei_m(:)                !  eigenvalue no. i. (hartree)

  real(REAL64), allocatable          ::  ei_so_p(:), ei_so_m(:)          !  eigenvalue no. i. (hartree)

  REAL(REAL64), allocatable          ::  xmass(:,:), xgrad(:)            !  effective mass and energy gradient

! local variables

  integer           ::  mxdscr                                           !  array dimension for screening potential

  integer           ::  ifail                                            !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  integer           ::  nder                                             !  order of derivative

  integer           ::  mxddim                                           !  array dimension for the hamiltonian
  integer           ::  mxdbnd                                           !  array dimension for the number of bands
  integer           ::  mxdwrk                                           !  array dimension for fft transform workspace

  integer           ::  mtxd                                             !  dimension of the hamiltonian
  integer           ::  neig                                             !  number of eigenvectors required (maybe modified on output)
  integer           ::  kmscr(7)                                         !  max value of kgv(i,n) used for the potential fft mesh
  integer           ::  idshift                                          !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)      ::  vmax, vmin                                       !  maximum and minimum values of vscr

  integer           ::  irk,nrka
  character(len=5)  ::  labelk
  integer           ::  ipr
  integer           ::  nsfft(3)
  real(REAL64)      ::  rkin(3), rkcar(3)

  real(REAL64)      ::  rk0(3)

  real(REAL64)      ::  delrk(3), rk_p(3), rk_m(3)
  real(REAL64)      ::  xrk, delta

  character(len=1)  ::  yesno, yesno_dir
  integer           ::  nocc

  integer           ::  icoor               !  coordinate system

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors
  real(REAL64)      ::  aconv(3,3)          !  conventional lattice vectors of the Niggli cell
  real(REAL64)      ::  avecnig(3,3)        !  primitive lattice vectors of the Niggli cell
  real(REAL64)      ::  bdot(3,3), vcell

  real(REAL64)      ::  adotconv(3,3)
  real(REAL64)      ::  bvecconv(3,3), tmp(3,3)

  integer           ::  jmax

  integer           ::  kdotpsize(10)                   !  suggested sizes for kdotp problem
  integer           ::  nk                              !  number of suggested values.

  integer           ::  nmodel                          !  size of the kdotp matrix

  logical           ::  lcorrect                        !  correct value


! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  ::  C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  ::  C_UM = cmplx(UM,ZERO,REAL64)
  real(REAL64), parameter     ::  PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-14_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer    ::  i, j, k, m


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
  write(6,'("  Enter initial number of bands (greater than ",i4,")")')  &
                nint(ztot/2)
  write(6,'("  for reference calculation and exploratory suggestion")')
  write(6,'("  of k.p matrix size")')
  write(6,'("  Suggested initial value ~ ",i4)') 3*nint(ztot)
  read(5,*) neig
  write(ioreplay,*) neig,'   number of initial bands'

  write(6,*)
  write(6,*) '  The reference calculation will include ',neig,' bands'
  write(6,*)

! gets the k-point, but first generates coordinate system

  call adot_to_avec_aconv(adot,avec,bvec,aconv,avecnig)
  call adot_to_bdot(adot,vcell,bdot)

  write(6,*)
  write(6,*) '   Which coordinate system do you want to use?'
  write(6,*) '   1) Primitive lattice coordinates.'
  write(6,*) '   2) Cartesian coordinates.'
  write(6,*) '   3) Conventional lattice coordinates.'
  write(6,*)
  write(6,*) '   Enter your choice (1-3).'
  write(6,*)

  read(5,*) icoor
  write(ioreplay,'(2x,i8,"   coordinate choice")') icoor

  if(icoor <1 .or. icoor > 3) then
    write(6,*) '   Wrong choice, using primitive lattice'
    icoor = 1
  endif

! prints information about the coordinate system

  if(icoor == 1) then

    write(6,*)
    write(6,*) '   using PRIMITIVE LATTICE COORDINATES '
    write(6,*)
    write(6,*) '   The primitive lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",f8.4,")")')   &
        avec(1,1),avec(1,2),avec(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
        avec(2,1),avec(2,2),avec(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
        avec(3,1),avec(3,2),avec(3,3)
     write(6,*)

  else

    if(icoor == 2) then

      write(6,*)
      write(6,*) '   using CARTESIAN COORDINATES (atomic units)'
      write(6,*)

    else

      write(6,*)
      write(6,*) '   using CONVENTIONAL LATTICE COORDINATES '
      write(6,*)

    endif

    write(6,*) '   The conventional lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",f8.4,")")')   &
         aconv(1,1),aconv(1,2),aconv(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
         aconv(2,1),aconv(2,2),aconv(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",f8.4,")")')  &
         aconv(3,1),aconv(3,2),aconv(3,3)
     write(6,*)

  endif

  write(6,*)
  write(6,*) ' enter k-point '
  write(6,*)

  read(5,*) rkin(1),rkin(2),rkin(3)
  write(ioreplay,'(5x,3f20.10,5x,"k-point")') rkin(1),rkin(2),rkin(3)

! transforms to primitive coordinates

  if(icoor == 1) then

    do j = 1,3
      rk0(j) = rkin(j)
    enddo

    do j = 1,3
      rkcar(j) = ZERO
      do k = 1,3
        rkcar(j) = rkcar(j) + bvec(j,k)*rkin(k)
      enddo
    enddo

  elseif(icoor == 2) then

    do j = 1,3
      rk0(j) = ZERO
      do k = 1,3
        rk0(j) = rk0(j) + rkin(k)*avec(k,j)
      enddo
      rk0(j) = rk0(j) / (2*PI)
    enddo

    do j = 1,3
      rkcar(j) = rkin(j)
    enddo

  elseif(icoor == 3) then

    do j = 1,3
    do k = 1,3
      adotconv(k,j) = ZERO
      do m = 1,3
        adotconv(k,j) = adotconv(k,j) + aconv(m,j)*aconv(m,k)
      enddo
    enddo
    enddo

    call adot_to_avec(adotconv,tmp,bvecconv)

    do j = 1,3
      rkcar(j) = ZERO
      do k = 1,3
        rkcar(j) = rkcar(j) + bvecconv(j,k)*rkin(k)
      enddo
    enddo

    do j = 1,3
      rk0(j) = ZERO
      do k = 1,3
        rk0(j) = rk0(j) + rkcar(k)*avec(k,j)
      enddo
      rk0(j) = rk0(j) / (2*PI)
    enddo

  endif

  write(6,*)
  write(6,*) '  coordinates of the chosen k-point:'
  write(6,*) '      lattice coord.                  cartesian coord.'
  write(6,*)
  write(6,'(4x,3f9.4,5x,3f9.4)') (rk0(j),j=1,3), (rkcar(j),j=1,3)
  write(6,*)

! Tries to get good suggestions for the k.p model

  call kdotp_suggest_size(ztot, adot, ng, kgv, rk0, kdotpsize, nk,       &
       mxdgve)

  write(6,*)
  write(6,*) '  From a free-electron model, the suggested'
  write(6,*) '  sizes for the kdopt matrix are:'
  write(6,'(10i6)') (kdotpsize(j),j=1,nk)
  write(6,*)
  write(6,*) '  Enter one of those values (first suggestions should be better)'

  read(5,*) nmodel
  write(ioreplay,'(3x,i5,5x,"k.p matrix size")') nmodel

  lcorrect = .FALSE.
  do j = 1,nk
    if(nmodel == kdotpsize(j)) then
      lcorrect = .TRUE.
      exit
    endif
  enddo

  if(.NOT. lcorrect) then
    write(6,*)
    write(6,*) '  Value is not on the list. Enter it again.'
    write(6,*) '  It will not be checked.'
    write(6,*)

    read(5,*) nmodel
    write(ioreplay,'(3x,i5,5x,"k.p matrix size again")') nmodel
  endif

  write(6,*)
  write(6,*) '  The k.p model will include ',nmodel,' bands'
  write(6,*)

! finds mxddim, mxdbnd

  call size_mtxd(emax, rk0, adot, ng, kgv, mxddim)
  mxdbnd = nmodel+1

! allocates arrays

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))


  nocc = nmodel

  call h_kb_dia_all('pw  ', emax, rk0, nmodel+1, nocc,                   &
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

  if(abs(ei(nmodel+1)-ei(nmodel)) < 0.1) then
    write(6,*)
    write(6,*) '  The gap beteen the last level included in the model'
    write(6,'("  and the first discarded level is: ",f10.3," eV")')      &
         (ei(nmodel+1)-ei(nmodel))*HARTREE
    write(6,*)
    write(6,*) '  You should consider another model size'
    write(6,*)
  endif

  write(6,*)
  write(6,*) '  Do you want to analyze the results WITH'
  write(6,*) '  spin-orbit interaction? (y/n).'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   with SO'

  if(yesno /= 'y' .and. yesno /= 'Y') then

    ipr = 1

    nrka = -1
    call print_eig(ipr, 1, labelk, nrka, rk0,                            &
        mtxd, icmplx, neig, psi,                                         &
        adot, ei, ekpsi, isort, kgv,                                     &
        mxddim, mxdbnd, mxdgve)

    allocate(h0(mxdbnd,mxdbnd))
    allocate(dh0drk(mxdbnd,mxdbnd,3))
    allocate(d2h0drk2(mxdbnd,mxdbnd,3,3))
    allocate(ei_p(mxdbnd),ei_m(mxdbnd))

    allocate(xmass(mxdbnd,3),xgrad(mxdbnd))

    nder = 2

    call kdotp_matrix(mtxd, nmodel, psi, ei, rk0, isort, nder,           &
        h0, dh0drk, d2h0drk2,                                            &
        ng, kgv,                                                         &
        ntype, natom, rat, adot,                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

  else

    write(6,*)
    write(6,*) '  Spin-orbit included as perturbation'
    write(6,*)

    allocate(ei_so(2*mxdbnd))
    allocate(psi_so(2*mxddim,2*mxdbnd))
    allocate(ekpsi_so(2*mxdbnd))

    call spin_orbit_perturb(rk0, mtxd, isort,                            &
        nmodel, psi, ei, ei_so, psi_so, .TRUE.,                          &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

    ipr = 1

    nrka = -1
    call print_eig_so(ipr, 1, labelk, nrka, rk0,                         &
        mtxd, neig, psi_so,                                              &
        adot, ei_so, ekpsi_so, isort, kgv,                               &
        mxddim, mxdbnd, mxdgve)


    allocate(hso0(2*mxdbnd,2*mxdbnd))
    allocate(dhso0drk(2*mxdbnd,2*mxdbnd,3))
    allocate(d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3))
    allocate(ei_so_p(2*mxdbnd),ei_so_m(2*mxdbnd))

    allocate(xmass(2*mxdbnd,3),xgrad(2*mxdbnd))

    nder = 2

    call kdotp_matrix_so_pert(mtxd, nmodel, psi, ei, rk0, isort, nder,   &
        hso0, dhso0drk, d2hso0drk2,                                      &
        ng, kgv,                                                         &
        ntype, natom, rat, adot,                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

    call kdotp_matrix_so_convert(nmodel, hso0, dhso0drk, d2hso0drk2,     &
         nder,                                                           &
         mxdbnd)

  endif

! loop over directions

  do i = 1,1000

    write(6,*)
    write(6,*) ' enter k-direction '
    write(6,*)

    read(5,*) rkin(1),rkin(2),rkin(3)
    write(ioreplay,'(5x,3f20.10,5x,"k-direction")') rkin(1),rkin(2),rkin(3)

    write(6,*)
    write(6,'("  The chosen direction is: ",3f12.4)') rkin(1),rkin(2),rkin(3)
    write(6,*)

!   transforms to primitive coordinates

    if(icoor == 1) then

      do j = 1,3
        delrk(j) = rkin(j)
      enddo

    elseif(icoor == 2) then

      do j = 1,3
        delrk(j) = ZERO
        do k = 1,3
          delrk(j) = delrk(j) + rkin(k)*avec(k,j)
        enddo
      enddo

    elseif(icoor == 3) then

      do j = 1,3
        rkcar(j) = ZERO
        do k = 1,3
          rkcar(j) = rkcar(j) + bvecconv(j,k)*rkin(k)
        enddo
      enddo

      do j = 1,3
        delrk(j) = ZERO
        do k = 1,3
          delrk(j) = delrk(j) + rkcar(k)*avec(k,j)
        enddo
      enddo

    endif

!   renormalizes delrk

    xrk = ZERO
    do j = 1,3
    do k = 1,3
      xrk = xrk + delrk(k)*bdot(k,j)*delrk(j)
    enddo
    enddo
    if(xrk < EPS) THEN
      write(6,*)
      write(6,*) '  vector zero not allowed, using (1,0,0)'
      write(6,*)
      delrk(1) = UM
      delrk(2) = ZERO
      delrk(3) = ZERO
    else
      do j = 1,3
        delrk(j) = delrk(j) / sqrt(xrk)
      enddo
    endif

    delta = UM/10

    do m = 1,3

      delta = delta / 10

      do j = 1,3
        rk_p(j) = rk0(j) + delta*delrk(j)
        rk_m(j) = rk0(j) - delta*delrk(j)
      enddo

      if(yesno /= 'y' .and. yesno /= 'Y') then

        call kdotp_diag_nopsi(rk_p, nmodel, ei_p,                        &
        rk0, h0, dh0drk, d2h0drk2,                                       &
        mxdbnd)

        call kdotp_diag_nopsi(rk_m, nmodel, ei_m,                        &
        rk0, h0, dh0drk, d2h0drk2,                                       &
        mxdbnd)

        do j = 1,nmodel
          xmass(j,m) = (ei_p(j)+ei_m(j)-2*ei(j)) / (delta*delta)
          xmass(j,m) = UM / xmass(j,m)
          xgrad(j) = (ei_p(j)-ei_m(j)) / (2*delta)
        enddo

      else

        call kdotp_diag_nopsi(rk_p, 2*nmodel, ei_so_p,                   &
        rk0, hso0, dhso0drk, d2hso0drk2,                                 &
        2*mxdbnd)

        call kdotp_diag_nopsi(rk_m, 2*nmodel, ei_so_m,                   &
        rk0, hso0, dhso0drk, d2hso0drk2,                                 &
        2*mxdbnd)

        do j=1,2*nmodel
          xmass(j,m) = (ei_so_p(j)+ei_so_m(j)-2*ei_so(j)) / (delta*delta)
          xmass(j,m) = UM / xmass(j,m)
          xgrad(j) = (ei_so_p(j)-ei_so_m(j)) / (2*delta)
        enddo

      endif

    enddo

    write(6,*)
    write(6,'("   using delta = ",e12.3)') delta

    write(6,*)
    write(6,*) '        energy(eV)   effective masses                gradient'
    write(6,*) '                     100*delta  10*delta   delta'
    write(6,*)
    if(yesno /= 'y' .and. yesno /= 'Y') then
      do j = 1,nmodel
        write(6,'(i5,f12.6,3(3x,f8.4),5x,f8.3)')  j, ei(j)*HARTREE,      &
                  xmass(j,1), xmass(j,2), xmass(j,3), xgrad(j)
      if(j == nint(ztot/2)) then
        write(6,*)
        write(6,*)
      endif
    enddo
    else
      do j = 1,2*nmodel
        write(6,'(i5,f12.6,3(3x,f8.4),5x,f8.3)')  j, ei_so(j)*HARTREE,   &
                  xmass(j,1), xmass(j,2), xmass(j,3), xgrad(j)
      if(j == nint(ztot)) then
        write(6,*)
        write(6,*)
      endif
      enddo
    endif
    write(6,*)

    write(6,*)
    write(6,*) '  Do you want another direction? (y/n)'
    write(6,*)

    read(5,*) yesno_dir
    write(ioreplay,*) yesno_dir,'      New direction'


    if(yesno_dir /= 'y' .and. yesno_dir /= 'Y') exit

  enddo

  if(yesno /= 'y' .and. yesno /= 'Y') then

    deallocate(h0)
    deallocate(dh0drk)
    deallocate(d2h0drk2)
    deallocate(ei_p,ei_m)

  else

    deallocate(hso0)
    deallocate(dhso0drk)
    deallocate(d2hso0drk2)


    deallocate(ei_so)
    deallocate(psi_so)
    deallocate(ekpsi_so)

  endif

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(ekpsi)

  deallocate(xmass,xgrad)

  return
end subroutine out_effective_mass
