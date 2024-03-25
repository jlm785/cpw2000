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
!>  using a finite difference interpolation.
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         8 November 2023. 25 March 2024
!>  \copyright    GNU Public License v2

subroutine out_mass_fd(ioreplay,                                         &
    emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,                   &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)


! Adapted from old out_effective_mass and out_psi_test.  8 November 2023. JLM
! imethod, 25 March 2025. JLM

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

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

  integer, allocatable               ::  levdeg(:)                       !  degeneracy of energy level
  integer, allocatable               ::  leveigs(:,:)                    !  states belonging to level

  real(REAL64), allocatable          ::  deidk_fd(:)                     !  derivative of energy
  real(REAL64), allocatable          ::  d2eidk2_fd(:)                   !  second derivative of energy

  real(REAL64), allocatable          ::  ei_so(:)                        !  energy with so
  real(REAL64), allocatable          ::  deidk_fd_so(:)                  !  derivative of energy with so
  real(REAL64), allocatable          ::  d2eidk2_fd_so(:)                !  second derivative of energy with so

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

  real(REAL64)      ::  rk0(3)
  character(len=20) ::  typeofk

  real(REAL64)      ::  xk(3)
  real(REAL64)      ::  xrk, delta          !  delta is interpolation step,
  integer           ::  npt                 !  2*npt+1 is order of interpolatiopn

  character(len=1)  ::  yesno_dir
  character(len=1)  ::  yesno_so

  real(REAL64)      ::  bdot(3,3), vcell

  integer           ::  nlevel              !  number of energy levels (not counting degeneracies)
  integer           ::  maxdeg              !  maximum number of degeneracies

  integer           ::  neigin              !  initial value of neig

  logical           ::  lsoinfo             !  pseudopotential includes spin-orbit components
  logical           ::  lso                 !  calculates with spin-orbit in perturbation
  integer           ::  imethod             !  method for treating spin-orbit

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64 , UM = 1.0_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-14_REAL64
  real(REAL64), parameter     ::  TOL = 1.0E-8_REAL64
  real(REAL64), parameter     ::  HARTREE = 27.21138386_REAL64

! counters

  integer    ::  i, j, k, l, n


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

! checks if there is spin-orbit information

  lsoinfo = .FALSE.
  do n = 1,ntype
    do l = 0,3
      if(nkb(n,-1,n) /=0 .or. nkb(n,-1,n) /=0) then
       lsoinfo = .TRUE.
       exit
      endif
    enddo
    if(lsoinfo) exit
  enddo

! finds number of bands

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

  call cpw_pp_get_k_vector(rk0, rkcar, adot, typeofk, ioreplay)

  write(6,*)
  write(6,*) '  coordinates of the chosen k-point:'
  write(6,*) '      lattice coord.                  cartesian coord.'
  write(6,*)
  write(6,'(4x,3f9.4,5x,3f9.4)') (rk0(j),j=1,3), (rkcar(j),j=1,3)
  write(6,*)

! Tries to get good suggestions for the k.p model


  call size_mtxd(emax, rk0, adot, ng, kgv, mxddim)

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
  allocate(ekpsi(mxdbnd))

  call h_kb_dia_all('pw  ', emax, rk0, neigin, neigin,                   &
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

! energy levels.  first finds dimensions

  allocate(levdeg(1))
  allocate(leveigs(1,1))

  call berry_degeneracy(.TRUE., neigin, neig, ei, TOL,                   &
         nlevel, maxdeg, levdeg, leveigs,                                &
         mxdbnd, 1, 1)

  mxdlev = nlevel
  mxddeg = maxdeg

  deallocate(levdeg)
  deallocate(leveigs)

  allocate(levdeg(mxdlev))
  allocate(leveigs(mxdlev,mxddeg))

! fills the information

  call berry_degeneracy(.FALSE., neigin, neig, ei, TOL,                  &
         nlevel, maxdeg, levdeg, leveigs,                                &
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

  allocate(deidk_fd(mxdbnd))
  allocate(d2eidk2_fd(mxdbnd))

  allocate(ei_so(2*mxdbnd))
  allocate(deidk_fd_so(2*mxdbnd))
  allocate(d2eidk2_fd_so(2*mxdbnd))

! loop over directions

  do i = 1,1000

    typeofk = 'direction in k-space'

    call cpw_pp_get_k_vector(xk, rkcar, adot, typeofk, ioreplay)

    write(6,*)
    write(6,*) '  coordinates of the chosen k-direction:'
    write(6,*) '      lattice coord.                  cartesian coord.'
    write(6,*)
    write(6,'(4x,3f9.4,5x,3f9.4)') (xk(j),j=1,3), (rkcar(j),j=1,3)
    write(6,*)

!   renormalizes xk

    xrk = ZERO
    do j = 1,3
    do k = 1,3
      xrk = xrk + xk(k)*bdot(k,j)*xk(j)
    enddo
    enddo
    if(xrk < EPS) THEN
      write(6,*)
      write(6,*) '  vector zero not allowed, using (1,0,0)'
      write(6,*)
      xk(1) = UM
      xk(2) = ZERO
      xk(3) = ZERO
    endif

    write(6,*)
    write(6,*) '  choose step for and order (2 n + 1) of finite differences:'
    write(6,*)
    write(6,*) '  enter delta and n (suggested 0.0001 and 3)'
    write(6,*)
    read(5,*) delta, npt
    write(ioreplay,'(g16.8,5x,i5,10x,"delta,npt")') delta, npt

    if(delta <= ZERO) then
      write(6,*)
      write(6,*) '   negative values of delta not allowed'
      write(6,*) '   setting delta to 0.0001'
      write(6,*)
      delta = 0.0001
    endif

    if(delta < TOL) then
      write(6,*)
      write(6,*) '   unreasonably small value of delta'
      write(6,*) '   expect disaster'
      write(6,*)
    endif

    if(npt < 1 .or. npt > 10) then
      write(6,*)
      write(6,*) '   unreasonably value of n'
      write(6,*) '   setting n to 3'
      write(6,*)
      npt = 3
    endif

    write(6,*)
    write(6,'("   using delta = ",e12.3,"  and order ",i5)') delta, 2*npt+1
    write(6,*)

    lso = .FALSE.
    if(lsoinfo) then
      write(6,*)
      write(6,*) '  Do you want the results with spin-orbit in perturbation (y/n)?'
      write(6,*)
      read(5,*) yesno_so
      write(ioreplay,*) yesno_so,'      with spin-orbit'
      if(yesno_so == 'y' .or. yesno_so == 'Y') lso = .TRUE.
      if(lso) then

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
          write(6,*) '   Wrong value, using iterative diagonalization'
          write(6,*)
          imethod = 2
        endif
      endif
    endif

    call out_mass_fd_xk(rk0, xk, neig, npt, delta, lso, imethod,         &
        deidk_fd, d2eidk2_fd,                                            &
        ei_so, deidk_fd_so, d2eidk2_fd_so,                               &
        emax, flgdal, flgpsd, epspsi, icmax,                             &
        adot, ntype, natom, rat,                                         &
        ng, kgv, phase, conj,                                            &
        ns, inds, kmax, indv, ek,                                        &
        sfact, icmplx,                                                   &
        veff,                                                            &
        nqnl, delqnl, vkb, nkb,                                          &
        latorb, norbat, nqwf, delqwf, wvfao, lorb,                       &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao,          &
        mxddim, mxdbnd)


    write(6,*)
    write(6,*) '   Results from finite differences'
    if(lso) write(6,*) '   without spin-orbit'
    write(6,*)

    call out_mass_print(neig, ei, deidk_fd, d2eidk2_fd,                  &
        mxdbnd)

    if(lso) then

      call out_mass_print(2*neig, ei_so, deidk_fd_so, d2eidk2_fd_so,     &
           2*mxdbnd)


    endif

    write(6,*)
    write(6,*) '  Do you want another direction? (y/n)'
    write(6,*)

    read(5,*) yesno_dir
    write(ioreplay,*) yesno_dir,'      New direction'


    if(yesno_dir /= 'y' .and. yesno_dir /= 'Y') exit

  enddo

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(hpsi)
  deallocate(ekpsi)

  deallocate(levdeg)
  deallocate(leveigs)

  deallocate(deidk_fd)
  deallocate(d2eidk2_fd)

  deallocate(ei_so)
  deallocate(deidk_fd_so)
  deallocate(d2eidk2_fd_so)

  return

end subroutine out_mass_fd
