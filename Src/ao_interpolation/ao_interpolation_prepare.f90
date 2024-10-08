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

!>  Calculates the Hamiltonian and overlap for a grid of k-points in
!>  the Brillouin zone with non-orthogonal atomic orbitals.
!>
!>  Having computed H(k) and S(k) these quantities are dft transformed
!>  to a Wigner-Seitz cell in direct space H(R), S(R).
!>  With the latter it is then possible to obtain H(k') and S(k') in
!>  an arbitrary k' point with idft.
!>
!>  \author       Carlos Loia Reis
!>  \version      5.11
!>  \date         2014. 3 October 2024.
!>  \copyright    GNU Public License v2


subroutine ao_interpolation_prepare(ioreplay, noiData,                   &
        emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,               &
        adot, ntype, natom, rat, ntrans, mtrx,                           &
        ng, kgv, phase, conj,                                            &
        ns, inds, kmax, indv, ek,                                        &
        sfact, icmplx,                                                   &
        veff,                                                            &
        nqnl,delqnl, vkb, nkb,                                           &
        norbat, nqwf, delqwf, wvfao, lorb,                               &
        mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

! Written by Carlos Loia Reis adapting previous code.
! First stable version in 4.54.1 2014.
! Modified 25 November 2015. JLM
! Modified, mxdlao, December 1 2015.  JLM.
! Calculation and Interpolation of Hamiltonian derivatives and
! spectral weights in unfolding. November 2018. CLR.
! Fixed a major bug in IrredBz and
! cleaned up unfolding routine, October 2019. CLR.
! Modified Wigner-Seitz cell routine, October 2019. CLR.
! Modified May 2020, for 4.96. CLR
! Modified documentation May 2020, nocc, icmax, June 12 2020.jlm
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, indentation, more comments. 2 October 2024. JLM


  use NonOrthoInterp
  use IrredBZM

  implicit none

  integer, parameter            :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

  type(noiData_t)                    ::  noiData                         !<  see NonOrthoInterp

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension of orbital per atom type

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4)                   ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  integer, intent(in)                ::  iguess                          !<  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)            !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

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
  real(REAL64), intent(in)  ::   vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)          !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. normalized to vcell, hartree
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  (1/q**l) * wavefunction for atom k, ang. mom. l (unnormalized to vcell)

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  component j of eigenvector i (guess on input)
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh

  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)

  real(REAL64), allocatable          ::  ev_pw(:,:), ev_pw_irred(:,:), ev_interp(:)
  complex(REAL64), allocatable       ::  Hao(:,:), S(:,:), dh0drk_AO(:,:,:)

  real(REAL64), allocatable:: rk_grid(:,:)

! local variables

  integer                            ::  mxdorb                          !  array dimension for orbitals
  integer                            ::  mxdscr                          !  array dimension for screening potential

  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace

  integer                            ::  mtxd                            !  dimension of the hamiltonian
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point
  integer                            ::  kmscr(7)                        !  max value of kgv(i,n) used for the potential fft mesh
  integer                            ::  idshift                         !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)                       ::  vmax, vmin                      !  maximum and minimum values of vscr

  integer                            ::  irk
  integer                            ::  ipr,nd
  integer                            ::  nsfft(3)
  integer                            :: iband, jband, idir

  integer                            ::  ifail                           !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions


!--------------------------------------------------------------------
  integer           ::  nk1,nk2,nk3
  integer           ::  ws_n1, ws_n2, ws_n3

  integer           ::  nkpt,norb, nequal, lmax
  integer           ::  nocc

  type(IrredBZ_t)   ::  IrredBZ

  integer           ::  lso
  integer           ::  loptical
  character(len=1)  ::  yesno

! parameters

  real(REAL64), parameter :: UM = 1.0_REAL64

! counters

  integer           ::  i, j, k
  integer           ::  ix,iy,iz

!--------------------------------------------------------------------
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

  call size_fft(kmscr,nsfft,mxdscr,mxdwrk)

  allocate(vscr(mxdscr))

  ipr = 1

  idshift = 0
  call pot_local(ipr, vscr, vmax, vmin, veff, kmscr, idshift,            &
      ng, kgv, phase, conj, ns, inds,                                    &
      mxdscr, mxdgve, mxdnst)

  ipr = 2


!-------------------------------------------------------------------------
! initializes interpolation routines

  write(6,*)
  write(6,*) '  initializing ao interpolation'
  write(6,*)

  write(6,*) '  please enter grid of k-points'
  write(6,*) '  recomended is 4 4 4 minimum'

  read(5,*) nk1,nk2,nk3    ! this is to be removed
  write(ioreplay,*) nk1,nk2,nk3,'     k-point grid'

  write(6,'( "  interpolation in a mesh of", i3," x ",i3," x ",i3, " k-points"  )') nk1,nk2,nk3

  write(6,*) '  do you want to modify Wigner-Seitz cell search size ?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'     ws_search_size_yesno'

  if(yesno == 'y' .or. yesno == 'Y') then
     write(6,*) '  please enter Wigner Seitz cell search size'
     write(6,*) '  recomended is 2 2 2 minimum'

     read(5,*) ws_n1,ws_n2,ws_n3
     write(ioreplay,*) ws_n1,ws_n2,ws_n3,'     ws_search_size'
  else
   ws_n1 = 2
   ws_n2 = 2
   ws_n3 = 2
  endif

! Allocates the arrays on "this" and fills the "this%grid" with a call to int_pnt
! It is uniform nk1 nk2 nk3 grid with zer o shift.

! Maybe should be inlined to simplify structure

  call IrredBZInit(IrredBZ, adot, ntrans, mtrx, nk1, nk2, nk3, 0.0D0, 0.0D0, 0.0D0)

!----------------------------------------------------

  write(6,*)
  write(6,*) '  do you want a calculation including Spin-Orbit (y/n) ?'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     spin-orbit'

  if(yesno == 'y' .or. yesno == 'Y') then
    lso = 1
  else
    lso = 0
  endif

  if(lso==1) then
    write(6,*)
    write(6,*) '  calculation will be done using SO perturbation.'
    write(6,*)
  endif

!----------------------------------------------------

  write(6,*)
  write(6,*) '  do you want a calculation including Hamiltonian Derivatives (y/n) ?'
  write(6,*) '  this will be needed to compute the dielectric function '
  write(6,*) '  using interpolation (it is memory intensive...).'
  write(6,*)
  read(5,*) yesno
  write(ioreplay,*) yesno,'     optical'

  if(yesno == 'y' .or. yesno == 'Y') then
   loptical = 1
  else
   loptical = 0
  endif
  if(loptical==1) then
    write(6,*)
    write(6,*) '  calculation will also include dh0drk computation.'
    write(6,*)
  endif

  allocate(rk_grid(3,nk1*nk2*nk3))

   irk=1
     do ix=0,nk1-1
       do iy=0,nk2-1
         do iz=0,nk3-1
         rk_grid(1,irk) = UM*ix/nk1
         rk_grid(2,irk) = UM*iy/nk2
         rk_grid(3,irk) = UM*iz/nk3
!dbg     write(*,'(f8.5,f8.5,f8.5)') rk_grid(1,irk),rk_grid(2,irk),rk_grid(3,irk)
         irk=irk+1
       enddo
     enddo
   enddo

!--------------------------------------------------------------------

! finds mxddim
  mxddim = 1
  do irk=1,nk1*nk2*nk3
    rkpt(1) = rk_grid(1,irk)
    rkpt(2) = rk_grid(2,irk)
    rkpt(3) = rk_grid(3,irk)
    call size_mtxd(emax,rkpt,adot,ng,kgv,nd)
    if(nd > mxddim) mxddim = nd
  enddo

!--------------------------------------------------------------------

! finds mxdorb
  mxdorb = 0
  lmax = 0
  do k=1,ntype
    do j=1,norbat(k)
      mxdorb = mxdorb + (2*lorb(j,k)+1)*natom(k)
      lmax = max(lorb(j,k),lmax)
    enddo
  enddo

! Figure out a better heuristic...

  NEQUAL = int(ztot)+4

  MXDORB = MAX(MXDORB,NEQUAL)

! allocates arrays


  allocate(ei(mxdorb))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdorb))
  allocate(hpsi(mxddim,mxdorb))
  allocate(ekpsi(mxdorb))

  allocate(ei_so(2*mxdorb))
  allocate(psi_so(2*mxddim,2*mxdorb))


  nkpt = nk1*nk2*nk3
  norb   = mxdorb

  write(6,*)
  write(6,*) '  number k-points:',nkpt
  write(6,*) '  number of atomic orbitals:',mxdorb
  write(6,*) '  number of eigenvalues equal to full precision:',nequal
  write(6,*) '  total charge:',ztot,nint(ztot)
  write(6,*)

!--------------------------------------------------------------------

! Make a single call with nspin = 1,2...

  if (lso==1) then

    call NonOrthoInterpInit(noiData, lso, loptical, nk1,nk2,nk3,         &
        ws_n1,ws_n2,ws_n3, adot, nkpt, 2*norb, 2*nequal)
    allocate(ev_pw(2*norb,nkpt))
    allocate(ev_pw_irred(2*norb,nkpt))
    allocate(ev_interp(2*norb))
    allocate(Hao(2*norb,2*norb))
    allocate(S(2*norb,2*norb))
    allocate(dh0drk_AO(2*norb,2*norb,3))

  else

    call NonOrthoInterpInit(noiData, lso, loptical, nk1,nk2,nk3,         &
        ws_n1,ws_n2,ws_n3, adot, nkpt, norb, nequal)
    allocate(ev_pw(norb,nkpt))
    allocate(ev_pw_irred(norb,nkpt))
    allocate(ev_interp(norb))
    allocate(Hao(norb,norb))
    allocate(S(norb,norb))
    allocate(dh0drk_AO(norb,norb,3))

  endif

!-------------------------------------------------------------------------

! The call to h_kb_dia_all is the same. Invert order of if/then/else and do loop...

  if (lso==1) then

    write(6,*)
    write(6,*) '  Computing full precision eigenvalues ev(i) including Spin-Orbit'
    write(6,*)

    do irk=1, IrredBZ % nrk

      rkpt(1) = IrredBZ %rk(1,irk)
      rkpt(2) = IrredBZ %rk(2,irk)
      rkpt(3) = IrredBZ %rk(3,irk)

      write(6,'(2i5,3f8.5)') IrredBZ % nrk, irk, rkpt(1),rkpt(2),rkpt(3)

      lkpg = .FALSE.
      ipr = 0

      nocc = nequal

      call h_kb_dia_all('pw  ', emax, rkpt, nequal, nocc, flgpsd,        &
          ipr, ifail, icmax, iguess, epspsi,                             &
          ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                &
          sfact, veff, icmplx,                                           &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                       &
          psi, hpsi, ei,                                                 &
          vscr, kmscr,                                                   &
          .TRUE., norbat, nqwf, delqwf, wvfao, lorb,                     &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,        &
          mxdorb, mxdscr, mxdlao)

      lpsiso = .false.

      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
          nequal, psi, ei, ei_so, psi_so, lpsiso,                        &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve)


      do i=1,2*nequal
        ev_pw_irred(i,irk) = ei_so(i)
      enddo

    enddo

!   now unpack the eigenvalues to

!   Maybe should be inlined to simplify structure

    call UnpackEv(IrredBZ, ev_pw_irred, ev_pw, 2*norb, nkpt)

  else

    write(6,*)
    write(6,*) '  Computing full precision eigenvalues ev(i)'
    write(6,*)

    do irk=1, IrredBZ % nrk

      rkpt(1) = IrredBZ %rk(1,irk)
      rkpt(2) = IrredBZ %rk(2,irk)
      rkpt(3) = IrredBZ %rk(3,irk)

      write(*,'(2i5,3f8.5)') IrredBZ % nrk, irk, rkpt(1),rkpt(2),rkpt(3)

      lkpg = .FALSE.
      ipr = 0

      nocc = nequal

      call h_kb_dia_all('pw  ', emax, rkpt, nequal, nocc, flgpsd,        &
          ipr, ifail, icmax, iguess, epspsi,                             &
          ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                &
          sfact, veff, icmplx,                                           &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                       &
          psi, hpsi, ei,                                                 &
          vscr, kmscr,                                                   &
          .TRUE., norbat, nqwf, delqwf, wvfao, lorb,                     &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,        &
          mxdorb, mxdscr, mxdlao)

      do i=1,nequal
        ev_pw_irred(i,irk) = ei(i)
      enddo

    enddo

!   now unpack the eigenvalues to

!   Maybe should be inlined to simplify structure

    call UnpackEv(IrredBZ, ev_pw_irred, ev_pw, norb, nkpt)

  endif


!-------------------------------------------------------------------------

! Now let's get to business...

  if (lso==1) then

    write(6,*)
    write(6,*) '  Computing atomic orbitals: |psi> and overlaps: S=<psi|psi>, H=<psi|H|psi>'
    write(6,*) '  SO version'
    write(6,*)

    do irk=1, nkpt

      rkpt(1) = rk_grid(1,irk)
      rkpt(2) = rk_grid(2,irk)
      rkpt(3) = rk_grid(3,irk)

      write(*,'(2i5,3f8.5)') nkpt,irk, rkpt(1),rkpt(2),rkpt(3)

      call ao_h_and_s_spin_orbit(emax, rkpt, norb,                       &
          flgpsd,                                                        &
          ng, kgv,                                                       &
          veff, icmplx,                                                  &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mtxd, hdiag, isort, qmod, ekpg,                                &
          norbat, nqwf, delqwf, wvfao, lorb,                             &
          psi, hpsi,                                                     &
          Hao, S, dh0drk_AO,                                             &
          vscr, kmscr,                                                   &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxddim, mxdorb,        &
          mxdscr, mxdorb, mxdlao)

      call NonOrthoInterpSetGridData(noiData, irk, Hao, S, ev_pw(:,irk))

      if (loptical==1) then

        do idir=1, 3
          do iband=1, 2*norb
          do jband=1, 2*norb
            noiData%fiData%dh0drk_GridK(iband,jband,idir,irk) = dh0drk_AO(iband,jband,idir)
          enddo
          enddo
        enddo

      endif

    enddo

  else

    write(6,*)
    write(6,*) '  Computing atomic orbitals: |psi> and overlaps: S=<psi|psi>, H=<psi|H|psi>'
    write(6,*)

    do irk=1, nkpt

      rkpt(1) = rk_grid(1,irk)
      rkpt(2) = rk_grid(2,irk)
      rkpt(3) = rk_grid(3,irk)

      write(*,'(2i5,3f8.5)') nkpt,irk, rkpt(1),rkpt(2),rkpt(3)

!     hdiag etc... are not needed here.  Should simplify call

      call ao_h_and_s(emax, rkpt, norb, flgpsd,                           &
          ng, kgv,                                                        &
          veff, icmplx,                                                   &
          nqnl, delqnl, vkb, nkb,                                         &
          ntype, natom, rat, adot,                                        &
          mtxd, hdiag, isort, qmod, ekpg,                                 &
          norbat, nqwf, delqwf, wvfao, lorb,                              &
          psi, hpsi,                                                      &
          Hao, S, dh0drk_AO,                                              &
          vscr, kmscr,                                                    &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxddim, mxdorb,         &
          mxdscr, mxdorb, mxdlao)

      call NonOrthoInterpSetGridData(noiData, irk, Hao, S, ev_pw(:,irk))

      if (loptical==1) then

        do idir=1, 3
          do iband=1, norb
          do jband=1, norb
            noiData%fiData%dh0drk_GridK(iband,jband,idir,irk) = dh0drk_AO(iband,jband,idir)
          enddo
          enddo
        enddo

      endif

    enddo

  endif

!-------------------------------------------------------------------------

  write(6,*)
  write(6,*) '  Finishing Setup'
  call NonOrthoInterpFinishSetup(noiData)
  write(6,*) '  Writing to file for later use.'
  call NonOrthoInterpWriteToFile(noiData)
  write(6,*)

  deallocate(vscr)
  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(ekpsi)

  deallocate(ei_so)
  deallocate(psi_so)

  return

end subroutine ao_interpolation_prepare
