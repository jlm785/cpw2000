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

!> Provides orbital information to post-processing companion program
!>
!>  \author       Carlos, Loia Reis, Jose Luis Martins
!>  \version      5.11
!>  \date         8 may 2004, 26 July 2024.
!>  \copyright    GNU Public License v2

subroutine out_band_atom_info_fold(diag_type, lworkers,                  &
      pwline, title, subtitle,                                           &
      emax, flgdal, flgpsd,  epspsi, icmax, ztot, efermi,                &
      adot, ntype, natom, nameat, rat,                                   &
      ng, kgv, phase, conj,                                              &
      ns, inds, kmax, indv, ek,                                          &
      sfact, icmplx,                                                     &
      veff,                                                              &
      nqnl, delqnl, vkb, nkb,                                            &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

!  Adapted 2 June 2019. JLM.
!  Modified for unfolding June 2019. CLR.
!  Spin-Orbit and Lowdin Symmetric Orthogonalization, October 2019. CLR.
!  several bug fixes February 2020. CLR.
!  Modified, documentation, 29 May 2020. JLM
!  Modified, latorb, 7 June 2020. JLM
!  Modified, vmax, vmin, 27 November 2020. JLM
!  Modified for QtBandViewer, July 2021. CLR.
!  Modified, efermi, 29 November 2021. JLM
!  Modified annoying warning pkn, iguess. 10 November 2023. JLM
!  Modified, ztot in out_band_circuit_size. 26 July 2024. JLM
!  Modified, ao_int_, improved indentation. 8 October 2024. JLM

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

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'
  logical, intent(in)                ::  lworkers                        !<  use lworkers

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots
  character(len=60), intent(in)      ::  pwline                          !<  identifier of the calculation.  May contain miscellaneous information!


  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntype                           !<  number of types of atoms
  integer, intent(in)                ::  natom(mxdtyp)                   !<  number of atoms of type i
  character(len=2)                   ::  nameat(mxdtyp)                  !<  chemical symbol for the type i


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
  real(REAL64), intent(in)  ::   vkb(-2:mxdlqp,0:3,-1:1,mxdtyp)          !<  (1/q**l) * KB nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)      ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)       !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k

  ! input and output


! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !  number of lines in reciprocal space
  integer, allocatable               ::  nkstep(:)                       !  number of steps in line
  logical, allocatable               ::  ljump(:)                        !  indicates if the new line contains a jump from the preceeding
  integer                            ::  nvert                           !  number of vertical lines in plot
  real(REAL64), allocatable          ::  xcvert(:)                       !  x coordinate of vertical line
  real(REAL64), allocatable          ::  xk(:)                           !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk(:,:)                         !  x coordinate of k-point in plot
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot
  character(len=6), allocatable      ::  label(:)                        !  label of symmetry k-points
  real(REAL64), allocatable          ::  xklab(:)                        !  x coordinate of label

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  ev(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  integer, allocatable               ::  isort_so(:)                     !  g-vector associated with row/column i of spin-orbit hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  |psi> component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi>
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i (guess on input)

! variables for local orbitals

  integer                            ::  mxdorb                          !  array dimension of number of local orbitals
  integer                            ::  nbaslcao                        !  number of atomic orbitals
  complex(REAL64), allocatable       ::  baslcao(:,:)                    !  atomic orbitals in plane-wave basis
  complex(REAL64), allocatable       ::  baslcao_aux(:,:)                !  atomic orbitals in plane-wave basis
  integer, allocatable               ::  infolcao(:,:)                   !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  integer, allocatable               ::  infolcao_aux(:,:)               !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  real(REAL64), allocatable          ::  basxpsi(:,:,:)                  !  |<bas|psi>|^2 for each k
  complex(REAL64), allocatable       ::  prod(:,:)                       !  <bas|psi>


  integer, allocatable               ::  infolcao_so(:,:)                !  information about the original atomic orbital.  (type of atom, atom of that type, n,l,m)
  complex(REAL64), allocatable       ::  prod_so(:,:)                    !  <bas|psi>
  real(REAL64), allocatable          ::  basxpsi_so(:,:,:)               !  |<bas|psi>|^2 for each k
  complex(REAL64), allocatable       ::  psi_in(:,:)                     !  component j of eigenvector i (guess on input)

! local variables

  integer                            ::  mxdscr                          !  array dimension for screening potential

  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdbnd                          !  array dimension for the number of bands
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace

  integer                            ::  iguess                          !  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0

  integer                            ::  mtxd                            !  dimension of the hamiltonian
  integer                            ::  neig                            !  number of eigenvectors required (maybe modified on output)
  real(REAL64)                       ::  rkpt(3)                         !  j-th component in lattice coordinates of the k-point
  integer                            ::  kmscr(7)                        !  max value of kgv(i,n) used for the potential fft mesh
  integer                            ::  idshift                         !  shift of the fft mesh, used /= 0 only in highly banked memory.

  real(REAL64)                       ::  vmax, vmin                      !  maximum and minimum values of vscr

  real(REAL64)                       ::  eref                            !  reference energy for plot
  integer                            ::  nocc                            !  number of occupied states (different color) or recycled
  integer                            ::  nstyle                          !  choice of plot style

  integer                            ::  irk,nrka
  integer                            ::  iotape
  character(len=5)                   ::  labelk
  integer                            ::  ipr,nrk2,nd
  integer                            ::  nsfft(3)
  integer                            ::  nextline                        !  indicates next line on band circuit
  integer                            ::  nl                              !  keeps track of next lines
  real(REAL64)                       ::  sq

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  :: C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  i, j, n

! unfolding

  integer                            ::  iMinv(3,3)
  integer                            ::  idet
  real(REAL64)                       ::  avec(3,3)
  real(REAL64)                       ::  bvec(3,3)
  real(REAL64)                       ::  adot_pc(3,3)

  real(REAL64),allocatable           ::  pkn(:,:)
  real(REAL64),allocatable           ::  pkn_so(:,:)
  real(REAL64),allocatable           ::  pkn_tmp(:)
  real(REAL64),allocatable           ::  pkn_tmp_so(:)

  integer                            ::  idum(3,3)
  character(len=6)                   ::  fdum
  integer                            ::  ioerr
  character(len=60)                  ::  pwlinloc                        !  local version of pwline

  integer                            ::  io62                            !  tape numbers
  logical                            ::  lex

  real(REAL64), allocatable          ::  rk_fld(:,:)                     !  x coordinate of k-point in plot

  real(REAL64)                       ::  veffr1

  integer, allocatable               ::  irow(:)
  complex(REAL64), allocatable       ::  hxvec(:)

! orbital information with Lowdin orthogonalization

  complex(REAL64), allocatable       ::  S(:,:)
  complex(REAL64), allocatable       ::  S12(:,:)
  complex(REAL64), allocatable       ::  S12_inv(:,:)
  complex(REAL64), allocatable       ::  Swrk(:,:)
  real(REAL64),    allocatable       ::  ev_wrk(:)

! workers and restart

  integer                            :: irec_err
  integer                            :: irk_rd, irk_start
  integer                            :: iworker, nworker, cworker
  integer                            :: ir_size, pp_flag

  real(REAL64) :: t1, t2
  logical                            :: lmyjob

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions
  integer                            ::  ifail                           !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.

!-----------------------------------------------------------------------

  io62 = 62

!  calculates local potential in fft mesh

  if(flgdal == 'DUAL') then
    kmscr(1) = kmax(1)/2 + 4
    kmscr(2) = kmax(2)/2 + 4
    kmscr(3) = kmax(3)/2 + 4
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


!-----------------------------------------------------------------------
!!               Reads unfold information from pwline
!-----------------------------------------------------------------------
! pwline from old PW_RHO_V.DAT does not work due to a bug in out_rho_v
! this is a workaround

  fdum = '      '
  read(pwline,'(3(3i4,2x),2x,a6)',IOSTAT=ioerr) ((idum(i,j),i=1,3),j=1,3),fdum

  if(ioerr == 0 .and. adjustl(trim(fdum)) == 'fcc SL') then
    pwlinloc = pwline
  else
    write(6,*)
    write(6,*) '  Trying to use PW.DAT for the file identifier'
    write(6,*)
    open(unit=11,file='PW.DAT',status='OLD',form='FORMATTED',IOSTAT=ioerr)
    if(ioerr == 0) then
      read(11,'(a60)',IOSTAT=ioerr) pwlinloc
      close(unit=11)
    else
      write(6,*)
      write(6,*)  '  Not using unfolding'
      write(6,*)
      do i = 1,60
        pwlinloc(i:i) = ' '
      enddo
      write(pwlinloc,'(3(3i4,2x),2x,a6)') 1,0,0,  0,1,0,  0,0,1,  'fcc SL'
    endif
  endif

!!!!        PrepFold needs avec
  call adot_to_avec_sym(adot,avec,bvec)
!!!!        ploting routines need metric of primitive cell !
  call Fold_Get_adot_pc(pwlinloc, avec, adot_pc)
!-----------------------------------------------------------------------

  iotape = 13
  call out_band_circuit_size('BAND_LINES.DAT', iotape, 1, adot_pc, ztot, &
      neig, nrk2, nlines, nvert)

  allocate(xk(nrk2))
  allocate(rk(3,nrk2))
  allocate(rk_fld(3,nrk2))
  allocate(xcvert(nvert))
  allocate(ljump(nlines))
  allocate(nkstep(nlines))
  allocate(label(nvert+nlines))
  allocate(xklab(nvert+nlines))

  call out_band_get_circuit('BAND_LINES.DAT', iotape, 1, adot_pc,        &
      xk, rk, xcvert, ljump, nkstep, label, xklab,                       &
      neig, nrk2, nlines, nvert)

  allocate(e_of_k(neig,nrk2))
  allocate(e_of_k_so(2*neig,nrk2))

!-----------------------------------------------------------------------
  call Fold_Prep(pwlinloc, avec,rk, iMinv,idet, rk_fld, nrk2)
!-----------------------------------------------------------------------

! finds mxddim, mxdbnd

  mxdbnd = neig
  mxddim = 1

  irk = 1
  do j=1,3
     rkpt(j) = rk_fld(j,irk)  ! computed in folded k-point
  enddo

  call size_mtxd(emax, rkpt, adot, ng, kgv, nd)
  if(nd > mxddim) mxddim = int(1.05*nd)


! allocates arrays

!-----------------------------------------------------------------------
  allocate(pkn(nrk2,neig))
  allocate(pkn_so(nrk2,2*neig))
  allocate(pkn_tmp(neig))
  allocate(pkn_tmp_so(2*neig))
!-----------------------------------------------------------------------

  allocate(ei(mxdbnd))
  allocate(ei_so(2*mxdbnd))

  allocate(ev(mxdbnd))

  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))
  allocate(ekpsi_so(mxdbnd))

  allocate(isort_so(2*mxddim))
  allocate(psi_so(2*mxddim,2*mxdbnd))

  call size_nbaslcao(ntype, natom, norbat, lorb, nbaslcao,               &
      mxdtyp, mxdlao)


  mxdorb = nbaslcao
  allocate(baslcao(mxddim,mxdorb))
  allocate(baslcao_aux(mxddim,mxdorb))
  allocate(infolcao(5,mxdorb))
  allocate(infolcao_aux(5,mxdorb))

  allocate(basxpsi(mxdorb,mxdbnd,nrk2))
  allocate(prod(mxdorb,mxdbnd))

  allocate(infolcao_so(5,2*mxdorb))
  allocate(prod_so(2*mxdorb,2*mxdbnd))
  allocate(psi_in(2*mxddim,2*mxdorb))
  allocate(basxpsi_so(2*mxdorb,2*mxdbnd,nrk2))

  allocate(S(mxdorb,mxdorb))
  allocate(S12(mxdorb,mxdorb))
  allocate(S12_inv(mxdorb,mxdorb))
  allocate(Swrk(mxdorb,mxdorb))
  allocate(ev_wrk(mxdorb))

  allocate(hxvec(mxdorb))
  allocate(irow(mxdorb))

  iguess = 0
  nextline = 1
  nl = 1

  veffr1 = real(veff(1),REAL64)

  if(lworkers) then
    write(6,*) "Using Workers and Restart Capablities"

    write(6,*) 'how many workers ?'
    read(5,*) nworker

    write(6,*) 'current worker ?'
    read(5,*) iworker
  else
    nworker = 1
    iworker = 1
  endif

  irk =1
  inquire(iolength=ir_size) irk_rd,                                      &
                            e_of_k(:,irk),                               &
                            pkn_tmp(:),  e_of_k_so(:,irk),               &
                            pkn_tmp_so(:),                               &
                            basxpsi(:,:,irk),                            &
                            basxpsi_so(:,:,irk),                         &
                            infolcao,                                    &
                            infolcao_so

!   if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io62, exist=lex)
    if(lex) close(unit = io62)
    call execute_command_line("rm  band_info_rec.dat 2> /dev/null ")
  endif

  open(unit = io62, file ="band_info_rec.dat", access="direct", recl=ir_size)

  irk_start=1

  pp_flag = 0

  call zeelap(t1)

! loop over k-points

  do irk=1,nrk2

    cworker = mod(irk-1,nworker) +1
    lmyjob = .TRUE.
    if (cworker== iworker) then
      irk_rd = -10
      read(io62,rec=irk, iostat=irec_err) irk_rd,                        &
                                          e_of_k(:,irk),                 &
                                          pkn_tmp(:), e_of_k_so(:,irk),  &
                                          pkn_tmp_so(:),                 &
                                          basxpsi(:,:,irk),              &
                                          basxpsi_so(:,:,irk),           &
                                          infolcao,                      &
                                          infolcao_so

      pkn(irk,:) = pkn_tmp(:)
      pkn_so(irk,:) = pkn_tmp_so(:)


      if(irk_rd ==irk) then
        write(*,'("not computing k-point #",i5)') irk
        lmyjob = .FALSE.
      endif
    else
      write(*,'(" not computing k-point # ",i5)') irk
      lmyjob = .FALSE.
    endif

    if(lmyjob) then

      do j=1,3
        rkpt(j) = rk_fld(j,irk)
      enddo

      lkpg = .FALSE.
      ipr = 0

      nocc = neig

      call h_kb_dia_all(diag_type, emax, rkpt, neig, nocc,               &
          flgpsd, ipr, ifail, icmax, iguess, epspsi,                     &
          ng, kgv, phase, conj, ns, inds, kmax, indv, ek,                &
          sfact, veff, icmplx,                                           &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mtxd, hdiag, isort, qmod, ekpg, lkpg,                          &
          psi, hpsi, ei,                                                 &
          vscr, kmscr,                                                   &
          latorb, norbat, nqwf, delqwf, wvfao, lorb,                     &
          mxdtyp, mxdatm, mxdgve, mxdnst, mxdcub, mxdlqp, mxddim,        &
          mxdbnd, mxdscr, mxdlao)

      do j=1, neig
        ev(j) = ei(j)
      enddo

!-----------------------------------------------------------------------
     call Fold_GetPkn(pkn, iMinv, idet, irk, kgv, isort, psi, nrk2,      &
         neig, mtxd, ng, mxdgve, mxddim, mxdbnd)
!-----------------------------------------------------------------------

      call atomic_orbital_c16(rkpt, mtxd, isort, 1,                      &
          nbaslcao, baslcao_aux, infolcao,                               &
          ng, kgv,                                                       &
          norbat, nqwf, delqwf, wvfao, lorb,                             &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdorb, mxdgve, mxdlao)

      call zgemm('C', 'N', nbaslcao, nbaslcao, mtxd, C_UM, baslcao_aux,  &
          mxddim, baslcao_aux, mxddim, C_ZERO, S, mxdorb)

      call  ao_int_GetS12(S, S12, S12_inv, Swrk, ev_wrk, nbaslcao)

      call zgemm('n', 'n', mtxd, nbaslcao, nbaslcao, C_UM, baslcao_aux,  &
          mxddim, S12, mxdorb, C_ZERO, baslcao, mxddim)

      call zgemm('C', 'N', nbaslcao, nbaslcao, mtxd, C_UM, baslcao,      &
          mxddim, baslcao, mxddim, C_ZERO, S, mxdorb)

      do n = 1,nbaslcao
        sq = ZERO
        do j = 1,mtxd
          sq = sq + real(baslcao(j,n)*conjg(baslcao(j,n)),REAL64)
        enddo
        sq = UM/sqrt(sq)
        do j = 1,mtxd
          baslcao(j,n) = sq * baslcao(j,n)
        enddo
      enddo

      call zgemm('C', 'N', nbaslcao, neig, mtxd, C_UM, baslcao, mxddim,  &
          psi, mxddim, C_ZERO, prod, mxdorb)

      do n = 1,neig
      do j = 1,nbaslcao
        basxpsi(j,n,irk) = real(prod(j,n)*conjg(prod(j,n)),REAL64)
      enddo
      enddo

      lpsiso = .true.
      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
          neig, psi, ei, ei_so, psi_so, lpsiso ,                         &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      do i=1,nbaslcao
       infolcao_so(:,2*i)   = infolcao(:,i)
       infolcao_so(:,2*i-1) = infolcao(:,i)
      enddo

      do i=1,nbaslcao
      do j=1,mtxd
        psi_in(2*j-1,2*i  ) = baslcao(j,i)
        psi_in(2*j  ,2*i  ) = C_ZERO
        psi_in(2*j-1,2*i-1) = C_ZERO
        psi_in(2*j  ,2*i-1) = baslcao(j,i)
      enddo
      enddo

      call zgemm('c', 'n', 2*nbaslcao, 2*neig, 2*mtxd, C_UM, psi_in,     &
          2*mxddim, psi_so, 2*mxddim, C_ZERO, prod_so, 2*mxdorb)

      do n = 1,2*neig
      do j = 1,2*nbaslcao
        basxpsi_so(j,n,irk) = real(prod_so(j,n)*conjg(prod_so(j,n)),REAL64)
      enddo
      enddo

!-----------------------------------------------------------------------
      call Fold_GetPknSO(pkn_so,iMinv,idet,irk,kgv,isort,psi_so, nrk2,   &
          2*neig, mtxd, ng, mxdgve,2*mxddim,2*mxdbnd)
!-----------------------------------------------------------------------

      call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                  &
          mxddim,mxdbnd)

      ipr = 1
      nrka = -1
      call print_eig(ipr, irk, labelk, nrka, rkpt,                       &
          mtxd, icmplx, neig, psi,                                       &
          adot, ei, ekpsi, isort, kgv,                                   &
          mxddim, mxdbnd, mxdgve)

      if(ipr == 2) then
        call kinetic_energy_so(neig, mtxd, ekpg, psi_so, ekpsi_so,       &
            mxddim, mxdbnd)
      endif

      call print_eig_so(ipr, irk, labelk, nrka, rkpt,                    &
          mtxd, neig, psi_so,                                            &
          adot, ei_so, ekpsi_so, isort, kgv,                             &
          mxddim, mxdbnd, mxdgve)


      write(6,'( "iworker #",i5, "   writing in irk # "                  &
          &        ,i5, "   of ", i5)') iworker, irk,nrk2


      if (cworker== iworker) then

        pkn_tmp(:) = pkn(irk,:)
        pkn_tmp_so(:) = pkn_so(irk,:)

        write(io62,rec=irk) irk, ev(:),                                  &
                            pkn_tmp(:),ei_so(:),                         &
                            pkn_tmp_so(:),                               &
                            basxpsi(:,:,irk),                            &
                            basxpsi_so(:,:,irk),                         &
                            infolcao,                                    &
                            infolcao_so

      endif

    endif

  enddo

  do irk = 1, nrk2
    read(io62,rec=irk, iostat=irec_err) irk_rd,                          &
                                        e_of_k(:,irk),                   &
                                        pkn_tmp(:),  e_of_k_so(:,irk),   &
                                        pkn_tmp_so(:),                   &
                                        basxpsi(:,:,irk),                &
                                        basxpsi_so(:,:,irk),             &
                                        infolcao,                        &
                                        infolcao_so

      pkn(irk,:) = pkn_tmp(:)
      pkn_so(irk,:) = pkn_tmp_so(:)


    if(irk_rd /=irk) then
      pp_flag = 0
      exit
    endif
    pp_flag = 1
  enddo

  call zeelap(t2)
  write(*,*)
  write(*,'(" elapsed time (s):", 2f14.3)') (t2-t1)

  if (pp_flag ==1) then
    write(*,*) 'Generating band structure files!'
  else
    close(io62)
    stop "this worker is done. run again when all workers are done to see results."
  endif

  n = min(nint(0.5*ztot + 0.01),neig)
  eref = e_of_k(n,1)
  do irk = 1,nrk2
  do j=1,n
    if(e_of_k(j,irk) > eref) eref = e_of_k(j,irk)
  enddo
  enddo

  nocc = n

! writes the output files for xmgrace and gnuplot

  iotape = 15
  nstyle = 2

  call out_band_eref(neig,nrk2,ztot,efermi,2,1,e_of_k,eref,nocc)

!-----------------------------------------------------------------------
  call out_band_fold_xmgrace('band_fld_ref.agr', iotape,                 &
      title, subtitle, nstyle,                                           &
      pkn, neig, nrk2, xk, e_of_k, eref, nocc,                           &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab)
!-----------------------------------------------------------------------

  call out_band_info_write('BAND.DAT', iotape,                           &
      title, subtitle, nstyle,                                           &
      neig, nrk2, rk, rk_fld, xk, e_of_k, eref, nocc,                    &
      nbaslcao, infolcao, basxpsi, pkn,                                  &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab, ntype, nameat)

  call out_band_eref(neig, nrk2, ztot, efermi, 1, 1, e_of_k_so, eref, nocc)

  n = min(nint(ztot + 0.01),2*neig)
  eref = e_of_k_so(n,1)
  do irk = 1,nrk2
  do j=1,n
    if(e_of_k_so(j,irk) > eref) eref = e_of_k_so(j,irk)
  enddo
  enddo

  nocc = n

!-----------------------------------------------------------------------
  call out_band_fold_xmgrace('band_fld_ref_so.agr', iotape,              &
      title, subtitle, nstyle,                                           &
      pkn_so, 2*neig, nrk2, xk, e_of_k_so, eref, nocc,                   &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab)
!-----------------------------------------------------------------------

  call out_band_info_write('BAND_SO.DAT', iotape,                        &
      title, subtitle, nstyle,                                           &
      2*neig, nrk2, rk, rk_fld, xk, e_of_k_so, eref, nocc,               &
      2*nbaslcao, infolcao_so, basxpsi_so, pkn_so,                       &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab, ntype, nameat)

!------------------------------------------------------------------
  deallocate(pkn)
  deallocate(pkn_so)
  deallocate(pkn_tmp)
  deallocate(pkn_tmp_so)
!------------------------------------------------------------------

  deallocate(nkstep)
  deallocate(ljump)

  deallocate(xcvert)
  deallocate(xk)
  deallocate(rk)
  deallocate(e_of_k)
  deallocate(label)
  deallocate(xklab)

  deallocate(vscr)

  deallocate(ei)
  deallocate(hdiag)
  deallocate(isort)
  deallocate(qmod)
  deallocate(ekpg)
  deallocate(psi)
  deallocate(ekpsi)

  deallocate(psi_so)

  deallocate(baslcao)
  deallocate(infolcao)

  deallocate(basxpsi)
  deallocate(prod)

  deallocate(psi_in)
  deallocate(infolcao_so)

  deallocate(basxpsi_so)
  deallocate(prod_so)

  return

end subroutine out_band_atom_info_fold
