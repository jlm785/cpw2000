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

!>  Calculates the band structure along a path
!>  in the Brillouin zone of the parent structure.
!>
!>
!>  \author       Carlos Loia reis, Jose Luis Martins
!>  \version      5.11
!>  \date         8 May 2004. 26 July 2024.
!>  \copyright    GNU Public License v2

subroutine out_band_fold(diag_type, lworkers,                            &
      pwline, title, subtitle,                                           &
      emax, flgdal, flgpsd, epspsi, icmax, ztot,                         &
      adot, ntype, natom, rat,                                           &
      ng, kgv, phase, conj,                                              &
      ns, inds, kmax, indv, ek,                                          &
      sfact, icmplx,                                                     &
      veff,                                                              &
      nqnl, delqnl, vkb, nkb,                                            &
      latorb, norbat, nqwf, delqwf, wvfao, lorb,                         &
      mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)


! version 4.42. 8 may 2004. jlm
! modified 11 february 2008 (read from file). JLM
! modified September 5, 2012 (f90, circuit). jlm
! modified October 18, 2013 (band index permutations). jlm
! modified January 8, 2013 (new interface). jlm
! modified, dimensions vkb, March 31, 2014. jlm
! modified, unfolding version, July 9, 2014. clr
! modified title, 6 August 2014. JLM
! parallel execution and restart capabilities, October 2019, clr
! adpated for version 4.96, clr
! Modified, documentation 29 May 2020. JLM
! Modified, latorb, 7 June 2020. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, iguess, annoying pkn warning and presentation. 10 November 2023. JLM
! Modified, ztot in out_band_circuit_size. 26 July 2024. JLM
! Modified, length of labels, 24 September 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdtyp                          !<  array dimension of types of atoms
  integer, intent(in)                ::  mxdatm                          !<  array dimension of number of atoms of a given type
  integer, intent(in)                ::  mxdgve                          !<  array dimension for g-space vectors
  integer, intent(in)                ::  mxdnst                          !<  array dimension for g-space stars
  integer, intent(in)                ::  mxdlqp                          !<  array dimension for local potential
  integer, intent(in)                ::  mxdcub                          !<  array dimension for 3-index g-space
  integer, intent(in)                ::  mxdlao                          !<  array dimension for 3-index g-space

  character(len=4), intent(in)       ::  diag_type                       !<  selects diagonalization, 'pw  ','ao  ','aojc'
  logical, intent(in)                ::  lworkers                        !<  use lworkers

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots
  character(len=250), intent(in)     ::  pwline                          !<  identifier of the calculation.  May contain miscellaneous information!

  real(REAL64), intent(in)           ::  emax                            !<  largest kinetic energy included in hamiltonian diagonal. (hartree).
  character(len=4), intent(in)       ::  flgdal                          !<  dual approximation if equal to 'DUAL'
  character(len=6), intent(in)       ::  flgpsd                          !<  type of pseudopotential
  real(REAL64), intent(in)           ::  epspsi                          !<  requested precision of the eigenvectors
  integer, intent(in)                ::  icmax                           !<  maximum number of iterations for diagonalization
  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

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
  integer, intent(in)                ::  nkb(0:3,-1:1,mxdtyp)            !<   KB pseudo.  normalization for atom k, ang. mom. l

  logical, intent(in)                ::  latorb                          !<  indicates if all atoms have information about atomic orbitals
  integer, intent(in)                ::  norbat(mxdtyp)                  !<  number of atomic orbitals for atom k
  integer, intent(in)                ::  nqwf(mxdtyp)                    !<  number of points for wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  delqwf(mxdtyp)                  !<  step used in the wavefunction interpolation for atom k
  real(REAL64), intent(in)           ::  wvfao(-2:mxdlqp,mxdlao,mxdtyp)  !<  wavefunction for atom k, ang. mom. l (normalized to vcell)
  integer, intent(in)                ::  lorb(mxdlao,mxdtyp)             !<  angular momentum of orbital n of atom k


! allocatable arrays for Brillouin zone path

  integer                            ::  nlines                          !   number of lines in reciprocal space
  integer, allocatable               ::  nkstep(:)                       !   number of steps in line
  logical, allocatable               ::  ljump(:)                        !   indicates if the new line contains a jump from the preceeding
  integer                            ::  nvert                           !   number of vertical lines in plot
  real(REAL64), allocatable          ::  xcvert(:)                       !   x coordinate of vertical line
  real(REAL64), allocatable          ::  xk(:)                           !   x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk(:,:)                         !   x coordinate of k-point in plot
  real(REAL64), allocatable          ::  rk_fld(:,:)                     !   x coordinate of k-point in plot
  real(REAL64), allocatable          ::  e_of_k(:,:)                     !   band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !   spin-orbit band energies of k-point in plot
  character(len=10), allocatable     ::  label(:)                        !   label of symmetry k-points
  real(REAL64), allocatable          ::  xklab(:)                        !   x coordinate of label

! allocatable arrays with larger scope

  real(REAL64), allocatable          ::  ei(:)                           !   eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  ev(:)                           !   eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  ei_so(:)                        !   spin-orbit eigenvalue (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !   hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !   g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !   length of k+g-vector of row/column i
  real(REAL64), allocatable          ::  ekpg(:)                         !   kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !   component j of eigenvector i (guess on input)
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !   H | psi>
  real(REAL64), allocatable          ::  ekpsi(:)                        !   kinetic energy of eigenvector i. (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !   kinetic energy of eigenvector i. (hartree)

  real(REAL64), allocatable          ::  vscr(:)                         !   screened potential in the fft real space mesh
  complex(REAL64), allocatable       ::  psi_so(:,:)                     !   component j of eigenvector i (guess on input)

! local variables

  integer                            ::  mxdscr                          !  array dimension for screening potential

  logical                            ::  lkpg                            !  If true use the previous G-vectors (same mtxd and isort)
  logical                            ::  lpsiso                          !  If true calculates the spin-orbit perturbed wave-functions

  integer                            ::  iguess                          !  if guess eigenvectors are available, iguess = 1, otherwise iguess = 0
  integer                            ::  ifail                           !  if ifail=0 the ditsp_c16 was successfull. Otherwise ifail indicates the number of correct digits.


  integer                            ::  mxddim                          !  array dimension for the hamiltonian
  integer                            ::  mxdbnd                          !  array dimension for the number of bands
  integer                            ::  mxdwrk                          !  array dimension for fft transform workspace

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

  integer                            ::  io63                            !  tape numbers
  logical                            ::  lex

! workers and restart

  integer    :: irec_err
  integer    :: irk_rd, irk_start
  integer    :: iworker, nworker, cworker
  integer    :: ir_size, pp_flag

  real(REAL64) :: t1, t2
  logical                            :: lmyjob

! counters

  integer    ::  i, j, n


  io63 = 63

! calculates local potential in fft mesh
! be more generous for fold as k-points are far from Gamma point

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
  ng, kgv, phase, conj, ns, inds,                                        &
  mxdscr, mxdgve, mxdnst)

!------------------------------------------------------------------
!!   pwline from old PW_RHO_V.DAT does not work due to a bug in out_rho_v
!!   this is a workaround

! checks it is the debugged version

  fdum = '      '
  read(pwline,'(3(3i4,2x),2x,a6)',IOSTAT=ioerr)                          &
            ((idum(i,j),i=1,3),j=1,3),fdum

  if(ioerr == 0 .and. adjustl(trim(fdum)) == 'fcc SL') then

    pwlinloc = pwline(1:60)

  else

    write(6,*)
    write(6,*) '  Trying to use PW.DAT for the file identifier'
    write(6,*)

    open(unit=11, file='PW.DAT', status='OLD', form='FORMATTED', IOSTAT=ioerr)

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

!!
!!     PrepFold needs avec
  call adot_to_avec_sym(adot,avec,bvec)
!!     ploting routines need metric of primitive cell !
  call Fold_Get_adot_pc(pwlinloc, avec, adot_pc)
!------------------------------------------------------------------

  iotape = 13
  call out_band_circuit_size('BAND_LINES.DAT', iotape, 1, adot_pc, ztot, &      ! note call with adot_pc
       neig, nrk2, nlines, nvert)

  allocate(xk(nrk2))
  allocate(rk(3,nrk2))
  allocate(rk_fld(3,nrk2))
  allocate(xcvert(nvert))
  allocate(ljump(nlines))
  allocate(nkstep(nlines))
  allocate(label(nvert+nlines))
  allocate(xklab(nvert+nlines))


  call out_band_get_circuit('BAND_LINES.DAT', iotape, 1, adot_pc,        &       ! note call with adot_pc
       xk, rk, xcvert, ljump, nkstep, label, xklab,                      &
       neig, nrk2, nlines, nvert)

  allocate(e_of_k(neig,nrk2))
  allocate(e_of_k_so(2*neig,nrk2))

!------------------------------------------------------------------
  call Fold_Prep(pwlinloc, avec,rk, iMinv,idet, rk_fld, nrk2)
!------------------------------------------------------------------

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

!------------------------------------------------------------------
  allocate(pkn(nrk2,neig))
  allocate(pkn_so(nrk2,2*neig))
  allocate(pkn_tmp(neig))
  allocate(pkn_tmp_so(2*neig))
!------------------------------------------------------------------

  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))

  allocate(ekpsi(mxdbnd))
  allocate(ekpsi_so(2*mxdbnd))

  allocate(ei_so(2*mxdbnd))

  allocate(ev(mxdbnd))

  allocate(psi_so(2*mxddim,2*mxdbnd))


!-

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

  irk=1
  inquire(iolength = ir_size) irk_rd, e_of_k(:,irk),                     &
      pkn_tmp(:), e_of_k_so(:,irk),                                      &
      pkn_tmp_so(:)

! if run by a human do not restart

  if(.not. lworkers) then
    inquire(unit = io63, exist=lex)
    if(lex) close(unit = io63)
    call execute_command_line("rm  band_fold_rec.dat 2> /dev/null ")
  endif

  open(unit = io63, file ="band_fold_rec.dat", access="direct", recl=ir_size)

  irk_start=1
  pp_flag = 0

  call zeelap(t1)

! loop over k-points
  do irk=1,nrk2

    cworker = mod(irk-1,nworker) +1
    lmyjob = .TRUE.
    if(cworker== iworker) then
      irk_rd = -10
      read(io63,rec=irk, iostat=irec_err) irk_rd,                        &
      e_of_k(:,irk), pkn_tmp(:),                                         &
      e_of_k_so(:,irk),pkn_tmp_so(:)

      pkn(irk,:) = pkn_tmp(:)
      pkn_so(irk,:) = pkn_tmp_so(:)

      if(irk_rd ==irk) then
        write(6,'("not computing k-point #",i5)') irk
        lmyjob = .FALSE.
      endif
    else
      write(6,'("not computing k-point #",i5)') irk
      lmyjob = .FALSE.
    endif

    if(lmyjob) then

      do j=1,3
        rkpt(j) = rk_fld(j,irk)                                          ! computed in folded k-point
      enddo

      call size_mtxd(emax,rkpt,adot,ng,kgv,nd)

      if(nd > mxddim) then
        write(*,'("mtxd has changed from ", i5, " to", i5)')             &
            mxddim, int(nd*1.01)

        mxddim = int(nd*1.01)

        deallocate(hdiag)
        deallocate(isort)
        deallocate(qmod)
        deallocate(ekpg)
        deallocate(psi)
        deallocate(hpsi)
        deallocate(ekpsi)
        deallocate(ekpsi_so)

        deallocate(psi_so)

        allocate(hdiag(mxddim))
        allocate(isort(mxddim))
        allocate(qmod(mxddim))
        allocate(ekpg(mxddim))
        allocate(psi(mxddim,mxdbnd))
        allocate(hpsi(mxddim,mxdbnd))
        allocate(ekpsi(mxdbnd))
        allocate(ekpsi_so(2*mxdbnd))

        allocate(psi_so(2*mxddim,2*mxdbnd))

      endif

      lkpg = .FALSE.
      ipr = 0

      iguess = 0
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


      if(ifail /= 0) then
        if(ifail < 3) then
          write(6,'("   stopped in out_band_dos:  failed diagon.",       &
             &     " number of accurate digits = ",i5)') ifail
          stop
        else
          write(6,*)
          write(6,*) "  WARNING,  number of accurate digits = ",ifail
          write(6,*)
        endif
      endif

!--------------------------------------------------------------------
      call Fold_GetPkn(pkn, iMinv, idet, irk, kgv, isort, psi,           &
          nrk2, neig, mtxd, ng, mxdgve ,mxddim, mxdbnd)
!--------------------------------------------------------------------

      lpsiso = .true.
      call spin_orbit_perturb(rkpt, mtxd, isort,                         &
          neig, psi, ei, ei_so, psi_so, lpsiso,                          &
          ng, kgv,                                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

!--------------------------------------------------------------------
      call Fold_GetPknSO(pkn_so, iMinv, idet, irk, kgv, isort, psi_so,   &
          nrk2, 2*neig, mtxd, ng, mxdgve, 2*mxddim, 2*mxdbnd)
!------------------------------------------------------------------

      call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                  &
          mxddim, mxdbnd)

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
              &    ,i5, "   of ", i5)') iworker, irk,nrk2
      write(6,*)


      if (cworker== iworker) then

        pkn_tmp(:) = pkn(irk,:)
        pkn_tmp_so(:) = pkn_so(irk,:)

        write(io63,rec=irk) irk, ev(:),                                  &
            pkn_tmp(:),ei_so(:),                                         &
            pkn_tmp_so(:)

      endif

    endif

  enddo

  do irk = 1, nrk2

    read(io63,rec=irk, iostat=irec_err) irk_rd,                          &
        e_of_k(:,irk),   pkn_tmp(:),                                     &
        e_of_k_so(:,irk),pkn_tmp_so(:)

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
    close(io63)
    write(*,*) "this worker is done. run again when all workers ",       &
                "are done to see results.", iworker
    stop
  endif

! writes the output files for xmgrace

  iotape = 15
  nstyle = 2

  n = min(nint(0.5*ztot + 0.01),neig)
  eref = e_of_k(n,1)
  do irk = 1,nrk2
  do j=1,n
    if(e_of_k(j,irk) > eref) eref = e_of_k(j,irk)
  enddo
  enddo

  nocc = n
!------------------------------------------------------------------
  call out_band_fold_xmgrace('band.agr', iotape,                         &
      title, subtitle, nstyle,                                           &
      pkn, neig, nrk2, xk, e_of_k, eref, nocc,                           &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab)
!------------------------------------------------------------------

  n = min(nint(ztot + 0.01),2*neig)
  eref = e_of_k_so(n,1)
  do irk = 1,nrk2
  do j=1,n
    if(e_of_k_so(j,irk) > eref) eref = e_of_k_so(j,irk)
  enddo
  enddo

  nocc = n
!------------------------------------------------------------------
  call out_band_fold_xmgrace('band_so.agr', iotape,                      &
      title, subtitle, nstyle,                                           &
      pkn_so, 2*neig, nrk2, xk, e_of_k_so, eref,nocc,                    &
      nvert, xcvert, nlines, ljump, nkstep, label, xklab)
!------------------------------------------------------------------


!------------------------------------------------------------------
  deallocate(rk_fld)
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

  return

end subroutine out_band_fold
