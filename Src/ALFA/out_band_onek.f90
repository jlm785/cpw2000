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

!>  Calculates bands, hamiltonian components, kdotp matrices,
!>  oscillator strengths, for a given k-vector
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         8 May 2004. 4 November 2025.
!>  \copyright    GNU Public License v2

subroutine out_band_onek(ioreplay,                                       &
    emax, flgdal, flgpsd, iguess, epspsi, icmax, ztot,                   &
    adot, ntype, natom, rat,                                             &
    ng, kgv, phase, conj,                                                &
    ns, inds, kmax, indv, ek,                                            &
    sfact, icmplx,                                                       &
    veff,                                                                &
    nqnl, delqnl, vkb, nkb,                                              &
    latorb, norbat, nqwf, delqwf, wvfao, lorb,                           &
    mxdtyp, mxdatm, mxdgve, mxdnst, mxdlqp, mxdcub, mxdlao)

! version 4.42. 8 may 2004. jlm
! modified 11 february 2008 (read from file). JLM
! modified September 5, 2012 (f90, circuit). jlm
! modified October 18, 2013 (band index permutations). jlm
! modified January 8, 2013 (new interface). jlm
! modified, vkb dimensions, March 31, 2014. jlm
! modified 4.7X November 2015. JLM
! Modified, documentation, 21 February 2020. JLm
! Modified h_kb_dia_all, icmax, 7-14 June 2020. JLM
! Modified, qmod-->ekpg in hk_psi, psi_h_psi. 13 February 2021. JLM
! Modified, vmax, vmin, 27 November 2020. JLM
! Modified, name oscillator_strength. 16 May 2024. JLM
! Modified, more flexibility in oscillator strength. 14 May 2025. JLM
! Modified, input of desired k-point in cpw_pp_get_k_vector, 24 September 2025. JLM
! Modified, oscillator strength on a given direction. 25 October 2025. JLM
! Modified, correct for dpin degeneracy. 3 November 2025. JLM
! Print spin hamiltonian. 4 November 2025. JLM

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

  complex(REAL64), allocatable       ::  psi_tmp(:,:)                    !  component j of eigenvector i
  complex(REAL64), allocatable       ::  hloc(:,:)                       !  local component of reduced hamiltonian
  complex(REAL64), allocatable       ::  hkin(:,:)                       !  kinetic component of reduced hamiltonian
  complex(REAL64), allocatable       ::  hnl(:,:)                        !  non-local pseudopotential component of reduced hamiltonian

  complex(REAL64), allocatable       ::  h0(:,:)                         !  <Psi|H|Psi> without spin-orbit
  complex(REAL64), allocatable       ::  dh0drk(:,:,:)                   !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2h0drk2(:,:,:,:)               !  d^2 <Psi|H|Psi> d k^2

  complex(REAL64), allocatable       ::  hso0(:,:)                       !  <Psi|H|Psi> with spin-orbit
  complex(REAL64), allocatable       ::  dhso0drk(:,:,:)                 !  d <Psi|H|Psi> d k
  complex(REAL64), allocatable       ::  d2hso0drk2(:,:,:,:)             !  d^2 <Psi|H|Psi> d k^2

  complex(REAL64), allocatable       ::  vec_so(:,:)                     !  spin-orbit vectors
  complex(REAL64), allocatable       ::  dh_so(:,:,:)                    !  dhso0drk in the vec_so representation
  complex(REAL64), allocatable       ::  tmp1(:,:),tmp2(:,:)             !  temporary

  real(REAL64), allocatable          ::  vscr_sp(:,:)

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

  real(REAL64)      ::  rk0(3), rkcar0(3)                                !  k-point of interest
  character(len=20) ::  typeofk                                          !  label for type of k-point

  integer           ::  natot                                            !  total number of atoms
  integer           ::  neltot                                           !  total number of electrons
  real(REAL64)      ::  emidgap                                          !  rough estimate of the mid-gap

  character(len=1)  ::  yesno
  character(len=1)  ::  yesno2, yesno3
  integer           ::  nsmall

  integer           ::  info
  integer           ::  nocc

  integer           ::  ninitbeg, ninitend                               !  begin and end of initial state index for oscillator strength
  integer           ::  nfinalbeg, nfinalend                             !  begin and end of final state index

  logical           ::  lpair, lexcit, lxyz                              !  display format of oscillator strength
  character(len=20) ::  typeofr                                          !  label for type of r-point                                         !
  real(REAL64)      ::  rdir(3)                                          !  choice of r-point (lattice coordinates)
  real(REAL64)      ::  rdircar(3)                                       !  choice of r-point (cartesian coordinates)

! constants

  real(REAL64), parameter     :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  complex(REAL64), parameter  :: C_ZERO = cmplx(ZERO,ZERO,REAL64)
  complex(REAL64), parameter  :: C_UM = cmplx(UM,ZERO,REAL64)

! counters

  integer    ::  i, j, k, n


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


! finds mxddim, mxdbnd

  write(6,*)
  write(6,*) ' enter number of bands '
  read(5,*) neig
  write(ioreplay,*) neig,'   number of bands'

  mxdbnd = neig

! gets k-vector

  typeofk = 'reference k-point   '

  call cpw_pp_get_k_vector(rk0, rkcar0, adot, typeofk, ioreplay)


  call size_mtxd(emax, rk0, adot, ng, kgv, mxddim)

! allocates arrays


  allocate(ei(mxdbnd))
  allocate(hdiag(mxddim))
  allocate(isort(mxddim))
  allocate(qmod(mxddim))
  allocate(ekpg(mxddim))
  allocate(psi(mxddim,mxdbnd))
  allocate(hpsi(mxddim,mxdbnd))
  allocate(ekpsi(mxdbnd))


  nocc = neig

  call h_kb_dia_all('pw  ', emax, rk0, neig, nocc,                       &
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

! orient the wave-functions to simplify printing

  call psi_orient_xyz(rk0, adot, mtxd, neig, isort,                      &
      psi, ei,                                                           &
      ng, kgv,                                                           &
      mxddim, mxdbnd, mxdgve)

  write(6,*)
  write(6,*) '  Do you want to analyze the results WITHOUT'
  write(6,*) '  spin-orbit interaction? (y/n)'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   without SO'

! Without spin-orbit

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*)
    write(6,*) '  Do you want to see the eigenvalues (0,1,2,3)'
    write(6,*)
    write(6,*) '  0) No'
    write(6,*) '  1) Just the eigenvalues (bands)'
    write(6,*) '  2) Eigenvalues and kinetic energy'
    write(6,*) '  3) Eigenvalues and eigenfunctions'
    write(6,*)

    read(5,*) ipr
    write(ioreplay,*) ipr,'   print level'

    if(ipr < 0 .or. ipr > 3) then
      ipr = 1
      write(6,*)
      write(6,*) ' Wrong answer, will show only eigenvalues'
      write(6,*)
    endif

    if(ipr /= 0) then

      if(ipr == 2) then
        call kinetic_energy(neig, mtxd, ekpg, psi, ekpsi,                &
            mxddim, mxdbnd)
      endif

      nrka = -1
      call print_eig(ipr, 0, labelk, nrka, rk0,                          &
          mtxd, icmplx, neig, psi,                                       &
          adot, ei, ekpsi, isort, kgv,                                   &
          mxddim, mxdbnd, mxdgve)

    endif

    write(6,*)
    write(6,*) '  Do you want to print hamiltonian components (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   hamiltonian components'

    if(yesno == 'y' .or. yesno == 'Y') then

      write(6,*)
      write(6,*) '  How many components do you want to see?'
      write(6,*) '  Suggestion: about 20.'
      write(6,*)

      read(5,*) nsmall
      write(ioreplay,*) nsmall,'   number of components'

      if(nsmall < 1) nsmall = 20
      if(nsmall > mtxd) nsmall = mtxd

      allocate(psi_tmp(mxddim,nsmall))

      allocate(hloc(nsmall,nsmall))
      allocate(hkin(nsmall,nsmall))
      allocate(hnl(nsmall,nsmall))

      do n = 1,nsmall
        do j = 1,mxddim
          psi_tmp(j,n) = C_ZERO
        enddo
        psi_tmp(n,n) = C_UM
      enddo

      allocate(vscr_sp(mxdscr,1))
      vscr_sp(:,1) = vscr(:)
      call psi_h_psi(rk0, nsmall, psi_tmp, mtxd, isort, ekpg, 1, 1,      &
          hloc, hkin, hnl,                                               &
          ng, kgv,                                                       &
          vscr_sp, kmscr, nqnl, delqnl, vkb, nkb,                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdgve, mxddim, nsmall, mxdlqp, mxdscr)
      deallocate(vscr_sp)

      write(6,*)
      write(6,*)  '  Local component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hloc, nsmall, 1,                    &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      write(6,*)
      write(6,*)  '  Kinetic component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hkin, nsmall, 1,                    &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      write(6,*)
      write(6,*)  '  Non-Local component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hnl, nsmall, 1,                     &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      deallocate(psi_tmp)

      deallocate(hloc)
      deallocate(hkin)
      deallocate(hnl)

    endif

    write(6,*)
    write(6,*) '  Do you want to create a file for k.p calculations? (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   create k.p file'

    if(yesno == 'y' .or. yesno == 'Y') then

      allocate(h0(mxdbnd,mxdbnd))
      allocate(dh0drk(mxdbnd,mxdbnd,3))
      allocate(d2h0drk2(mxdbnd,mxdbnd,3,3))

      nder = 2

      call kdotp_matrix(mtxd, neig, psi, ei, rk0, isort, nder,           &
          h0, dh0drk, d2h0drk2,                                          &
          ng, kgv,                                                       &
          ntype, natom, rat, adot,                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      natot = 0
      do n = 1,ntype
        natot = natot + natom(n)
      enddo
      neltot = nint(ztot)
      j = neltot/2
      emidgap = (h0(j,j) + h0(j+1,j+1)) / 2

      call kdotp_silvaco_out('matrix_kdotp.dat', 23, neig, adot, rk0,    &
          h0, dh0drk, d2h0drk2, natot, neltot, emidgap, .FALSE.,         &
          mxdbnd)

      deallocate(h0)
      deallocate(dh0drk)
      deallocate(d2h0drk2)

    endif

    write(6,*)
    write(6,*) '  Do you want to calculate oscillator strengths? (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   oscillator strengths'

    if(yesno == 'y' .or. yesno == 'Y') then

      allocate(h0(mxdbnd,mxdbnd))
      allocate(dh0drk(mxdbnd,mxdbnd,3))
      allocate(d2h0drk2(mxdbnd,mxdbnd,3,3))

      nder = 2

      call kdotp_matrix(mtxd, neig, psi, ei, rk0, isort, nder,           &
          h0, dh0drk, d2h0drk2,                                          &
          ng, kgv,                                                       &
          ntype, natom, rat, adot,                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      deallocate(h0)
      deallocate(d2h0drk2)

      call out_band_oscillator_range(ioreplay, .FALSE.,                  &
            neig, ei, ztot,                                              &
            ninitbeg, ninitend, nfinalbeg, nfinalend, lpair, lexcit,     &
            mxdbnd)

      write(6,*)
      write(6,*) '  Do you want the oscillator strengths on x y z directions? (y/n)'
      write(6,*) '  Otherwise you will be asked for the desired direction.'
      write(6,*)

      read(5,*) yesno2
      write(ioreplay,*) yesno2,'   x y z directions'

      if(yesno2 == 'y' .or. yesno2 == 'Y') then

        lxyz = .TRUE.

        call out_band_oscillator_strength(neig, ei, dh0drk, adot,        &
            .FALSE., lpair, lexcit, lxyz, rdircar,                       &
            ninitbeg, ninitend, nfinalbeg, nfinalend,                    &
            mxdbnd)

      else

        lxyz = .FALSE.
        typeofr = 'direction'

! loop over directions

        do i = 1,1000

          call cpw_pp_get_r_point(rdir, rdircar, adot, typeofr, ioreplay)

          call out_band_oscillator_strength(neig, ei, dh0drk, adot,      &
             .FALSE., lpair, lexcit, lxyz, rdircar,                      &
             ninitbeg, ninitend, nfinalbeg, nfinalend,                   &
             mxdbnd)

          write(6,*)
          write(6,*) '  Do you want the oscillator strengths in a new direction? (y/n)'
          write(6,*)


          read(5,*) yesno3
          write(ioreplay,*) yesno3,'     new direction'

          if(yesno3 /= 'y' .and. yesno3 /= 'Y') exit

        enddo

      endif

      deallocate(dh0drk)

    endif

  endif

  write(6,*)
  write(6,*) '  Do you want to analyze the results with'
  write(6,*) '  spin-orbit interaction (perturbation)? (y/n)'
  write(6,*)

  read(5,*) yesno
  write(ioreplay,*) yesno,'   with SO'

  if(yesno == 'y' .or. yesno == 'Y') then

    allocate(ei_so(2*mxdbnd))
    allocate(psi_so(2*mxddim,2*mxdbnd))
    allocate(ekpsi_so(2*mxdbnd))

    call spin_orbit_perturb(rk0, mtxd, isort,                            &
        neig, psi, ei, ei_so, psi_so, .TRUE.,                            &
        ng, kgv,                                                         &
        nqnl, delqnl, vkb, nkb,                                          &
        ntype, natom, rat, adot,                                         &
        mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

    call psi_orient_spin(mtxd, neig, psi_so, ei_so,                      &
        mxddim, mxdbnd)

    write(6,*)
    write(6,*) '  Do you want to see the eigenvalues (0,1,2,3)'
    write(6,*)
    write(6,*) '  0) No'
    write(6,*) '  1) Just the eigenvalues'
    write(6,*) '  2) Eigenvalues and kinetic energy'
    write(6,*) '  3) Eigenvalues and eigenfunctions'
    write(6,*)

    read(5,*) ipr
    write(ioreplay,*) ipr,'   print level'

    if(ipr < 0 .or. ipr > 3) then
      ipr = 1
      write(6,*)
      write(6,*) ' Wrong answer, will show only eigenvalues'
      write(6,*)
    endif

    if(ipr /= 0) then
      if(ipr == 2) then
        call kinetic_energy_so(neig, mtxd, ekpg, psi_so, ekpsi_so,       &
            mxddim, mxdbnd)
      endif

      nrka = -1
      irk = 0
      call print_eig_so(ipr, irk, labelk, nrka, rk0,                     &
          mtxd, neig, psi_so,                                            &
          adot, ei_so, ekpsi_so, isort, kgv,                             &
          mxddim, mxdbnd, mxdgve)

    endif

    write(6,*)
    write(6,*) '  Do you want to print hamiltonian components (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   hamiltonian components'

    if(yesno == 'y' .or. yesno == 'Y') then

      write(6,*)
      write(6,*) '  How many components do you want to see?'
      write(6,*) '  Suggestion: about 20.'
      write(6,*)

      read(5,*) nsmall
      write(ioreplay,*) nsmall,'   number of components'

      if(nsmall < 1) nsmall = 20
      if(nsmall > mtxd) nsmall = mtxd

      allocate(psi_tmp(mxddim,nsmall))

      allocate(hloc(nsmall,nsmall))
      allocate(hkin(nsmall,nsmall))
      allocate(hnl(nsmall,nsmall))

      do n = 1,nsmall
        do j = 1,mxddim
          psi_tmp(j,n) = C_ZERO
        enddo
        psi_tmp(n,n) = C_UM
      enddo

      allocate(vscr_sp(mxdscr,1))
      vscr_sp(:,1) = vscr(:)
      call psi_h_psi(rk0, nsmall, psi_tmp, mtxd, isort, ekpg, 2, 1,      &
          hloc, hkin, hnl,                                               &
          ng, kgv,                                                       &
          vscr_sp, kmscr, nqnl, delqnl, vkb, nkb,                        &
          ntype, natom, rat, adot,                                       &
          mxdtyp, mxdatm, mxdgve, mxddim, nsmall, mxdlqp, mxdscr)
      deallocate(vscr_sp)

      write(6,*)
      write(6,*)  '  Local component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hloc, nsmall, 2,                    &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      write(6,*)
      write(6,*)  '  Kinetic component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hkin, nsmall, 2,                    &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      write(6,*)
      write(6,*)  '  Non-Local component of hamiltonian'
      write(6,*)

      call print_hamilt(icmplx, rk0, hnl, nsmall, 2,                     &
          adot, isort, kgv,                                              &
          mxddim, mxdgve, nsmall)

      deallocate(psi_tmp)

      deallocate(hloc)
      deallocate(hkin)
      deallocate(hnl)

    endif

    write(6,*)
    write(6,*) '  Do you want to create a file for k.p calculations? (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   create k.p file'

    if(yesno == 'y' .or. yesno == 'Y') then

      allocate(hso0(2*mxdbnd,2*mxdbnd))
      allocate(dhso0drk(2*mxdbnd,2*mxdbnd,3))
      allocate(d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3))

      nder = 2

      call kdotp_matrix_so_pert(mtxd, neig, psi, ei, rk0, isort, nder,   &
          hso0, dhso0drk, d2hso0drk2,                                    &
          ng, kgv,                                                       &
          ntype, natom, rat, adot,                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      call kdotp_matrix_so_convert(neig, hso0, dhso0drk, d2hso0drk2,     &
           nder,                                                         &
           mxdbnd)

      natot = 0
      do n = 1,ntype
        natot = natot + natom(n)
      enddo
      neltot = nint(ztot)
      j = neltot
      emidgap = (hso0(j,j) + hso0(j+1,j+1)) / 2

      call kdotp_silvaco_out('matrix_kdotp_so.dat', 23,                  &
          2*neig, adot, rk0,                                             &
          hso0, dhso0drk, d2hso0drk2, natot, neltot, emidgap, .TRUE.,    &
          2*mxdbnd)

      deallocate(hso0)
      deallocate(dhso0drk)
      deallocate(d2hso0drk2)

    endif

    write(6,*)
    write(6,*) '  Do you want to calculate oscillator strengths? (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,*) yesno,'   oscillator strengths'

    if(yesno == 'y' .or. yesno == 'Y') then

      allocate(hso0(2*mxdbnd,2*mxdbnd))
      allocate(dhso0drk(2*mxdbnd,2*mxdbnd,3))
      allocate(d2hso0drk2(2*mxdbnd,2*mxdbnd,3,3))

      nder = 2

      call kdotp_matrix_so_pert(mtxd, neig, psi, ei, rk0, isort, nder,   &
          hso0, dhso0drk, d2hso0drk2,                                    &
          ng, kgv,                                                       &
          ntype, natom, rat, adot,                                       &
          nqnl, delqnl, vkb, nkb,                                        &
          mxdtyp, mxdatm, mxdlqp, mxddim, mxdbnd, mxdgve)

      call kdotp_matrix_so_convert(neig, hso0, dhso0drk, d2hso0drk2,     &
           nder,                                                         &
           mxdbnd)

      deallocate(d2hso0drk2)

      allocate(vec_so(2*mxdbnd,2*mxdbnd))
      allocate(dh_so(2*mxdbnd,2*mxdbnd,3))
      allocate(tmp1(2*mxdbnd,2*mxdbnd))
      allocate(tmp2(2*mxdbnd,2*mxdbnd))

      call diag_c16(2*neig, hso0, ei_so, vec_so, 2*mxdbnd, info)


      if(info /= 0) stop


      do k = 1,3

        do j = 1,2*neig
        do i = 1,2*neig
          tmp1(i,j) = dhso0drk(i,j,k)
        enddo
        enddo

        call zgemm('n','n', 2*neig, 2*neig, 2*neig, C_UM, tmp1,          &
           2*mxdbnd, vec_so, 2*mxdbnd, C_ZERO, tmp2, 2*mxdbnd)

        call zgemm('c','n', 2*neig, 2*neig, 2*neig, C_UM, vec_so,        &
           2*mxdbnd, tmp2, 2*mxdbnd, C_ZERO, tmp1, 2*mxdbnd)

        do j = 1,2*neig
        do i = 1,2*neig
          dh_so(i,j,k) = tmp1(i,j)
        enddo
        enddo

      enddo

      call out_band_oscillator_range(ioreplay, .TRUE.,                   &
            2*neig, ei_so, ztot,                                         &
            ninitbeg, ninitend, nfinalbeg, nfinalend, lpair, lexcit,     &
            2*mxdbnd)

      write(6,*)
      write(6,*) '  Do you want the oscillator strengths on x y z directions? (y/n)'
      write(6,*) '  Otherwise you will be asked for the desired direction.'
      write(6,*)

      read(5,*) yesno2
      write(ioreplay,*) yesno2,'   x y z directions'

      if(yesno2 == 'y' .or. yesno2 == 'Y') then

        lxyz = .TRUE.

        call out_band_oscillator_strength(2*neig, ei_so, dh_so, adot,    &
            .TRUE., lpair, lexcit, lxyz, rdircar,                        &
            ninitbeg, ninitend, nfinalbeg, nfinalend,                    &
            2*mxdbnd)

      else

        lxyz = .FALSE.
        typeofr = 'direction'

        do i = 1,1000

          call cpw_pp_get_r_point(rdir, rdircar, adot, typeofr, ioreplay)

          call out_band_oscillator_strength(2*neig, ei_so, dh_so, adot,  &
              .TRUE., lpair, lexcit, lxyz, rdircar,                      &
              ninitbeg, ninitend, nfinalbeg, nfinalend,                  &
              2*mxdbnd)

          write(6,*)
          write(6,*) '  Do you want the oscillator strengths in a new direction? (y/n)'
          write(6,*)


          read(5,*) yesno3
          write(ioreplay,*) yesno3,'     new direction'

          if(yesno3 /= 'y' .and. yesno3 /= 'Y') exit

        enddo

      endif

      deallocate(vec_so)
      deallocate(dh_so)
      deallocate(tmp1)
      deallocate(tmp2)

      deallocate(hso0)
      deallocate(dhso0drk)

    endif



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

  return

end subroutine out_band_onek
