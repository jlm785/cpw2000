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

!>  This subroutines reads a file with the atomic structure and
!>  the effective potential, calculates the wave-functions for a
!>  k-point and makes 1D, 2D, and 3D plots.
!>
!>  \author       Jose Luis Martins
!>  \version      5.12
!>  \date         16 February 2018, 22 October 2025.
!>  \copyright    GNU Public License v2


subroutine plot_psi_sub(ioreplay)

! written 16 February 2018 based on cpw_analysis_sub, out_band_onek,
! and rho_v_plot_sub.
! Modernized, documentation, APIs, 2 February 2021. JLM
! Modified, efermi, 29 November 2021. JLM
! Modified, size of author, 13 January 2024.
! Modified, enter method for k-point. 22 October 2025. JLM


  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)         :: ioreplay                                !  tape number for reproducing calculations

  type(dims_t)                       ::  dims_                           !<  array dimensions

! dimensions

!  integer                            ::  mxdtyp                          !  array dimension of types of atoms
!  integer                            ::  mxdatm                          !  array dimension of number of atoms of a given type

!  integer                            ::  mxdgve                          !  array dimension for g-space vectors
!  integer                            ::  mxdnst                          !  array dimension for g-space stars
!  integer                            ::  mxdcub                          !  array dimension for 3-index g-space

!  integer                            ::  mxdlqp                          !  array dimension for local potential
!  integer                            ::  mxdlao                          !  array dimension of orbital per atom type

  type(dims_t)                       ::  dims_in_                        !<  input array dimensions

! flags

  type(flags_t)                      ::  flags_                          !<  computational flags

!  character(len=6)                   ::  flgpsd
!  character(len=4)                   ::  flgdal                          !  dual approximation if equal to 'DUAL'
!  character(len=6)                   ::  flgscf                          !  type of self consistent field and diagonalizatioN

  character(len=4)                   ::  flgdal_in                       !  whether the dual approximation is used

! atomic structure variables

  type(crys_t)                       ::  crys_                           !<  crystal structure

!  real(REAL64)                       ::  alatt                           !  lattice constant

!  real(REAL64)                       ::  adot(3,3)                       !  metric in direct space
!  integer                            ::  ntype                           !  number of types of atoms
!  integer,allocatable                ::  natom(:)                        !  number of atoms of type i
!  character(len=2),allocatable       ::  nameat(:)                       !  chemical symbol for the type i
!  real(REAL64),allocatable           ::  rat(:,:,:)                      !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! reciprocal space

  type(recip_t)                      ::  recip_                          !<  reciprocal space information

!  integer                            ::  ng                              !  total number of g-vectors with length less than gmax
!  integer, allocatable               ::  kgv(:,:)                        !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
!  complex(REAL64), allocatable       ::  phase(:)                        !  phase factor of G-vector n
!  real(REAL64), allocatable          ::  conj(:)                         !  is -1 if one must take the complex conjugate of x*phase
!  integer                            ::  ns                              !  number os stars with length less than gmax
!  integer, allocatable               ::  inds(:)                         !  star to which g-vector n belongs
!  integer                            ::  kmax(3)                         !  max value of kgv(i,n)
!  integer, allocatable               ::  indv(:)                         !  kgv(i,indv(jadd)) is the g-vector associated with jadd jadd is defined by the g-vector components and kmax
!  integer, allocatable               ::  mstar(:)                        !  number of g-vectors in the j-th star
!  real(REAL64), allocatable          ::  ek(:)                           !  kinetic energy (Hartree) of g-vectors in star j
!  integer, allocatable               ::  izstar(:)                       !  is 0 if the phase=0

  type(recip_t)                      ::  recip_in_                       !<  reciprocal space information

!  integer                            ::  ngin                            !  input value of ng
!  integer, allocatable               ::  kgvin(:,:)                      !  input values of kgv
!  complex(REAL64), allocatable       ::  phasein(:)                      !  phase factor of G-vector n
!  real(REAL64), allocatable          ::  conjin(:)                       !  is -1 if one must take the complex conjugate of x*phase
!  integer                            ::  nsin                            !  input value of ns
!  integer                            ::  kmaxin(3)                       !  max value of kgvin(i,n)
!  integer, allocatable               ::  mstarin(:)                      !  input values of kgv

! space group information

  type(spaceg_t)                     ::  spaceg_                         !<  space group information

!  integer                            ::  ntrans                          !  number of symmetry operations in the factor group
!  integer                            ::  mtrx(3,3,48)                    !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
!  real(REAL64)                       ::  tnp(3,48)                       !  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!   pseudo-potential (Kleinman-Bylander)


  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)

!  real(REAL64)                       ::  ealraw                          !  G=0 contrib. to the total energy. (non norm. to vcell, Hartree)
!  real(REAL64), allocatable          ::  vloc  (:,:)                     !  local pseudopotential for atom k (hartree)
!  real(REAL64), allocatable          ::  dcor(:,:)                       !  core charge density for atom k
!  real(REAL64), allocatable          ::  dval (:,:)                      !  valence charge density for atom k

!  real(REAL64), allocatable          ::  zv(:)                           !  valence of atom with type i

!  real(REAL64)                       ::  ztot                            !  total charge density (electrons/cell)
!  integer, allocatable               ::  nqnl(:)                         !  number of points for pseudo interpolation for atom k
!  real(REAL64), allocatable          ::  delqnl(:)                       !  step used in the pseudo interpolation for atom k

!  real(REAL64), allocatable          ::  vkb(:,:,:,:)                    !  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
!  integer, allocatable               ::  nkb(:,:,:)                      !   kb pseudo.  normalization for atom k, ang. mom. l

! k-point data

  type(kpoint_t)                     ::  kpoint_                         !<  k-point data

!  integer                            ::  nx,ny,nz                        !  divisions of Brillouin zone for integration (Monkhorst-Pack)
!  real(REAL64)                       ::  sx,sy,sz                        !  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)

! atomic orbitals in G-space

  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space

!  logical                            ::  latorb                          !  indicates if all atoms have information about atomic orbitals

!  integer, allocatable               ::  norbat(:)                       !  number of atomic orbitals for atom k
!  integer, allocatable               ::  lorb(:,:)                       !  angular momentum of orbital n of atom k
!  real(REAL64), allocatable          ::  wvfao(:,:,:)                    !  wavefunction for atom k, ang. mom. l

!  integer, allocatable               ::  nqwf(:)                         !  number of points for wavefunction interpolation for atom k
!  real(REAL64), allocatable          ::  delqwf(:)                       !  step used in the wavefunction interpolation for atom k

  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices

!  real(REAL64)                       ::  emax                            !  kinetic energy cutoff of plane wave expansion (hartree).
!  real(REAL64)                       ::  teleck                          !  electronic temperature (in kelvin)

!  integer                            ::  nbandin                         !  target for number of bands

  type(chdens_t)                     ::  chdens_in_                       !<  input charge densities

!  complex(REAL64), allocatable       ::  den(:)                          !  total charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  denc(:)                         !  core charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  dens(:)                         !  spherical atomic valence charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  dend(:)                         !  bonding charge density from previous md step
!  complex(REAL64), allocatable       ::  dend1(:)                        !  bonding charge density from second previous md step

  type(vcomp_t)                      ::  vcomp_in_                       !<  local potential contributions

!  complex(REAL64), allocatable       ::  vion(:)                         !  ionic potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vhar(:)                         !  Hartree potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vxc(:)                          !  Hartree+exchange+correlation potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  veff(:)                         !  Effective potential for the prototype G-vector in star j

  type(vcomp_t)                      ::  vcomp_                          !<  local potential contributions

!  complex(REAL64), allocatable       ::  vion(:)                         !  ionic potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vhar(:)                         !  Hartree potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vxc(:)                          !  Hartree+exchange+correlation potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  veff(:)                         !  Effective potential for the prototype G-vector in star j

  real(REAL64)                       ::  emax_in                         !  input kinetic energy cutoff of plane wave expansion (hartree).


  type(strfac_t)                     ::  strfac_                         !<  structure factors

!  integer                            ::  icmplx                          !  indicates if the structure factor is complex
!  complex(REAL64), allocatable       ::  strfac_%sfact(:,:)              !  structure factor for atom k and star i

! information about the calculation

  character(len=4)                   ::  author                          !  type of xc wanted (CA=PZ , PW92 , PBE)

  character(len=60)                  ::  pwline                          !  identifier of the calculation.  May contain miscellaneous information!
  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  title for plots
  character(len=250)                 ::  meta_cpw2000                    !  metadata from cpw2000

! allocatable arrays

  real(REAL64), allocatable          ::  vscr(:)                         !  screened potential in the fft real space mesh
  real(REAL64), allocatable          ::  ei(:)                           !  eigenvalue no. i. (hartree)
  real(REAL64), allocatable          ::  hdiag(:)                        !  hamiltonian diagonal
  integer, allocatable               ::  isort(:)                        !  g-vector associated with row/column i of hamiltonian
  real(REAL64), allocatable          ::  qmod(:)                         !  length of k+g-vector of row/column i

  real(REAL64), allocatable          ::  ekpg(:)                         !  kinetic energy (hartree) of k+g-vector of row/column i
  complex(REAL64), allocatable       ::  psi(:,:)                        !  component j of eigenvector i! other variables
  complex(REAL64), allocatable       ::  hpsi(:,:)                       !  H | psi>
  real(REAL64), allocatable          ::  ekpsi(:)                        !  kinetic energy of eigenvector i. (hartree)  integer           ::  iotape

  complex(REAL64), allocatable       ::  psi_so(:,:)                     !  component j of eigenvector i!      allocatable local arrays
  real(REAL64), allocatable          ::  ei_so(:)                        !  spin-orbit eigenvalue (hartree)
  real(REAL64), allocatable          ::  ekpsi_so(:)                     !  kinetic energy of eigenvector i. (hartree)

! other variables

  integer           ::  ipr

  integer           ::  neig                                             !  number of eigenvectors required (maybe modified on output)

  integer           ::  kmscr(7)                                         !  max value of kgv(i,n) used for the potential fft mesh
  integer           ::  idshift                                          !  shift of the fft mesh, used /= 0 only in highly banked memory.
  integer           ::  nsfft(3)

  integer           ::  mxddim                                           !  array dimension for the hamiltonian
  integer           ::  mxdbnd                                           !  array dimension for the number of bands
  integer           ::  mxdscr                                           !  array dimension for screening potential
  integer           ::  mxdwrk                                           !  array dimension for fft transform workspace

  integer           ::  mtxd

  real(REAL64)           ::  vmax, vmin                                  !  maximum and minimum values of vscr

  character(len=4)       ::  diag_type                                   !  selects diagonalization, 'pw  ','ao  ','aojc'
  integer                ::  nocc
  integer                ::  ifail
  integer                ::  icmax
  integer                ::  nrka
  character(len=5)       ::  labelk

  integer                ::  iguess
  real(real64)           ::  epspsi
  character(len=4)       ::  flgdal                                      !  dual approximation if equal to 'DUAL'
  real(REAL64)           ::  efermi                                      !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  real(REAL64)           ::  rk0(3)

  real(REAL64)           ::  rkcar(3)
  character(len=20)      ::  typeofk

  integer                ::  iotape
  character(len=60)      ::  filename

  character(len=1)   ::  yesno


! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counter

  integer           ::  i

! writes preamble to standard output

  write(6,*)
  write(6,'("  Plot of the wave-functions for a given k-point")')
  write(6,*)

! open file and reads data

  filename = 'PW_RHO_V.DAT'
  iotape = 21

  call cpw_pp_band_dos_init(filename, iotape,                            &
       dims_, crys_, spaceg_, flags_, pwexp_, pseudo_, kpoint_,          &
       atorb_, efermi,  author,                                          &
       pwline, title, subtitle ,meta_cpw2000,                            &
       dims_in_, recip_in_, chdens_in_, vcomp_in_, emax_in, flgdal_in)

! prepares calculation

  call cpw_pp_band_prepare(ioreplay,                                     &
       dims_, crys_, spaceg_, recip_, pwexp_, strfac_,  vcomp_,          &
       dims_in_, recip_in_, vcomp_in_, emax_in)

!   calculates the local potential

  flags_%flgpsd = 'PSEUKB'
  iguess = 0

! for plots precision does not need to be large.  Modify if you want higher precision

  epspsi = 0.0005

    flgdal = 'DUAL'

    kmscr(1) = recip_%kmax(1)/2 + 2
    kmscr(2) = recip_%kmax(2)/2 + 2
    kmscr(3) = recip_%kmax(3)/2 + 2

  call size_fft(kmscr,nsfft,mxdscr,mxdwrk)

  allocate(vscr(mxdscr))

  ipr = 1
  idshift = 0

  call pot_local(ipr, vscr, vmax, vmin, vcomp_%veff, kmscr, idshift,     &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                  &
      recip_%ns, recip_%inds,                                            &
      mxdscr, dims_%mxdgve, dims_%mxdnst)



  do i = 1,100

!   finds mxddim, mxdbnd

    write(6,*)
    write(6,*) '  enter number of bands'
    write(6,*)
    write(6,'("   if the number is zero or negative exits  the ",        &
              & " wave-function plot")')
    write(6,*)
    write(6,'("   the number of electrons is: ",i5)') nint(pseudo_%ztot)
    write(6,'("   in doubt enter: ",i5)') nint(pseudo_%ztot)/2 + 10
    write(6,*)
    read(5,*) neig
    write(ioreplay,'(2x,i10,"   number of bands")') neig

    if(neig < 1) then

      exit

    endif

    mxdbnd = neig

    typeofk = 'reference k-point'

    call cpw_pp_get_k_vector(rk0, rkcar, crys_%adot, typeofk, ioreplay)

    call size_mtxd(pwexp_%emax, rk0, crys_%adot, recip_%ng, recip_%kgv, mtxd)

    mxddim = mtxd

!   allocates arrays

    allocate(ei(mxdbnd))
    allocate(hdiag(mxddim))
    allocate(isort(mxddim))
    allocate(qmod(mxddim))
    allocate(ekpg(mxddim))
    allocate(psi(mxddim,mxdbnd))
    allocate(hpsi(mxddim,mxdbnd))
    allocate(ekpsi(mxdbnd))

!   diagonalizes the ahmiltoniana

    ipr = 1
    diag_type = 'pw  '
    nocc = neig
    icmax = 100

    call h_kb_dia_all(diag_type, pwexp_%emax, rk0, neig, nocc,           &
        flags_%flgpsd, ipr, ifail, icmax, iguess, epspsi,                &
        recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                &
        recip_%ns, recip_%inds, recip_%kmax,                             &
        recip_%indv, recip_%ek,                                          &
        strfac_%sfact, vcomp_%veff, strfac_%icmplx,                      &
        pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,              &
        crys_%ntype,crys_%natom,crys_%rat,crys_%adot,                    &
        mtxd, hdiag, isort, qmod, ekpg, .FALSE.,                         &
        psi, hpsi, ei,                                                   &
        vscr, kmscr,                                                     &
        atorb_%latorb, atorb_%norbat, atorb_%nqwf,                       &
        atorb_%delqwf, atorb_%wvfao, atorb_%lorb,                        &
        dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,          &
        dims_%mxdcub, dims_%mxdlqp, mxddim, mxdbnd,                      &
        mxdscr, dims_%mxdlao)


!   plots

    write(6,*)
    write(6,*) '  Do you want to plot the orbitals with'
    write(6,*) '  spin-orbit interaction (perturbation)? (y/n)'
    write(6,*)

    read(5,*) yesno
    write(ioreplay,'(2x,a1,"   with SO")') yesno

    if(yesno == 'y' .or. yesno == 'Y') then

      write(6,*)
      write(6,*) '  preparing plot with spin-orbit'
      write(6,*)

      allocate(ei_so(2*mxdbnd))
      allocate(psi_so(2*mxddim,2*mxdbnd))
      allocate(ekpsi_so(2*mxdbnd))

      call spin_orbit_perturb(rk0, mtxd, isort,                          &
          neig, psi, ei, ei_so, psi_so, .TRUE.,                          &
          recip_%ng, recip_%kgv,                                         &
          pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,            &
          crys_%ntype, crys_%natom, crys_%rat,crys_% adot,               &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdlqp,                      &
          mxddim, mxdbnd, dims_%mxdgve)

      if(ipr > 0) then
        nrka = -1
        ipr = 1
        labelk = '     '
        call print_eig_so(ipr, 0, labelk, nrka, rk0,                     &
        mtxd, neig, psi_so,                                              &
        crys_%adot, ei_so, ekpsi_so, isort, recip_%kgv,                  &
        mxddim, mxdbnd, dims_%mxdgve)
      endif


      call plot_psi_plotit(ioreplay, 2,                                  &
          crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,            &
          crys_%rat, pseudo_%zv,                                         &
          recip_%ng, recip_%kgv, recip_%kmax,                            &
          mtxd, neig, isort, psi_so, ei_so,                              &
          mxddim, mxdbnd, dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve)

      deallocate(ei_so)
      deallocate(psi_so)
      deallocate(ekpsi_so)

    else

      write(6,*)
      write(6,*) '  preparing plot without spin-orbit'
      write(6,*)

      nrka = -1
      ipr = 1
      labelk = '     '
      call print_eig(ipr, 0, labelk, nrka, rk0,                          &
      mtxd, strfac_%icmplx, neig, psi,                                   &
      crys_%adot, ei, ekpsi, isort, recip_%kgv,                          &
      mxddim, mxdbnd, dims_%mxdgve)

      call plot_psi_plotit(ioreplay, 1,                                  &
          crys_%adot, crys_%ntype, crys_%natom,crys_%nameat,             &
          crys_%rat, pseudo_%zv,                                         &
          recip_%ng, recip_%kgv, recip_%kmax,                            &
          mtxd, neig, isort, psi, ei,                                    &
          mxddim, mxdbnd, dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve)

    endif

    deallocate(ei)
    deallocate(hdiag)
    deallocate(isort)
    deallocate(qmod)
    deallocate(ekpg)
    deallocate(psi)
    deallocate(hpsi)
    deallocate(ekpsi)


  enddo

  deallocate(vscr)

  return

end subroutine plot_psi_sub

