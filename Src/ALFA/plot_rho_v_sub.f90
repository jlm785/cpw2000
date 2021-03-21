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

! This subroutine calculates the layer average and double average of
! the charge density, or any other scalar periodic quantity.

! For the double average see  PRL 61, 734 (1988).

! writen September 5, 2012.jlm
! Modified, mxdlao, December 1, 2015. JLM
! Modified, May 2020, cpw_variables. CLR
! copyright  Jose Luis Martins, Carlos Loia Reis/INESC-MN

! version 4.96

subroutine plot_rho_v_sub(ioreplay)



  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)         :: ioreplay                           !  tape number for reproducing calculations

  type(dims_t)                       ::  dims_                      !<  array dimensions

! dimensions

!  integer                            ::  mxdtyp                     !  array dimension of types of atoms
!  integer                            ::  mxdatm                     !  array dimension of number of atoms of a given type

!  integer                            ::  mxdgve                     !  array dimension for g-space vectors
!  integer                            ::  mxdnst                     !  array dimension for g-space stars
!  integer                            ::  mxdcub                     !  array dimension for 3-index g-space

!  integer                            ::  mxdlqp                     !  array dimension for local potential
!  integer                            ::  mxdlao                     !  array dimension of orbital per atom type

  integer                            ::  mxdgvein                   !  array dimension for input g-space vectors
  integer                            ::  mxdnstin                   !  array dimension for input g-space stars

! flags

  type(flags_t)                      ::  flags_                     !<  computational flags

!  character(len=6)                   ::  flgpsd
!  character(len=4)                   ::  flgdal                     !  dual approximation if equal to 'DUAL'
!  character(len=6)                   ::  flgscf                     !  type of self consistent field and diagonalizatioN

  character(len=4)                   ::  flgdalin                     !  whether the dual approximation is used

! atomic structure variables

  type(crys_t)                       ::  crys_                      !<  crystal structure

!  real(REAL64)                       ::  alatt                      !  lattice constant

!  real(REAL64)                       ::  adot(3,3)                  !  metric in direct space
!  integer                            ::  ntype                      !  number of types of atoms
!  integer,allocatable                ::  natom(:)                   !  number of atoms of type i
!  character(len=2),allocatable       ::  nameat(:)                  !  chemical symbol for the type i
!  real(REAL64),allocatable           ::  rat(:,:,:)                 !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

! reciprocal space

  type(recip_t)                      ::  recip_                     !<  reciprocal space information

!  integer                            ::  ng                         !  total number of g-vectors with length less than gmax
!  integer, allocatable               ::  kgv(:,:)                   !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
!  complex(REAL64), allocatable       ::  phase(:)                   !  phase factor of G-vector n
!  real(REAL64), allocatable          ::  conj(:)                    !  is -1 if one must take the complex conjugate of x*phase
!  integer                            ::  ns                         !  number os stars with length less than gmax
!  integer, allocatable               ::  inds(:)                    !  star to which g-vector n belongs
!  integer                            ::  kmax(3)                    !  max value of kgv(i,n)
!  integer, allocatable               ::  indv(:)                    !  kgv(i,indv(jadd)) is the g-vector associated with jadd jadd is defined by the g-vector components and kmax
!  integer, allocatable               ::  mstar(:)                   !  number of g-vectors in the j-th star
!  real(REAL64), allocatable          ::  ek(:)                      !  kinetic energy (Hartree) of g-vectors in star j
!  integer, allocatable               ::  izstar(:)                 !  is 0 if the phase=0

  type(recip_t)                      ::  recip_in_                    !<  reciprocal space information

!  integer                            ::  ngin                       !  input value of ng
!  integer, allocatable               ::  kgvin(:,:)                 !  input values of kgv
!  complex(REAL64), allocatable       ::  phasein(:)                 !  phase factor of G-vector n
!  real(REAL64), allocatable          ::  conjin(:)                  !  is -1 if one must take the complex conjugate of x*phase
!  integer                            ::  nsin                       !  input value of ns
!  integer                            ::  kmaxin(3)                  !  max value of kgvin(i,n)
!  integer, allocatable               ::  mstarin(:)                 !  input values of kgv

! space group information

  type(spaceg_t)                     ::  spaceg_                    !<  space group information

!  integer                            ::  ntrans                     !  number of symmetry operations in the factor group  
!  integer                            ::  mtrx(3,3,48)               !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
!  real(REAL64)                       ::  tnp(3,48)                  !  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!   pseudo-potential (Kleinman-Bylander)


  type(pseudo_t)                     ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)

!  real(REAL64)                       ::  ealraw                   !<  G=0 contrib. to the total energy. (non norm. to vcell, Hartree)
!  real(REAL64), allocatable          ::  vloc  (:,:)                !  local pseudopotential for atom k (hartree)
!  real(REAL64), allocatable          ::  dcor(:,:)                  !  core charge density for atom k
!  real(REAL64), allocatable          ::  dval (:,:)                 !  valence charge density for atom k

!  real(REAL64), allocatable          ::  zv(:)                      !  valence of atom with type i

!  real(REAL64)                       ::  ztot                       !  total charge density (electrons/cell)
!  integer, allocatable               ::  nqnl(:)                    !  number of points for pseudo interpolation for atom k
!  real(REAL64), allocatable          ::  delqnl(:)                  !  step used in the pseudo interpolation for atom k

!  real(REAL64), allocatable          ::  vkb(:,:,:,:)               !  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
!  integer, allocatable               ::  nkb(:,:,:)                 !   kb pseudo.  normalization for atom k, ang. mom. l

! k-point data

  type(kpoint_t)                     ::  kpoint_                    !<  k-point data

!  integer                            ::  nx,ny,nz                   !  divisions of Brillouin zone for integration (Monkhorst-Pack)
!  real(REAL64)                       ::  sx,sy,sz                   !  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)

! atomic orbitals in G-space

  type(atorb_t)                      ::  atorb_                     !<  atomic orbitals in G-space

!  logical                            ::  latorb                     !  indicates if all atoms have information about atomic orbitals

!  integer, allocatable               ::  norbat(:)                  !  number of atomic orbitals for atom k
!  integer, allocatable               ::  lorb(:,:)                  !  angular momentum of orbital n of atom k
!  real(REAL64), allocatable          ::  wvfao(:,:,:)               !  wavefunction for atom k, ang. mom. l

!  integer, allocatable               ::  nqwf(:)                    !  number of points for wavefunction interpolation for atom k
!  real(REAL64), allocatable          ::  delqwf(:)                  !  step used in the wavefunction interpolation for atom k

  type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices

!  real(REAL64)                       ::  emax                       !  kinetic energy cutoff of plane wave expansion (hartree).
!  real(REAL64)                       ::  teleck                     !  electronic temperature (in kelvin)

!  integer                            ::  nbandin                    !  target for number of bands

  type(chdens_t)                     ::  chdensin_                  !<  input charge densities

!  complex(REAL64), allocatable       ::  den(:)                   !<  total charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  denc(:)                  !<  core charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  dens(:)                  !<  spherical atomic valence charge density for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  dend(:)                  !<  bonding charge density from previous md step
!  complex(REAL64), allocatable       ::  dend1(:)                 !<  bonding charge density from second previous md step

  type(vcomp_t)                      ::  vcompin_                   !<  local potential contributions

!  complex(REAL64), allocatable       ::  vion(:)                  !<  ionic potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vhar(:)                  !<  Hartree potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vxc(:)                   !<  Hartree+exchange+correlation potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  veff(:)                  !<  Effective potential for the prototype G-vector in star j

  type(vcomp_t)                      ::  vcomp_                   !<  local potential contributions

!  complex(REAL64), allocatable       ::  vion(:)                  !<  ionic potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vhar(:)                  !<  Hartree potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  vxc(:)                   !<  Hartree+exchange+correlation potential for the prototype G-vector in star j
!  complex(REAL64), allocatable       ::  veff(:)                  !<  Effective potential for the prototype G-vector in star j

  real(REAL64)                       ::  emaxin                     !  kinetic energy cutoff of plane wave expansion (hartree).


  type(strfac_t)                     ::  strfac_                    !<  structure factors

!  integer                            ::  icmplx                     !  indicates if the structure factor is complex
!  complex(REAL64), allocatable       ::  strfac_%sfact(:,:)                 !  structure factor for atom k and star i
 
! information about the calculation

  character(len=3)                   ::  author                     !  type of xc wanted (CA=PZ , PW92 , PBE)

  character(len=60)                  ::  pwline                     !  identifier of the calculation.  May contain miscellaneous information!
  character(len=50)                  ::  title                      !  title for plots
  character(len=140)                 ::  subtitle                   !  title for plots
  character(len=250)                 ::  meta_cpw2000               !  metadata from cpw2000


  complex(REAL64), allocatable       ::  qplot(:)                   !  Quantity to be plotted (in stars)
  complex(REAL64), allocatable       ::  qunsym(:)                  !  Quantity to be plotted (in G-vectors)

! other variables

  integer           ::  iotape
  integer           ::  ifunc
  integer           ::  istar,istop  
  
  integer iplot
  
  
  character(len=60)   ::  filename


! constants

  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counter

  integer           ::  i, j, n

! writes preamble to standard output

  write(6,*)
  write(6,'("  Analysis of the electronic structure ")')
  write(6,*)

! open file and reads data

  filename = 'PW_RHO_V.DAT'
  iotape = 21

  call cpw_pp_band_dos_init(filename, iotape,                            &
     dims_, spaceg_, flags_, crys_, recip_in_, pseudo_, kpoint_,         &
     pwexp_, chdensin_, vcompin_, atorb_,                                &
     emaxin, flgdalin, author,                                           &
     pwline, title, subtitle ,meta_cpw2000,                              &
     mxdgvein, mxdnstin)

  
  call cpw_pp_band_prepare(ioreplay,                                     &
    dims_, crys_, spaceg_, recip_, recip_in_, pwexp_, strfac_,           &
    vcomp_, vcompin_,                                                    &
    emaxin)


  allocate(qplot(dims_%mxdnst))
  allocate(qunsym(dims_%mxdgve))


  do i=1,100
  
    write(6,*)
    write(6,*) '  What kind of plot do you want?'
    write(6,*)
    write(6,*) '  0)  Exit'
    write(6,*) '  1)  One-dimensional planar and z-average ',            &
       'charge and electrostatic potential'
    write(6,*) '  2)  Simple one dimensional plot'
    write(6,*) '  3)  Two dimensional (contour) plot'
    write(6,*) '  4)  Three dimensional plot'
    write(6,*)
    
    read(5,*) iplot
    write(ioreplay,'(2x,i8,"   plot dimension")') iplot
    
    if(iplot < 0 .or. iplot > 4) then

      iplot = 0
      write(6,*)
      write(6,*) '  Wrong answer '
      write(6,*)

    endif


    if(iplot == 0) then

      exit

    elseif(iplot == 1) then

      call plot_rho_v_average(ioreplay,                                  &
      crys_%adot, crys_%ntype, crys_%natom,                              &
      crys_%nameat, crys_%rat, pseudo_%zv,                               &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                  &
      recip_%ns, recip_%mstar, chdensin_%den,                            &
      dims_%mxdtyp,dims_%mxdatm,dims_%mxdgve,dims_%mxdnst)
  
    else

      write(6,*)
      write(6,*) '  What function do you want to plot?'
      write(6,*)
      write(6,*) '  1)  Charge density'
      write(6,*) '  2)  Bonding charge density'
      write(6,*) '  3)  Effective potential V_H + V_xc + V_ion'
      write(6,*) '  4)  Hartree potential V_H'
      write(6,*) '  5)  Exchange correlation potential V_xc'
      write(6,*) '  6)  Ionic potential V_ion'
    
      read(5,*) ifunc
      write(ioreplay,'(2x,i8,"   function to plot")') ifunc
    
      if(ifunc < 0 .or. ifunc > 6) then

        ifunc = 0
        write(6,*)
        write(6,*) '  Wrong answer'
        write(6,*)

      endif

      if(ifunc == 0) then

        exit
  
      elseif(ifunc == 1) then
        do j = 1,recip_%ns
          qplot(j) = chdensin_%den(j)
        enddo
      elseif(ifunc == 2) then
        do j = 1,recip_%ns
          qplot(j) = chdensin_%dend(j)
        enddo
      elseif(ifunc == 3) then
        do j = 1,recip_%ns
          qplot(j) = vcompin_%veff(j)
        enddo
      elseif(ifunc == 4) then
        do j = 1,recip_%ns
          qplot(j) = vcompin_%vhar(j)
        enddo
      elseif(ifunc == 5) then
        do j = 1,recip_%ns
          qplot(j) = vcompin_%vxc(j)
        enddo
      elseif(ifunc == 6) then
        do j = 1,recip_%ns
          qplot(j) = vcompin_%vion(j)
        enddo
      endif


      istop = 0
      do n = 1,recip_%ns
        istar = istop+1
        istop = istar+recip_%mstar(n)-1
        do j=istar,istop
          qunsym(j) = qplot(n)*conjg(recip_%phase(j))
          if(recip_%conj(j) < ZERO) qunsym(j) = conjg(qunsym(j))
        enddo
      enddo
      

      if(iplot == 2) then

        call plot_average_simple(ioreplay, qunsym,                       &
            crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,          &
            crys_%rat, pseudo_%zv,                                       &
            recip_%ng, recip_%kgv,                                       &
            dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve)

      elseif(iplot == 3) then

        call plot_contour(ioreplay, qunsym,                              &
            crys_%adot, recip_%ng, recip_%kgv,                           &
            dims_%mxdgve)

      elseif(iplot == 4) then

        call plot_contour3D(ioreplay, qunsym,                            &
          crys_%adot, crys_%ntype, crys_%natom, crys_%nameat, crys_%rat, &
          recip_%ng, recip_%kgv,                                         &
          dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve)

      endif

    endif

  enddo

! deallocates the stuff

  
  deallocate(qplot)
  deallocate(qunsym)


  return
end subroutine plot_rho_v_sub

