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
!>  the effective potential and charge density and writes file
!>  with the information for the atomic orbital interpolation
!>  MTB (Modified Tight Binding) method.
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      5.03
!>  \date         December 18, 2013, 29 November 2021.
!>  \copyright    GNU Public License v2

       subroutine ao_interpolation_sub(ioreplay)

!      written December 18, 2013. jlm
!      modified March 31, 2014.
!      modified May 22, 2014.
!      modified, mxdlao, December 1, 2015.  JLM
!      Modified, ioreplay, 8 January 2017. JLM
!      Modified, cpw_variables. 7 January 2020. JLM
!      Modified, epspsi, icmax, 14 June 2020. JLM
!      Modified, efermi, 29 November 2021. JLM
!      Modified, size of author, 13 January 2024.


       use cpw_variables
       use NonOrthoInterp

       implicit none

       type(noiData_t)                   ::  noiData
       logical                           ::  hsr_exists


!       integer, parameter          :: REAL64 = selected_real_kind(12)

       integer, intent(in)         :: ioreplay                           !<  tape number for reproducing calculations

       type(dims_t)                       ::  dims_                      !<  array dimensions

!      dimensions

!       integer                            ::  mxdtyp                     !  array dimension of types of atoms
!       integer                            ::  mxdatm                     !  array dimension of number of atoms of a given type

!       integer                            ::  mxdgve                     !  array dimension for g-space vectors
!       integer                            ::  mxdnst                     !  array dimension for g-space stars
!       integer                            ::  mxdcub                     !  array dimension for 3-index g-space

!       integer                            ::  mxdlqp                     !  array dimension for local potential
!       integer                            ::  mxdlao                     !  array dimension of orbital per atom type

       integer                            ::  mxdgvein                   !  array dimension for input g-space vectors
       integer                            ::  mxdnstin                   !  array dimension for input g-space stars

!      flags

       type(flags_t)                      ::  flags_                     !<  computational flags

!       character(len=6)                   ::  flgpsd
!       character(len=4)                   ::  flgdal                     !  dual approximation if equal to 'DUAL'
!       character(len=6)                   ::  flgscf                     !  type of self consistent field and diagonalizatioN

       character(len=4)                   ::  flgdalin                     !  whether the dual approximation is used

!      atomic structure variables

       type(crys_t)                       ::  crys_                      !<  crystal structure

!       real(REAL64)                       ::  alatt                      !  lattice constant

!       real(REAL64)                       ::  adot(3,3)                  !  metric in direct space
!       integer                            ::  ntype                      !  number of types of atoms
!       integer,allocatable                ::  natom(:)                   !  number of atoms of type i
!       character(len=2),allocatable       ::  nameat(:)                  !  chemical symbol for the type i
!       real(REAL64),allocatable           ::  rat(:,:,:)                 !  k-th component (in lattice coordinates) of the position of the n-th atom of type i

!      reciprocal space

       type(recip_t)                      ::  recip_                     !<  reciprocal space information

!       integer                            ::  ng                         !  total number of g-vectors with length less than gmax
!       integer, allocatable               ::  kgv(:,:)                   !  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
!       complex(REAL64), allocatable       ::  phase(:)                   !  phase factor of G-vector n
!       real(REAL64), allocatable          ::  conj(:)                    !  is -1 if one must take the complex conjugate of x*phase
!       integer                            ::  ns                         !  number os stars with length less than gmax
!       integer, allocatable               ::  inds(:)                    !  star to which g-vector n belongs
!       integer                            ::  kmax(3)                    !  max value of kgv(i,n)
!       integer, allocatable               ::  indv(:)                    !  kgv(i,indv(jadd)) is the g-vector associated with jadd jadd is defined by the g-vector components and kmax
!       integer, allocatable               ::  mstar(:)                   !  number of g-vectors in the j-th star
!       real(REAL64), allocatable          ::  ek(:)                      !  kinetic energy (Hartree) of g-vectors in star j
!       integer, allocatable               ::  izstar(:)                 !  is 0 if the phase=0

       type(recip_t)                      ::  recip_in_                    !<  reciprocal space information

!       integer                            ::  ngin                       !  input value of ng
!       integer, allocatable               ::  kgvin(:,:)                 !  input values of kgv
!       complex(REAL64), allocatable       ::  phasein(:)                 !  phase factor of G-vector n
!       real(REAL64), allocatable          ::  conjin(:)                  !  is -1 if one must take the complex conjugate of x*phase
!       integer                            ::  nsin                       !  input value of ns
!       integer                            ::  kmaxin(3)                  !  max value of kgvin(i,n)
!       integer, allocatable               ::  mstarin(:)                 !  input values of kgv

!      space group information

       type(spaceg_t)                     ::  spaceg_                    !<  space group information

!       integer                            ::  ntrans                     !  number of symmetry operations in the factor group
!       integer                            ::  mtrx(3,3,48)               !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
!       real(REAL64)                       ::  tnp(3,48)                  !  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

!        pseudo-potential (Kleinman-Bylander)


       type(pseudo_t)                     ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)

!       real(REAL64)                       ::  ealraw                   !  G=0 contrib. to the total energy. (non norm. to vcell, Hartree)
!       real(REAL64), allocatable          ::  vloc  (:,:)                !  local pseudopotential for atom k (hartree)
!       real(REAL64), allocatable          ::  dcor(:,:)                  !  core charge density for atom k
!       real(REAL64), allocatable          ::  dval (:,:)                 !  valence charge density for atom k

!       real(REAL64), allocatable          ::  zv(:)                      !  valence of atom with type i

!       real(REAL64)                       ::  ztot                       !  total charge density (electrons/cell)
!       integer, allocatable               ::  nqnl(:)                    !  number of points for pseudo interpolation for atom k
!       real(REAL64), allocatable          ::  delqnl(:)                  !  step used in the pseudo interpolation for atom k

!       real(REAL64), allocatable          ::  vkb(:,:,:,:)               !  kb nonlocal pseudo. for atom k, ang. mom. l. (not normalized to vcell, hartree)
!       integer, allocatable               ::  nkb(:,:,:)                 !   kb pseudo.  normalization for atom k, ang. mom. l

!      k-point data

       type(kpoint_t)                     ::  kpoint_                    !<  k-point data

!       integer                            ::  nx,ny,nz                   !  divisions of Brillouin zone for integration (Monkhorst-Pack)
!       real(REAL64)                       ::  sx,sy,sz                   !  shift of points in division of Brillouin zone for integration (Monkhorst-Pack)

!      atomic orbitals in G-space

       type(atorb_t)                      ::  atorb_                     !<  atomic orbitals in G-space

!       logical                            ::  latorb                     !  indicates if all atoms have information about atomic orbitals

!       integer, allocatable               ::  norbat(:)                  !  number of atomic orbitals for atom k
!       integer, allocatable               ::  lorb(:,:)                  !  angular momentum of orbital n of atom k
!       real(REAL64), allocatable          ::  wvfao(:,:,:)               !  wavefunction for atom k, ang. mom. l

!       integer, allocatable               ::  nqwf(:)                    !  number of points for wavefunction interpolation for atom k
!       real(REAL64), allocatable          ::  delqwf(:)                  !  step used in the wavefunction interpolation for atom k

       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices

!       real(REAL64)                       ::  emax                       !  kinetic energy cutoff of plane wave expansion (hartree).
!       real(REAL64)                       ::  teleck                     !  electronic temperature (in kelvin)

!       integer                            ::  nbandin                    !  target for number of bands

       type(chdens_t)                     ::  chdensin_                  !<  input charge densities

!       complex(REAL64), allocatable       ::  den(:)                   !  total charge density for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  denc(:)                  !  core charge density for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  dens(:)                  !  spherical atomic valence charge density for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  dend(:)                  !  bonding charge density from previous md step
!       complex(REAL64), allocatable       ::  dend1(:)                 !  bonding charge density from second previous md step

       type(vcomp_t)                      ::  vcompin_                   !<  local potential contributions

!       complex(REAL64), allocatable       ::  vion(:)                  !  ionic potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  vhar(:)                  !  Hartree potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  vxc(:)                   !  Hartree+exchange+correlation potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  veff(:)                  !  Effective potential for the prototype G-vector in star j

       type(vcomp_t)                      ::  vcomp_                   !<  local potential contributions

!       complex(REAL64), allocatable       ::  vion(:)                  !  ionic potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  vhar(:)                  !  Hartree potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  vxc(:)                   !  Hartree+exchange+correlation potential for the prototype G-vector in star j
!       complex(REAL64), allocatable       ::  veff(:)                  !  Effective potential for the prototype G-vector in star j

       real(REAL64)                       ::  emaxin                     !  kinetic energy cutoff of plane wave expansion (hartree).
       real(REAL64)                       ::  efermi                     !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree


       type(strfac_t)                     ::  strfac_                    !<  structure factors

!       integer                            ::  icmplx                     !  indicates if the structure factor is complex
!       complex(REAL64), allocatable       ::  strfac_%sfact(:,:)                 !  structure factor for atom k and star i

!      information about the calculation

       character(len=4)                   ::  author                     !  type of xc wanted (CA=PZ , PW92 , PBE)

       character(len=60)                  ::  pwline                     !  identifier of the calculation.  May contain miscellaneous information!
       character(len=50)                  ::  title                      !  title for plots
       character(len=140)                 ::  subtitle                   !  title for plots
       character(len=250)                 ::  meta_cpw2000               !  metadata from cpw2000

!      other variables

!       integer           ::  ipr
       integer           ::  iotape


       integer           :: iguess                                       !  kept for compatibility
       real(real64)      :: epspsi                                       !  accuracy of eigenvalues
       integer           :: icmax                                        !  maximum number of iterations for diagonalization

       character(len=60)   ::  filename

       character(len=1)   ::  yesno

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64


!      hard coded values, change if you know what you are doing...

       flags_%flgpsd = 'PSEUKB'
       iguess = 0
       epspsi = 0.0001
       icmax = 40

!      writes preamble to standard output

       write(6,*)
       write(6,'("  Analysis of the electronic structure ")')
       write(6,*)

!      open file and reads data

       filename = 'PW_RHO_V.DAT'
       iotape = 21

       call cpw_pp_band_dos_init(filename, iotape,                       &
     &    dims_, spaceg_, flags_, crys_, recip_in_, pseudo_, kpoint_,    &
     &    pwexp_, chdensin_, vcompin_, atorb_,                           &
     &    emaxin, efermi, flgdalin, author,                              &
     &    pwline, title, subtitle ,meta_cpw2000,                         &
     &    mxdgvein, mxdnstin)


       call cpw_pp_band_prepare(ioreplay,                                &
     &   dims_, crys_, spaceg_, recip_, recip_in_, pwexp_, strfac_,      &
     &   vcomp_, vcompin_,                                               &
     &   emaxin)


       write(6,*)
       write(6,*) '  Do you want to use the dual approximation (y/n)?'
       read(5,*) yesno
       write(ioreplay,*) yesno,'   dual'

       flags_%flgdal = '    '
       if(yesno == 'y' .or. yesno == 'Y') then
         flags_%flgdal = 'DUAL'
         write(6,*)
         write(6,*) '  Will use DUAL approximation'
         write(6,*)
       else
         write(6,*)
         write(6,*) '  Will not use dual approximation'
         write(6,*)
       endif



!------------------------AO INTERPOLATION ------------------------------
!              write(*,*) '  deb ztot=', ztot
             inquire (file='hsr.dat', exist=hsr_exists)
             if(hsr_exists) then
                write(*,*) '  I found a file containing ao interpolation setup data.'
                call NonOrthoInterpFileInfo(noiData)
                write(*,*) '  Do you want to use it? (y/n)'
                read(5,*) yesno
                write(ioreplay,*) yesno,'   use hsr'
              if(yesno == 'y' .or. yesno == 'Y') then
                call NonOrthoInterpReadFromFile(noiData)
              else



             call ao_interpolation_prepare(ioreplay,noiData,             &
     &       pwexp_%emax, flags_%flgdal, flags_%flgpsd,                  &
     &       iguess, epspsi, icmax, pseudo_%ztot,                        &
     &       crys_%adot, crys_%ntype, crys_%natom, crys_%rat,            &
     &       spaceg_%ntrans,spaceg_%mtrx,                                &
     &       recip_%ng, recip_%kgv, recip_%phase, recip_%conj,           &
     &       recip_%ns, recip_%inds, recip_%kmax,                        &
     &       recip_%indv, recip_%ek,                                     &
     &       strfac_%sfact, strfac_%icmplx,                              &
     &       vcomp_%veff,                                                &
     &       pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,         &
     &       atorb_%norbat, atorb_%nqwf, atorb_%delqwf,                  &
     &       atorb_%wvfao, atorb_%lorb,                                  &
     &       dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,     &
     &       dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)


              endif

             else


             call ao_interpolation_prepare(ioreplay,noiData,             &
     &       pwexp_%emax, flags_%flgdal, flags_%flgpsd,                  &
     &       iguess, epspsi, icmax, pseudo_%ztot,                        &
     &       crys_%adot, crys_%ntype, crys_%natom, crys_%rat,            &
     &       spaceg_%ntrans,spaceg_%mtrx,                                &
     &       recip_%ng, recip_%kgv, recip_%phase, recip_%conj,           &
     &       recip_%ns, recip_%inds, recip_%kmax,                        &
     &       recip_%indv, recip_%ek,                                     &
     &       strfac_%sfact, strfac_%icmplx,                              &
     &       vcomp_%veff,                                                &
     &       pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,         &
     &       atorb_%norbat, atorb_%nqwf, atorb_%delqwf,                  &
     &       atorb_%wvfao, atorb_%lorb,                                  &
     &       dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,     &
     &       dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao)


             endif


             write(*,*) '  do you want a band structure with unfolding'
             write(*,*) '  using ao interpolation? (y/n)'
             read(5,*) yesno
             write(ioreplay,*) yesno,'   band unfold'

             if(yesno == 'y' .or. yesno == 'Y') then

              write(6,*) '     band structure with unfolding  '
              write(6,*) '         (ao interpolation)         '

              call ao_interpolation_out_band_fold_full(ioreplay,         &
     &        meta_cpw2000, title, subtitle, noiData,                    &
     &        pwexp_%emax, pseudo_%ztot,                                 &
     &        crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,        &
     &        crys_%rat,                                                 &
     &        recip_%ng, recip_%kgv,                                     &
     &        strfac_%icmplx,                                            &
     &        atorb_%norbat, atorb_%nqwf, atorb_%delqwf,                 &
     &        atorb_%wvfao, atorb_%lorb,                                 &
     &        dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve,dims_%mxdlqp,     &
     &        dims_%mxdlao)

              write(6,*) '     '
              write(6,*) '     completed'

             endif

             write(*,*) '  do you want a band structure '
             write(*,*) '  using ao interpolation? (y/n)'
             read(5,*) yesno
             write(ioreplay,*) yesno,'   band'

             if(yesno == 'y' .or. yesno == 'Y') then
              write(6,*) '     band structure (ao interpolation)'

              call ao_interpolation_out_band(title,subtitle,             &
     &            noiData, pseudo_%ztot, crys_%adot)
              write(6,*) '     '
              write(6,*) '     completed'
             endif

             write(*,*) '  do you want prepare file for density of states '
             write(*,*) '  using ao interpolation? (y/n)'
             read(5,*) yesno
             write(ioreplay,*) yesno,'   dos prepare'

             if(yesno == 'y' .or. yesno == 'Y') then
              write(6,*) '     files for density of states (ao interpolation)'
              call ao_interpolation_out_dos(noiData, pseudo_%ztot, crys_%adot, spaceg_%ntrans, spaceg_%mtrx)
              write(6,*) '     '
              write(6,*) '     completed'
             endif


             if(noiData % loptical == 1) then
              write(*,*) '  do you want prepare file for the dielectric function '
              write(*,*) '  using ao interpolation? (y/n)'
              read(5,*) yesno
              write(ioreplay,*) yesno,'   ie_prepare'

              if(yesno == 'y' .or. yesno == 'Y') then
                write(6,*) '     files for dielectric function calculation (ao interpolation)'
                call ao_interpolation_out_ie(noiData, pseudo_%ztot, crys_%adot, spaceg_%ntrans, spaceg_%mtrx)
                write(6,*) '     '
                write(6,*) '     completed'
              endif
             endif

              write(*,*) '  do you want the files for Wannier90 run? '
              read(5,*) yesno
              write(ioreplay,*) yesno,'   w90_interface'

              if(yesno == 'y' .or. yesno == 'Y') then

              call ao_interpolation_w90(noiData,                         &
     &        pwexp_%emax, pseudo_%ztot,                                 &
     &        crys_%adot,                                                &
     &        crys_%ntype, crys_%natom, crys_%nameat, crys_%rat,         &
     &        recip_%ng, recip_%kgv,                                     &
     &        strfac_%icmplx,                                            &
     &        atorb_%norbat, atorb_%nqwf, atorb_%delqwf,                 &
     &        atorb_%wvfao, atorb_%lorb,                                 &
     &        dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve,dims_%mxdlqp,     &
     &        dims_%mxdlao)
             endif






             call NonOrthoInterpFinish(noiData)

             write(*,*) '  do you want prepare file compatible with WannierTools? '
             read(5,*) yesno
             write(ioreplay,*) yesno,'   write hr'
             if(yesno == 'y' .or. yesno == 'Y') then
                call ao_interpolation_write_hr(crys_%adot)
                call ao_interpolation_write_wt_in(noiData,pseudo_%ztot,crys_%adot,crys_%ntype,crys_%natom, &
     &               crys_%nameat, crys_%rat, atorb_%norbat, atorb_%lorb,             &
     &               dims_%mxdtyp, dims_%mxdatm, dims_%mxdlao)


                write(6,*) '     '
                write(6,*) '     completed'
             endif


             write(*,*) '  All Done!'
             stop
!------------------------------------------------------------------------


       return
       end subroutine ao_interpolation_sub

