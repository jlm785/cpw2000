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

!>  Reads the file with the atomic structure, the effective potential
!>  and the pseudopotentials (PW_RHO_V.DAT), calculates the modified
!>  tight-binding (MTB) hamiltonian and writes it to the file hsr2.dat.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         24 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_interpolation_create(ioreplay)

! Written 24 October 2026 from ao_interpolation_sub (Carlos Loia Reis,
! Jose Luis Martins and others). claude

  use cpw_variables
  use NonOrthoInterp

  implicit none

! input

  integer, intent(in)                ::  ioreplay                        !<  tape number for reproducing calculations

! variables

  type(noiData_t)                    ::  noiData                         !  see NonOrthoInterp

  type(dims_t)                       ::  dims_                           !  array dimensions
  type(dims_t)                       ::  dims_in_                        !  input array dimensions
  type(flags_t)                      ::  flags_                          !  computational flags
  type(crys_t)                       ::  crys_                           !  crystal structure
  type(recip_t)                      ::  recip_                          !  reciprocal space information
  type(recip_t)                      ::  recip_in_                       !  input reciprocal space information
  type(spaceg_t)                     ::  spaceg_                         !  space group information
  type(pseudo_t)                     ::  pseudo_                         !  pseudo-potential (Kleinman-Bylander)
  type(kpoint_t)                     ::  kpoint_                         !  k-point data
  type(atorb_t)                      ::  atorb_                          !  atomic orbitals in G-space
  type(pwexp_t)                      ::  pwexp_                          !  plane-wave expansion choices
  type(chdens_t)                     ::  chdens_in_                      !  input charge densities
  type(vcomp_t)                      ::  vcomp_in_                       !  input local potential contributions
  type(vcomp_t)                      ::  vcomp_                          !  local potential contributions
  type(strfac_t)                     ::  strfac_                         !  structure factors

  character(len=4)                   ::  flgdal_in                       !  whether the dual approximation is used
  real(REAL64)                       ::  emax_in                         !  kinetic energy cutoff of plane wave expansion (hartree).
  real(REAL64)                       ::  efermi                          !  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  character(len=4)                   ::  author                          !  type of xc wanted (CA=PZ , PW92 , PBE)
  character(len=4)                   ::  xcbase                          !  meta-GGA functional used in xc_mgga
  character(len=4)                   ::  tausrc                          !  source of tau, 'PSI ' for generalized Kohn-Sham
  logical                            ::  lgks                            !  generalized Kohn-Sham meta-GGA (vtau term)

  character(len=60)                  ::  pwline                          !  identifier of the calculation.  May contain miscellaneous information!
  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  title for plots
  character(len=250)                 ::  meta_cpw2000                    !  metadata from cpw2000

  integer            ::  iotape
  integer            ::  iguess                                          !  kept for compatibility
  real(REAL64)       ::  epspsi                                          !  accuracy of eigenvalues
  integer            ::  icmax                                           !  maximum number of iterations for diagonalization
  integer            ::  imtb                                            !  MTB method, 1 projection, 2 matching
  character(len=60)  ::  filename
  character(len=1)   ::  yesno

! parameters

  integer, parameter :: IOHSR = 777                                      !  tape number of hsr2.dat


! hard coded values, change if you know what you are doing...

  flags_%flgpsd = 'PSEUKB'
  iguess = 0
  epspsi = 0.0001
  icmax = 40

! open file and reads data

  filename = 'PW_RHO_V.DAT'
  iotape = 21

  call cpw_pp_band_dos_init(filename, iotape,                            &
       dims_, crys_, spaceg_, flags_, pwexp_, pseudo_, kpoint_,          &
       atorb_, efermi,  author,                                          &
       pwline, title, subtitle ,meta_cpw2000,                            &
       dims_in_, recip_in_, chdens_in_, vcomp_in_, emax_in, flgdal_in)

! prepares calculation

  pwexp_%emax = emax_in

  call cpw_pp_band_prepare(                                              &
     dims_, crys_, spaceg_, recip_, pwexp_, strfac_,  vcomp_,            &
     dims_in_, recip_in_, vcomp_in_)

! generalized Kohn-Sham meta-GGA

  call xc_author_tau(author, xcbase, tausrc)
  lgks = .FALSE.
  if(tausrc == 'PSI ') lgks = .TRUE.


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

! choice of the atomic basis set

  call cpw_pp_atorb_questions(ioreplay,                                  &
      crys_%ntype, crys_%nameat,                                         &
      atorb_%n_bsets, atorb_%norbat, atorb_%lorb,                        &
      atorb_%jbset,                                                      &
      dims_%mxdtyp, dims_%mxdlao, dims_%mxdset)

! modified tight-binding hamiltonian

  call ao_interpolation_prepare(ioreplay, noiData, imtb,                 &
     pwexp_%emax, flags_%flgdal, flags_%flgpsd,                          &
     iguess, epspsi, icmax, pseudo_%ztot,                                &
     crys_%adot, crys_%ntype, crys_%natom, crys_%rat,                    &
     spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                          &
     recip_%ng, recip_%kgv, recip_%phase, recip_%conj,                   &
     recip_%ns, recip_%inds, recip_%kmax,                                &
     recip_%indv, recip_%ek,                                             &
     strfac_%sfact, strfac_%icmplx,                                      &
     vcomp_%veff, lgks, vcomp_%vtau,                                     &
     pseudo_%nqnl, pseudo_%delqnl, pseudo_%vkbraw, pseudo_%nkb,          &
     atorb_%n_bsets, atorb_%jbset,                                       &
     atorb_%norbat, atorb_%nqwf, atorb_%delqwf,                          &
     atorb_%wvfao, atorb_%lorb,                                          &
     dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst,             &
     dims_%mxdlqp, dims_%mxdcub, dims_%mxdlao, dims_%mxdset)

! writes hsr2.dat, header and then the hamiltonian

  write(6,*) '  Writing hsr2.dat for later use.'
  write(6,*)

  open(unit=IOHSR, file='hsr2.dat', form='unformatted', status='replace')

  call ao_int_hsr2_write(IOHSR, meta_cpw2000, title, subtitle,           &
      crys_%adot, pseudo_%ztot, pwexp_%emax,                             &
      crys_%ntype, crys_%natom, crys_%nameat, crys_%rat,                 &
      spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                         &
      atorb_%n_bsets, atorb_%jbset, atorb_%norbat, atorb_%lorb,          &
      atorb_%nqwf, atorb_%delqwf, atorb_%wvfao,                          &
      noiData%lso, noiData%loptical, imtb, noiData%fiData%mp_grid,       &
      noiData%nband,                                                     &
      dims_%mxdtyp, dims_%mxdatm, dims_%mxdlqp, dims_%mxdlao,            &
      dims_%mxdset)

  call NonOrthoInterpWriteToFile(noiData, IOHSR)

  close(unit=IOHSR)

  call NonOrthoInterpFinish(noiData)

  return

end subroutine ao_interpolation_create
