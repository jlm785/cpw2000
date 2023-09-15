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

!>  Reads a file (default PW_RHO_V.DAT set in filename)
!>  with the atomic structure and
!>  the effective potential and charge density.
!>
!>  \author       Jose Luis Martins
!>  \version      5.07
!>  \date         February 2020, 15 September 2023.
!>  \copyright    GNU Public License v2

subroutine cpw_pp_band_dos_init( filename, iotape,                       &
       dims_, spaceg_, flags_, crys_, recip_in_, pseudo_, kpoint_,       &
       pwexp_, chdensin_, vcompin_, atorb_,                              &
       emaxin, efermi, flgdalin, author,                                 &
       pwline, title, subtitle ,meta_cpw2000,                            &
       mxdgvein, mxdnstin)

! written February 1, 2020 from previous code. JLM
! Modified, efermi, 29 November 2021. JLM
! Deallocation, 15 September 2023. JLM

  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  iotape                          !<  tape number for input file filename
  character(len=60), intent(in)      ::  filename                        !<  by default PW_RHO_V.DAT

! output

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(flags_t)                      ::  flags_                          !<  computational flags
  type(spaceg_t)                     ::  spaceg_                         !<  space group information
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(recip_t)                      ::  recip_in_                       !<  reciprocal space information
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(kpoint_t)                     ::  kpoint_                         !<  k-point data
  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space
  type(pwexp_t)                      ::  pwexp_                          !<  plane-wave expansion choices
  type(chdens_t)                     ::  chdensin_                       !<  input charge densities
  type(vcomp_t)                      ::  vcompin_                        !<  local potential contributions

  character(len=3), intent(out)      ::  author                          !<  type of xc wanted (CA=PZ , PW92 , PBE)

  character(len=60), intent(out)     ::  pwline                          !<  identifier of the calculation.  May contain miscellaneous information!
  character(len=50), intent(out)     ::  title                           !<  title for plots
  character(len=140), intent(out)    ::  subtitle                        !<  title for plots
  character(len=250), intent(out)    ::  meta_cpw2000                    !<  metadata from cpw2000

  real(REAL64), intent(out)          ::  emaxin                          !<  kinetic energy cutoff of plane wave expansion (hartree).
  character(len=4), intent(out)      ::  flgdalin                        !<  whether the dual approximation is used
  real(REAL64), intent(out)          ::  efermi                          !<  eigenvalue of highest occupied state (T=0) or fermi energy (T/=0), Hartree

  integer, intent(out)               ::  mxdgvein                        !<  array dimension for input g-space vectors
  integer, intent(out)               ::  mxdnstin                        !<  array dimension for input g-space stars

! Pseudopotential variables not used elsewhere

  character(len=3), allocatable      ::  irel(:)                        !  type of calculation relativistic/spin
  character(len=4), allocatable      ::  icore(:)                       !  type of partial core correction
  character(len=2), allocatable      ::  icorr(:)
  character(len=60), allocatable     ::  iray(:)                        !  information about pseudopotential
  character(len=70), allocatable     ::  ititle(:)                      !  further information about pseudopotential

! other variables

  integer           ::  ipr


! writes preamble to standard output

  write(6,*)
  write(6,'("  Analysis of the electronic structure ")')
  write(6,*)

! open file and reads data

  call pw_rho_v_in_size(filename, iotape,                                &
     dims_%mxdtyp, dims_%mxdatm, mxdgvein, mxdnstin, dims_%mxdlqp,       &
     dims_%mxdlao)

  allocate(crys_%natom(dims_%mxdtyp))
  allocate(crys_%nameat(dims_%mxdtyp))
  allocate(crys_%rat(3,dims_%mxdatm,dims_%mxdtyp))

  allocate(recip_in_%kgv(3,mxdgvein))
  allocate(recip_in_%phase(mxdgvein))
  allocate(recip_in_%conj(mxdgvein))
  allocate(recip_in_%mstar(mxdnstin))
  allocate(chdensin_%den(mxdnstin))
  allocate(chdensin_%dend(mxdnstin))
  allocate(vcompin_%veff(mxdnstin))


  allocate(pseudo_%nq(dims_%mxdtyp))
  allocate(pseudo_%delq(dims_%mxdtyp))
  allocate(pseudo_%vkb(-2:dims_%mxdlqp,0:3,-1:1,dims_%mxdtyp))
  allocate(pseudo_%nkb(0:3,-1:1,dims_%mxdtyp))
  allocate(pseudo_%vloc(-1:dims_%mxdlqp,dims_%mxdtyp))
  allocate(pseudo_%dcor(-1:dims_%mxdlqp,dims_%mxdtyp))
  allocate(pseudo_%dval(-1:dims_%mxdlqp,dims_%mxdtyp))
  allocate(pseudo_%zv(dims_%mxdtyp))
  allocate(atorb_%norbat(dims_%mxdtyp))
  allocate(atorb_%lorb(dims_%mxdlao,dims_%mxdtyp))
  allocate(atorb_%wvfao(-2:dims_%mxdlqp,dims_%mxdlao,dims_%mxdtyp))
  allocate(atorb_%nqwf(dims_%mxdtyp))
  allocate(atorb_%delqwf(dims_%mxdtyp))

  allocate(irel(dims_%mxdtyp))
  allocate(icore(dims_%mxdtyp))
  allocate(icorr(dims_%mxdtyp))
  allocate(iray(dims_%mxdtyp))
  allocate(ititle(dims_%mxdtyp))

  ipr = 1

  call pw_rho_v_in(filename, iotape, ipr,                                &
         pwline, title, subtitle, meta_cpw2000,                          &
         author, flags_%flgscf, flgdalin, emaxin, pwexp_%teleck,         &
         kpoint_%nx, kpoint_%ny, kpoint_%nz,                             &
         kpoint_%sx, kpoint_%sy, kpoint_%sz, pwexp_%nbandin,             &
         crys_%alatt, efermi,                                            &
         crys_%adot, crys_%ntype, crys_%natom,                           &
         crys_%nameat, crys_%rat,                                        &
         spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                      &
         recip_in_%ng, recip_in_%kmax, recip_in_%kgv,                    &
         recip_in_%phase, recip_in_%conj, recip_in_%ns,                  &
         recip_in_%mstar,                                                &
         vcompin_%veff, chdensin_%den, chdensin_%dend,                   &
         irel, icore, icorr, iray, ititle,                               &
         pseudo_%ealraw, pseudo_%zv, pseudo_%ztot,                       &
         pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,             &
         pseudo_%vloc, pseudo_%dcor, pseudo_%dval,                       &
         atorb_%norbat, atorb_%nqwf, atorb_%delqwf, atorb_%wvfao,        &
         atorb_%lorb, atorb_%latorb,                                     &
         dims_%mxdtyp, dims_%mxdatm, mxdgvein, mxdnstin,                 &
         dims_%mxdlqp, dims_%mxdlao)

! processes crystal structure

  ipr = 1
  call print_crystal(ipr, crys_%adot, crys_%ntype, crys_%natom,     &
          crys_%nameat, crys_%rat,                                  &
          dims_%mxdtyp, dims_%mxdatm)


  deallocate(irel)
  deallocate(icore)
  deallocate(icorr)
  deallocate(iray)
  deallocate(ititle)

  return

end subroutine cpw_pp_band_dos_init

