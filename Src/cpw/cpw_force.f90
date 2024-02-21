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

!>  Manages the calculation of forces and stresses.
!>  it is an interface to force_stress_kb with the extra functions
!>  of adding the Keating corrections and printing the result.
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         20 October 93, 21 February 2024.
!>  \copyright    GNU Public License v2

subroutine cpw_force(iprglob,strxc, ealpha, deltentpy, errfrc,           &
    errstr,                                                              &
    dims_, crys_, total_, ewald_, flags_, spaceg_, recip_, pseudo_,      &
    vcomp_, chdens_, hamallk_, psiallk_, kpoint_, vcsdyn_)

! Written November 2019. JLM
! Modified upstream in January 2010. vff_add_keating. JLM
! Modified, error after keating, 18 February 2020. JLM
! Modified, indentation, another printing choices, 21 February 2024. JLM



  use cpw_variables

  implicit none

  type(dims_t)                       ::  dims_                           !<  array dimensions
  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(enfrst_t)                     ::  total_                          !<  Total energy force stress
  type(enfrst_t)                     ::  ewald_                          !<  Ewald energy force stress
  type(flags_t)                      ::  flags_                          !<  computational flags
  type(spaceg_t)                     ::  spaceg_                         !<  space group information
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential
  type(chdens_t)                     ::  chdens_                         !<  charge densities
  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points
  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points
  type(kpoint_t)                     ::  kpoint_                         !<  k-point data
  type(vcsdyn_t)                     ::  vcsdyn_                         !<  variational cell shape molecular dynamics variables

  integer, intent(in)                ::  iprglob                         !<  controls the amount of printing by subroutines
  real(REAL64), intent(inout)        ::  strxc(3,3)                      !<  contribution of xc to the stress tensor (contravariant,Hartree)
  real(REAL64), intent(in)           ::  ealpha                          !<  G=0 contribution to the total energy (Hartree)

  real(REAL64), intent(out)          ::  deltentpy                       !<  enthalpy difference from last iteration
  real(REAL64), intent(out)          ::  errfrc                          !<  maximum error in force (cartesian coordiantes) Hartree/Bohr
  real(REAL64), intent(out)          ::  errstr                          !<  maximum error in stress

  integer     ::  ipr
  integer     ::  iotape

  integer                ::  minrat                                      !  if =1 minimize with respect to atomic positions
  integer                ::  minstr                                      !  if =1 minimize with respect to all adot variables, if =2 minimize with respect to adot(3,3)

  call force_stress_kb(total_%force, total_%stress, total_%energy,       &
      ewald_%force, ewald_%stress, strxc, ealpha, flags_%flgpsd,         &
      crys_%ntype, crys_%natom, crys_%nameat, crys_%rat, crys_%adot,     &
      spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                         &
      recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%ns,       &
      recip_%mstar, recip_%ek,                                           &
      pseudo_%nq, pseudo_%delq, pseudo_%vkb, pseudo_%nkb,                &
      vcomp_%vion, vcomp_%vhar, vcomp_%vxc, chdens_%den,                 &
      hamallk_%mtxd_allk, hamallk_%isort_allk,                           &
      psiallk_%psi_allk, psiallk_%occ_allk,                              &
      pseudo_%vql, pseudo_%dnc, pseudo_%dvql, pseudo_%ddc,               &
      kpoint_%nrk, kpoint_%nband, kpoint_%rk, kpoint_%wgk,               &
      dims_%mxdtyp, dims_%mxdatm, dims_%mxdlqp, dims_%mxddim,            &
      dims_%mxdbnd, dims_%mxdgve, dims_%mxdnst, dims_%mxdnrk)


!-----------------------clr----------------------------------------


  ipr = 0
  if(iprglob > 1) ipr = 1

  call print_energy(ipr, 'Total',                                        &
      total_%energy, total_%force, total_%stress,                        &
      crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,                &
      dims_%mxdtyp, dims_%mxdatm)

! Adds keating correction

  if(flags_%flgkeat == 'KEATNG') then

    write(6,*)
    write(6,*) '  Using Keating Force Field Corrections'
    write(6,*)

    ipr = 0
    if(iprglob > 2) ipr = 1
    iotape = 6

    call vff_add_keating(ipr, iotape,                                    &
        crys_%ntype, crys_%natom, crys_%nameat,                          &
        crys_%rat, crys_%adot,                                           &
        total_%force, total_%energy, total_%stress,                      &
        dims_%mxdtyp, dims_%mxdatm)

    ipr = 1
    call print_energy(ipr, 'Total+Keating',                              &
        total_%energy, total_%force, total_%stress,                      &
        crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,              &
        dims_%mxdtyp, dims_%mxdatm)

  else

    ipr = 0
    if(iprglob < 2) ipr = 1

    call print_energy(ipr, 'Total',                                      &
        total_%energy, total_%force, total_%stress,                      &
        crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,              &
        dims_%mxdtyp, dims_%mxdatm)

  endif

  if(flags_%flgcal == 'LBFSYM' .or. flags_%flgcal == 'VCSLBF' .or.       &
           flags_%flgcal == 'EPILBF') then

    minrat = 1
    if(flags_%flgcal == 'LBFSYM') then
      minstr = 0
    elseif(flags_%flgcal == 'VCSLBF') then
      minstr = 1
    elseif(flags_%flgcal == 'EPILBF') then
      minstr = 2
    endif

    call force_stress_error(total_%energy, total_%force,                 &
        total_%stress, vcsdyn_%press, vcsdyn_%strext,                    &
        deltentpy, errfrc, errstr, minrat, minstr,                       &
        crys_%adot, crys_%ntype, crys_%natom,                            &
        dims_%mxdtyp, dims_%mxdatm)

  endif

  return

end subroutine cpw_force
