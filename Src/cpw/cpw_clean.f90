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
!>  Deallocates arrays before exiting cpw
!>
!>  \author       Jose Luis Martins
!>  \version      5.13
!>  \date         18 February 2026, 16 October 2026.
!>  \copyright    GNU Public License v2

subroutine cpw_clean(crys_, moldyn_, recip_, strfac_, chdens_, psding_,  &
      vcomp_, pseudo_, atorb_, total_, ewald_, kpoint_,                  &
      hamallk_, psiallk_)


! Written 11 February 2026. JLM
! Added new_atorb_ temporarily (only used in read pseudo). 18 February 2026.
! core kinetic energy densidty. 3 March 2026. JLM
! vcomp_%vtau, vtau for the generalized Kohn-Sham meta-GGA. 6 October 2026. JLM+claude
! pseudo_%tnc, pseudo_%dtauc. 6 October 2026. JLM+claude
! psding_ (pseudopotential quantities in G-space for the crystal). 14 October 2026. JLM+claude WARNING NEW API
! pseudo_%nq, pseudo_%delq, pseudo_%vkb renamed nqnl, delqnl, vkbraw. 14 October 2026. JLM+claude
! new_atorb_ removed, atorb_ has several basis sets. 16 October 2026. JLM+claude


  use cpw_variables

  implicit none

  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(moldyn_t)                     ::  moldyn_                         !<  molecular dynamics variables
  type(recip_t)                      ::  recip_                          !<  reciprocal space information
  type(strfac_t)                     ::  strfac_                         !<  structure factors
  type(chdens_t)                     ::  chdens_                         !<  charge densities
  type(vcomp_t)                      ::  vcomp_                          !<  Componemts of local potential
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(psding_t)                     ::  psding_                         !<  pseudopotential quantities in G-space for the crystal
  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space
  type(enfrst_t)                     ::  total_                          !<  Total energy force stress
  type(enfrst_t)                     ::  ewald_                          !<  Ewald energy force stress
  type(kpoint_t)                     ::  kpoint_                         !<  k-point data
  type(hamallk_t)                    ::  hamallk_                        !<  hamiltonian size and indexation for all k-points
  type(psiallk_t)                    ::  psiallk_                        !<  psi for all k-points


  deallocate(crys_%natom)
  deallocate(crys_%rat)
  deallocate(crys_%atmass)
  deallocate(crys_%nameat)

  deallocate(moldyn_%vat)
  deallocate(moldyn_%rat1)
  deallocate(moldyn_%frc1)



  deallocate(total_%force)

  deallocate(ewald_%force)

  deallocate(pseudo_%nqnl)
  deallocate(pseudo_%delqnl)
  deallocate(pseudo_%zv)
  deallocate(pseudo_%nkb)
  deallocate(pseudo_%vloc)
  deallocate(pseudo_%dcor)
  deallocate(pseudo_%dval)
  deallocate(pseudo_%tauc_q)

  deallocate(pseudo_%vkbraw)

  deallocate(atorb_%n_bsets)
  deallocate(atorb_%ibset)
  deallocate(atorb_%norbat)
  deallocate(atorb_%nqwf)
  deallocate(atorb_%delqwf)
  deallocate(atorb_%wvfao)
  deallocate(atorb_%lorb)

  deallocate(kpoint_%rk)
  deallocate(kpoint_%wgk)
  deallocate(kpoint_%nband)
  deallocate(kpoint_%indk)
  deallocate(kpoint_%kmap)

  deallocate(hamallk_%mtxd_allk)
  deallocate(hamallk_%isort_allk)

  deallocate(psiallk_%eig_allk)
  deallocate(psiallk_%occ_allk)
  deallocate(psiallk_%psi_allk)

  deallocate(recip_%kgv)
  deallocate(recip_%inds)
  deallocate(recip_%phase)
  deallocate(recip_%conj)
  deallocate(recip_%indv)
  deallocate(recip_%mstar)
  deallocate(recip_%izstar)
  deallocate(recip_%ek)

  deallocate(strfac_%sfact)

  deallocate(psding_%vql)
  deallocate(psding_%dnc)
  deallocate(psding_%dvql)
  deallocate(psding_%ddc)
  deallocate(psding_%tnc)
  deallocate(psding_%dtauc)

  deallocate(chdens_%den)
  deallocate(psding_%denc)
  deallocate(chdens_%dens)
  deallocate(chdens_%dend)
  deallocate(chdens_%dend1)
  deallocate(psding_%tauc_g)

  deallocate(vcomp_%vion)
  deallocate(vcomp_%vhar)
  deallocate(vcomp_%vxc)
  deallocate(vcomp_%vtau)
  deallocate(vcomp_%veff)


  return

end subroutine cpw_clean
