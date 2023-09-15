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

!>  Deallocates stuff to help debugging (with valgrind for example)
!>
!>  \author       Jose Luis Martins
!>  \version      5.07
!>  \date         February 2020, 29 November 2021.
!>  \copyright    GNU Public License v2

  subroutine cpw_pp_band_dos_clean(crys_, recip_in_, pseudo_,            &
        chdensin_, vcompin_, atorb_)

! written 15 September 2023. JLM

  use cpw_variables

  implicit none

!  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  type(crys_t)                       ::  crys_                           !<  crystal structure
  type(recip_t)                      ::  recip_in_                       !<  reciprocal space information
  type(pseudo_t)                     ::  pseudo_                         !<  pseudo-potential (Kleinman-Bylander)
  type(atorb_t)                      ::  atorb_                          !<  atomic orbitals in G-space
  type(chdens_t)                     ::  chdensin_                       !<  input charge densities
  type(vcomp_t)                      ::  vcompin_                        !<  local potential contributions


  deallocate(crys_%natom)
  deallocate(crys_%nameat)
  deallocate(crys_%rat)

  deallocate(recip_in_%kgv)
  deallocate(recip_in_%phase)
  deallocate(recip_in_%conj)
  deallocate(recip_in_%mstar)
  deallocate(chdensin_%den)
  deallocate(chdensin_%dend)
  deallocate(vcompin_%veff)


  deallocate(pseudo_%nq)
  deallocate(pseudo_%delq)
  deallocate(pseudo_%vkb)
  deallocate(pseudo_%nkb)
  deallocate(pseudo_%vloc)
  deallocate(pseudo_%dcor)
  deallocate(pseudo_%dval)
  deallocate(pseudo_%zv)
  deallocate(atorb_%norbat)
  deallocate(atorb_%lorb)
  deallocate(atorb_%wvfao)
  deallocate(atorb_%nqwf)
  deallocate(atorb_%delqwf)

  return

end subroutine cpw_pp_band_dos_clean

