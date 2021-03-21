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

       subroutine cpw_print_crystal_energy(iprglob, iconv, istep,        &
     &    elects,                                                        &
     &    dims_, crys_, flags_, total_, moldyn_, vcsdyn_)
     
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(flags_t)                      ::  flags_                     !<  computational flags
       type(enfrst_t)                     ::  total_                     !<  Total energy force stress
       type(moldyn_t)                     ::  moldyn_                    !<  molecular dynamics variables
       type(vcsdyn_t)                     ::  vcsdyn_                    !<  variational cell shape molecular dynamics variables

       integer, intent(in)                ::  iprglob                    !<  controls the amount of printing by subroutines
       integer, intent(in)                ::  iconv                      !<  iconv = 1, lbfgs converged of flgcal = one
       integer,intent(in)                 ::  istep                      !<  md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(in)           ::  elects                     !<  electronic temperature*entropy (hartree)

       integer     ::  ipr

       if(iconv /= 1 .AND. istep /= moldyn_%nstep) then

!        prints the geometry

         ipr = 0
         if(iprglob > 0) ipr = 1

         call print_crystal(ipr, crys_%adot, crys_%ntype, crys_%natom,   &
     &      crys_%nameat, crys_%rat,                                     &
     &      dims_%mxdtyp, dims_%mxdatm)

       endif

!      Print energies

       call move_print_mdenergy(flags_%flgcal, istep,                    &
     &    total_%energy, moldyn_%ekin, vcsdyn_%ekcell, elects)

        return

        end subroutine cpw_print_crystal_energy

