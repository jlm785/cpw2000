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

       subroutine cpw_move(newcel, newcalc, iconv, istep ,mxdlbf,        &
     &     nsave, chdsave,                                               &
     &     dims_, crys_, flags_, vcsdyn_, total_, moldyn_, spaceg_,      &
     &     chdens_, recip_)
     
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(flags_t)                      ::  flags_                     !<  computational flags
       type(vcsdyn_t)                     ::  vcsdyn_                    !<  variational cell shape molecular dynamics variables
       type(enfrst_t)                     ::  total_                     !<  Total energy force stress
       type(moldyn_t)                     ::  moldyn_                    !<  molecular dynamics variables
       type(spaceg_t)                     ::  spaceg_                    !<  space group information
       type(chdens_t)                     ::  chdens_                    !<  charge densities    
       type(recip_t)                      ::  recip_                     !<  reciprocal space information

       integer,intent(in)                 ::  istep                      !<  md step. Equal to 1 in first step of molecular dynamics
       logical, intent(out)               ::  newcel                     !<  indicates that adot has changed on output
       logical, intent(in)                ::  newcalc                    !<  indicates that it is a new calculation (equivalent iter = ist0+1)
       integer, intent(out)               ::  iconv                      !<  iconv = 1, lbfgs converged of flgcal = one

       integer, intent(in)                ::  mxdlbf                     !<  number of iterations remembered by lbfgs

       integer, intent(in)                ::  nsave(3)                   !<  dimensions of chdsave
       complex(REAL64), intent(out)       ::                             &
     & chdsave(-nsave(1):nsave(1),-nsave(2):nsave(2),-nsave(3):nsave(3)) !<  quantity in reciprocal point i,j,k

       call move_atom_cell(flags_%flgcal, newcel, iconv, istep,          &
     &   crys_%rat, moldyn_%vat, crys_%adot, vcsdyn_%vadot,              &
     &   total_%energy, total_%force, total_%stress,                     &
     &   moldyn_%tstep, moldyn_%ekin, vcsdyn_%ekcell, vcsdyn_%strext,    &
     &   vcsdyn_%press, vcsdyn_%celmas,                                  &
     &   moldyn_%beta, moldyn_%tempk, moldyn_%iseed,                     &
     &   mxdlbf, moldyn_%pgtol, moldyn_%dxmax,                           &
     &   moldyn_%rat1, moldyn_%frc1, vcsdyn_%adot1, vcsdyn_%frcel1,      &
     &   crys_%ntype, crys_%natom, crys_%nameat, crys_%atmass,           &
     &   spaceg_%ntrans, spaceg_%mtrx, spaceg_%tnp,                      &
     &   dims_%mxdtyp, dims_%mxdatm)


       call move_save_density(newcalc, flags_%flgcal, nsave, chdsave,    &
     &   chdens_%den, chdens_%dens, chdens_%dend, chdens_%dend1,         &
     &   recip_%ng, recip_%kgv, recip_%phase, recip_%conj, recip_%inds,  &
     &   recip_%ns, recip_%mstar,                                        &
     &   dims_%mxdgve, dims_%mxdnst)


!      Allows restart of the calculation

       call move_restart_out(flags_%flgcal,                              &
     &   crys_%ntype, crys_%natom, crys_%nameat, crys_%atmass,           &
     &   crys_%rat, moldyn_%vat, crys_%adot, vcsdyn_%vadot,              &
     &   istep, moldyn_%tstep, moldyn_%beta, moldyn_%tempk,              &
     &   moldyn_%iseed, vcsdyn_%strext, vcsdyn_%press, vcsdyn_%celmas,   &
     &   moldyn_%rat1, moldyn_%frc1, vcsdyn_%adot1, vcsdyn_%frcel1,      &
     &   dims_%mxdatm, dims_%mxdtyp)

        return

        end subroutine cpw_move
