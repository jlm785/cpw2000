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

       subroutine cpw_init_mov_print_sym(ist0,iprglob,symkip,symtol,     &
     &       crys_,flags_,dims_,spaceg_,pwexp_,acc_,xc_,moldyn_,vcsdyn_)
 
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(moldyn_t)                     ::  moldyn_                    !<  molecular dynamics variables
       type(vcsdyn_t)                     ::  vcsdyn_                    !<  variational cell shape molecular dynamics variables
       type(flags_t)                      ::  flags_                     !<  computational flags
       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices
       type(xc_t)                         ::  xc_                        !<  exchange and correlation choice
       type(acc_t)                        ::  acc_                       !<  accuracy parameters
       type(spaceg_t)                     ::  spaceg_                    !<  space group information

       integer, intent(inout)             ::  ist0                       !<  0 or last available dynamics step
       integer, intent(in)                ::  iprglob                    !<  controls the amount of printing by subroutines
       logical, intent(in)                ::  symkip                     !<  calculate symmetry
       real(REAL64), intent(in)           ::  symtol                     !<  tolerance for symmetry recognition

       integer              ::  ipr
       real(REAL64)         ::  tol
       integer              ::  isym

!      Random initial velocities

       call move_random_v_atom(moldyn_%iseed,moldyn_%tempinik,moldyn_%vat,vcsdyn_%vadot,                 &
     & crys_%ntype,crys_%natom,crys_%atmass,crys_%adot,                                          &
     & dims_%mxdtyp,dims_%mxdatm)

       if(flags_%flgcal == 'RSTRT ') then

         call move_restart_in(flags_%flgcal,                                    &
     &   crys_%ntype,crys_%natom,crys_%nameat,crys_%atmass,crys_%rat,moldyn_%vat,crys_%adot,vcsdyn_%vadot,                   &
     &   ist0,moldyn_%tstep,moldyn_%beta,moldyn_%tempk,moldyn_%iseed,vcsdyn_%strext,vcsdyn_%press,vcsdyn_%celmas,                &
     &   moldyn_%rat1,moldyn_%frc1,vcsdyn_%adot1,vcsdyn_%frcel1,                                         &
     &   dims_%mxdatm,dims_%mxdtyp)

       endif

       call print_parameters(flags_%flgcal,flags_%flgdal,flags_%flgscf,                       &
     & pwexp_%emax,xc_%author,xc_%tblaha,acc_%epscv,acc_%epscvao,acc_%epspsi,                          &
     & moldyn_%tempk,pwexp_%teleck,moldyn_%tempinik,moldyn_%nstep,moldyn_%tstep,moldyn_%beta,moldyn_%iseed,                     &
     & vcsdyn_%press,vcsdyn_%strext,vcsdyn_%celmas,pwexp_%lkplusg,pwexp_%epskplusg)


!      prints the geometry

       ipr = 1
       call print_crystal(ipr,crys_%adot,crys_%ntype,crys_%natom,crys_%nameat,crys_%rat,               &
     & dims_%mxdtyp,dims_%mxdatm)


!      calculates the symmetry operations

       isym = 0
       if(symkip) isym = 1
       ipr = 0
       tol = symtol
       if(iprglob .gt. 0) ipr = 1

       call sym_identify(isym,ipr,tol,                                   &
     & spaceg_%ntrans,spaceg_%mtrx,spaceg_%tnp,                          &
     & crys_%adot,crys_%ntype,crys_%natom,crys_%rat,                     &
     & dims_%mxdtyp,dims_%mxdatm)


       return

       end subroutine cpw_init_mov_print_sym
