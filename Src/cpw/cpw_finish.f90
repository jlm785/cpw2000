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

       subroutine cpw_finish(fname,iotape,meta_pwdat,meta_cpw2000,       &
     &     dims_,crys_,xc_,flags_,pwexp_,kpoint_,recip_,vcomp_,chdens_)
 
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(xc_t)                         ::  xc_                        !<  exchange and correlation choice
       type(flags_t)                      ::  flags_                     !<  computational flags
       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices
       type(kpoint_t)                     ::  kpoint_                    !<  k-point data
       type(recip_t)                      ::  recip_                     !<  reciprocal space information
       type(vcomp_t)                      ::  vcomp_                     !<  Componemts of local potential
       type(chdens_t)                     ::  chdens_                    !<  charge densities    

       character(len=*),intent(in)        ::  fname                      !<  file name, default PW_RHO_V.DAT
       integer, intent(in)                ::  iotape                     !<  tape number:-)
       character(len=250), intent(in)     ::  meta_pwdat                 !<  metadata from PW.DAT
       character(len=250), intent(in)     ::  meta_cpw2000               !<  metadata from cpw2000

       call pw_rho_v_out(fname, iotape, xc_%author, xc_%tblaha,          &
     &    flags_%flgscf, flags_%flgdal,                                  &
     &    meta_pwdat,meta_cpw2000,                                       &
     &    pwexp_%emax, pwexp_%teleck, kpoint_%nx, kpoint_%ny,            &
     &    kpoint_%nz, kpoint_%sx, kpoint_%sy, kpoint_%sz,                &
     &    pwexp_%nbandin, crys_%alatt, crys_%adot,                       &
     &    crys_%ntype,crys_%natom,crys_%nameat,crys_%rat,                &
     &    recip_%ng, recip_%kmax, recip_%kgv, recip_%phase,              &
     &    recip_%conj, recip_%ns, recip_%mstar,                          &
     &    vcomp_%veff, chdens_%den, chdens_%dens,                        &
     &    dims_%mxdtyp, dims_%mxdatm, dims_%mxdgve, dims_%mxdnst)
     
       call write_cpwout(meta_pwdat,                                     &
     &    crys_%adot, crys_%ntype, crys_%natom, crys_%nameat,            &
     &    crys_%rat, crys_%atmass, crys_%alatt,                          &
     &    pwexp_%emax, pwexp_%nbandin, kpoint_%nx, kpoint_%ny,           &
     &    kpoint_%nz, kpoint_%sx, kpoint_%sy, kpoint_%sz,                &
     &    .TRUE., .TRUE.,                                                &
     &    dims_%mxdtyp, dims_%mxdatm)

       end subroutine cpw_finish
