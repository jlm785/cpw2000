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

       subroutine cpw_bzint(iprglob,lnewnrk,                             &
     &    dims_,kpoint_,crys_,spaceg_,pwexp_)
 
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(kpoint_t)                     ::  kpoint_                    !<  k-point data
       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices
       type(spaceg_t)                     ::  spaceg_                    !<  space group information

       integer, intent(in)                ::  iprglob                    !<  controls the amount of printing by subroutines

       logical, intent(out)               ::  lnewnrk                    !<  new number of k-points

       integer      ::  ipr

       logical, save         ::  lfirst = .TRUE.                         !  first time run
       integer               ::  maxnrk                                  !  provisional array dimensions

       call size_mxdnrk(kpoint_%nx,kpoint_%ny,kpoint_%nz,kpoint_%sx,     &
     &    kpoint_%sy,kpoint_%sz,crys_%adot,spaceg_%ntrans,spaceg_%mtrx,  &
     &    maxnrk)
     
!      mxdpnt = kpoint_%nx*kpoint_%ny*kpoint_%nz

!      nbandin and kpoint_%nx*kpoint_%ny*kpoint_%nz do not change.

       if(lfirst) then

         dims_%mxdnrk = maxnrk

         lnewnrk = .TRUE.
         lfirst = .FALSE.

         allocate(kpoint_%rk(3,dims_%mxdnrk))
         allocate(kpoint_%wgk(dims_%mxdnrk))
         allocate(kpoint_%nband(dims_%mxdnrk))
         allocate(kpoint_%indk(6,dims_%mxdnrk))
         allocate(kpoint_%kmap(3,kpoint_%nx,kpoint_%ny,kpoint_%nz))

       else

         lnewnrk = .FALSE.

         if(maxnrk > dims_%mxdnrk) then
         
           dims_%mxdnrk = maxnrk

           lnewnrk = .TRUE.

           deallocate(kpoint_%rk)
           deallocate(kpoint_%wgk)
           deallocate(kpoint_%nband)
           deallocate(kpoint_%indk)
           deallocate(kpoint_%kmap)

           allocate(kpoint_%rk(3,dims_%mxdnrk))
           allocate(kpoint_%wgk(dims_%mxdnrk))
           allocate(kpoint_%nband(dims_%mxdnrk))
           allocate(kpoint_%indk(6,dims_%mxdnrk))
           allocate(kpoint_%kmap(3,kpoint_%nx,kpoint_%ny,kpoint_%nz))

         endif

       endif

       ipr = 0
       if(iprglob > 0) ipr = 1

       call int_pnt(pwexp_%nbandin,kpoint_%nx,kpoint_%ny,kpoint_%nz,     &
     &    kpoint_%sx,kpoint_%sy,kpoint_%sz,ipr,                          &
     &    crys_%adot,                                                    &
     &    spaceg_%ntrans,spaceg_%mtrx,                                   &
     &    kpoint_%nrk,kpoint_%rk,kpoint_%wgk,kpoint_%nband,              &
     &    kpoint_%indk, kpoint_%kmap,                                    &
     &    dims_%mxdnrk,dims_%mxdbnd)

       return

       end subroutine cpw_bzint
