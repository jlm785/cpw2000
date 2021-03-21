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

       subroutine cpw_size_alloc_hampsi(lnewnrk,                         &
     &    dims_,crys_,kpoint_,pwexp_,recip_,atorb_,hamallk_,psiallk_)
 
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(kpoint_t)                     ::  kpoint_                    !<  k-point data
       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices
       type(recip_t)                      ::  recip_                     !<  reciprocal space information
       type(atorb_t)                      ::  atorb_                     !<  atomic orbitals in G-space
       type(hamallk_t)                    ::  hamallk_                   !<  hamiltonian size and indexation for all k-points
       type(psiallk_t)                    ::  psiallk_                   !<  psi for all k-points

       logical, intent(in)                ::  lnewnrk                    !<  new number of k-points

       integer               ::  mtxd
       integer               ::  maxdim
       integer               ::  nbaslcao
       logical, save         ::  lfirst = .TRUE.                         !  first time run
       integer               ::  n


       maxdim = 0
       do n=1,kpoint_%nrk
         call size_mtxd(pwexp_%emax,kpoint_%rk(1,n),                     &
&            crys_%adot,recip_%ng,recip_%kgv,mtxd)
         maxdim = max(maxdim,mtxd)
       enddo
  
       if(lfirst) then

         lfirst = .FALSE.       
       
         maxdim = int(1.05*maxdim)

         if(maxdim < 200) then
           dims_%mxdsml = maxdim
         else
           dims_%mxdsml = 3*dims_%mxdbnd+20
         endif

         if(atorb_%latorb) then
           call size_nbaslcao(crys_%ntype,crys_%natom,atorb_%norbat,     &
     &         atorb_%lorb,nbaslcao,dims_%mxdtyp,dims_%mxdlao)
           dims_%mxdsml = max(dims_%mxdsml,nbaslcao)
         endif
         dims_%mxdsml = min(maxdim,dims_%mxdsml)

         dims_%mxddim = maxdim

         allocate(hamallk_%mtxd_allk(dims_%mxdnrk))
         allocate(hamallk_%isort_allk(dims_%mxddim,dims_%mxdnrk))
         allocate(psiallk_%eig_allk(dims_%mxdnrk*dims_%mxdbnd))
         allocate(psiallk_%occ_allk(dims_%mxdnrk*dims_%mxdbnd))
         allocate(psiallk_%psi_allk(dims_%mxddim,dims_%mxdbnd,           &
     &             dims_%mxdnrk))

       else

         if(maxdim > dims_%mxddim .OR. lnewnrk) then

           if(maxdim > dims_%mxddim) then
             dims_%mxddim = maxdim
           endif

           deallocate(hamallk_%mtxd_allk)
           deallocate(hamallk_%isort_allk)
           deallocate(psiallk_%eig_allk)
           deallocate(psiallk_%occ_allk)
           deallocate(psiallk_%psi_allk)

           allocate(hamallk_%mtxd_allk(dims_%mxdnrk))
           allocate(hamallk_%isort_allk(dims_%mxddim,dims_%mxdnrk))
           allocate(psiallk_%eig_allk(dims_%mxdnrk*dims_%mxdbnd))
           allocate(psiallk_%occ_allk(dims_%mxdnrk*dims_%mxdbnd))
           allocate(psiallk_%psi_allk(dims_%mxddim,dims_%mxdbnd,dims_%mxdnrk))

         endif

       endif

       return

       end subroutine cpw_size_alloc_hampsi
