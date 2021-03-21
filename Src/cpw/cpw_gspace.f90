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

       subroutine cpw_gspace(iprglob,kmscr,                              &
     &    dims_,crys_,spaceg_,pwexp_,recip_,strfac_,pseudo_,chdens_,     &
     &    vcomp_,flags_)
 
       use cpw_variables

       implicit none

       type(dims_t)                       ::  dims_                      !<  array dimensions
       type(crys_t)                       ::  crys_                      !<  crystal structure
       type(flags_t)                      ::  flags_                     !<  computational flags
       type(pwexp_t)                      ::  pwexp_                     !<  plane-wave expansion choices
       type(recip_t)                      ::  recip_                     !<  reciprocal space information
       type(strfac_t)                     ::  strfac_                    !<  structure factors
       type(chdens_t)                     ::  chdens_                    !<  charge densities    
       type(vcomp_t)                      ::  vcomp_                     !<  Componemts of local potential
       type(pseudo_t)                     ::  pseudo_                    !<  pseudo-potential (Kleinman-Bylander)
       type(spaceg_t)                     ::  spaceg_                    !<  space group information

       integer, intent(in)                ::  iprglob                    !<  controls the amount of printing by subroutines

       integer      ::  ipr
       integer, intent(out)               ::  kmscr(7)                   !<  max value of kgv(i,n) used for the potential FFT mesh (DUAL APPROXIMATION TYPE)

       logical, save         ::  lfirst = .TRUE.                         !  first time run
       integer               ::  maxgve,maxnst,maxcub                    !  provisional array dimensions

       call size_g_space(pwexp_%emax,crys_%adot,                         &
     &     spaceg_%ntrans,spaceg_%mtrx,                                  &
     &     maxgve,maxnst,maxcub)
  
       if(lfirst) then

         lfirst = .FALSE.

         dims_%mxdgve = maxgve
         dims_%mxdnst = maxnst
         dims_%mxdcub = maxcub

         allocate(recip_%kgv(3,dims_%mxdgve))
         allocate(recip_%inds(dims_%mxdgve))
         allocate(recip_%phase(dims_%mxdgve))
         allocate(recip_%conj(dims_%mxdgve))
         allocate(recip_%indv(dims_%mxdcub))
         allocate(recip_%mstar(dims_%mxdnst))
         allocate(recip_%izstar(dims_%mxdnst))
         allocate(recip_%ek(dims_%mxdnst))

         allocate(strfac_%sfact(dims_%mxdtyp,dims_%mxdnst))

         allocate(pseudo_%vql(dims_%mxdtyp,dims_%mxdnst))
         allocate(pseudo_%dnc(dims_%mxdtyp,dims_%mxdnst))
         allocate(pseudo_%dvql(dims_%mxdnst))
         allocate(pseudo_%ddc(dims_%mxdnst))

         allocate(chdens_%den(dims_%mxdnst))
         allocate(chdens_%denc(dims_%mxdnst))
         allocate(chdens_%dens(dims_%mxdnst))
         allocate(chdens_%dend(dims_%mxdnst))
         allocate(chdens_%dend1(dims_%mxdnst))

         allocate(vcomp_%vion(dims_%mxdnst))
         allocate(vcomp_%vhar(dims_%mxdnst))
         allocate(vcomp_%vxc(dims_%mxdnst))
         allocate(vcomp_%veff(dims_%mxdnst))

       else

         if(maxgve > dims_%mxdgve) then

           dims_%mxdgve = maxgve
         
           deallocate(recip_%kgv)
           deallocate(recip_%inds)
           deallocate(recip_%phase)
           deallocate(recip_%conj)
           
           allocate(recip_%kgv(3,dims_%mxdgve))
           allocate(recip_%inds(dims_%mxdgve))
           allocate(recip_%phase(dims_%mxdgve))
           allocate(recip_%conj(dims_%mxdgve))


         endif

         if(maxcub > dims_%mxdcub) then

           dims_%mxdcub = maxcub

           deallocate(recip_%indv)
           allocate(recip_%indv(dims_%mxdcub))

         endif
         
         if(maxnst > dims_%mxdnst) then

           dims_%mxdnst = maxnst

           deallocate(recip_%mstar)
           deallocate(recip_%izstar)
           deallocate(recip_%ek)


           deallocate(strfac_%sfact)

           deallocate(pseudo_%vql)
           deallocate(pseudo_%dnc)
           deallocate(pseudo_%dvql)
           deallocate(pseudo_%ddc)

           deallocate(chdens_%den)
           deallocate(chdens_%denc)
           deallocate(chdens_%dens)
           deallocate(chdens_%dend)
           deallocate(chdens_%dend1)

           deallocate(vcomp_%vion)
           deallocate(vcomp_%vhar)
           deallocate(vcomp_%vxc)
           deallocate(vcomp_%veff)

           allocate(recip_%mstar(dims_%mxdnst))
           allocate(recip_%izstar(dims_%mxdnst))
           allocate(recip_%ek(dims_%mxdnst))

           allocate(strfac_%sfact(dims_%mxdtyp,dims_%mxdnst))

           allocate(pseudo_%vql(dims_%mxdtyp,dims_%mxdnst))
           allocate(pseudo_%dnc(dims_%mxdtyp,dims_%mxdnst))
           allocate(pseudo_%dvql(dims_%mxdnst))
           allocate(pseudo_%ddc(dims_%mxdnst))

           allocate(chdens_%den(dims_%mxdnst))
           allocate(chdens_%denc(dims_%mxdnst))
           allocate(chdens_%dens(dims_%mxdnst))
           allocate(chdens_%dend(dims_%mxdnst))
           allocate(chdens_%dend1(dims_%mxdnst))

           allocate(vcomp_%vion(dims_%mxdnst))
           allocate(vcomp_%vhar(dims_%mxdnst))
           allocate(vcomp_%vxc(dims_%mxdnst))
           allocate(vcomp_%veff(dims_%mxdnst))

         endif

       endif

       ipr = 0
       if(iprglob > 0) ipr = 1
       if(iprglob == 4) ipr = 2

       call g_space(ipr,pwexp_%emax,                                     &
     &     crys_%adot,spaceg_%ntrans,spaceg_%mtrx,spaceg_%tnp,           &
     &     recip_%ng,recip_%kgv,recip_%phase,recip_%conj,                &
     &     recip_%inds,recip_%kmax,recip_%indv,recip_%ns,recip_%mstar,   &
     &     recip_%ek,recip_%izstar,                                      &
     &     dims_%mxdgve,dims_%mxdnst,dims_%mxdcub)



!      gspace size for dual space method

       if(flags_%flgdal == 'DUAL') then
         kmscr(1) = recip_%kmax(1)/2 + 2
         kmscr(2) = recip_%kmax(2)/2 + 2
         kmscr(3) = recip_%kmax(3)/2 + 2
       else
         kmscr(1) = recip_%kmax(1)
         kmscr(2) = recip_%kmax(2)
         kmscr(3) = recip_%kmax(3)
       endif

       return

       end subroutine cpw_gspace
