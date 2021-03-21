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

!>     This subroutine suggests a reference energy for a band structure plot

       subroutine out_band_eref(neig,nrk,ztot,ispin,ivc,e_of_k,          &
     &                          eref,nocc)

!      Written September 4, 2014. JLM
!      Modified, writes information, 7 November 2018. JLM
!      Modified, documentation, August 1, 2019. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      version 4.93

       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  neig                       !<  number of eigenvectors (without spin)
       integer, intent(in)                ::  nrk                        !<  number of k-points in path
       real(REAL64), intent(in)           ::  ztot                       !<  total charge density (electrons/cell)
       integer, intent(in)                ::  ispin                      !<  spin degeneracy (must be 1 or 2)
       integer, intent(in)                ::  ivc                        !<  1=CBM, 2=VBM, 3=midgap=(CBM+VBM)/2
       
       real(REAL64), intent(in)           ::  e_of_k(2*neig/ispin,nrk)   !<  band energies of k-point in plot

!      output

       real(REAL64), intent(out)          ::  eref                       !<  reference energy for plot
       integer, intent(out)               ::  nocc                       !<  number of occupied bands (in semiconductors)
       
!      allocatable arrays

       integer, allocatable               ::  indx(:)                    !  sorting index


!      local variables

       real(REAL64)          ::  evbm           !  energy of valence band maximum
       real(REAL64)          ::  ecbm           !  energy of conduction band minimum
       logical               ::  lmetal         !  metallic case
       character(len=16)     ::  cso            !  for spin-orbit

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64
       real(REAL64), parameter ::  EV = 27.2116_REAL64

!      counters

       integer    ::  j, n, k, irk


       if(ispin /=1 .and. ispin /= 2) then

         eref = ZERO
         nocc = 0

       else

         if(ispin == 2) then
           cso = ' no spin-orbit  '
         else
           cso = ' with spin-orbit'
         endif

         lmetal = .FALSE.
         if(ispin == 2 .and. mod(nint(ztot),2) == 1) lmetal = .TRUE.
         
         if(.not. lmetal) then

           n = min(nint(ztot/ispin + 0.01),2*neig/ispin-1)
           
           if(n == 0) then
           
             evbm = e_of_k(1,1)

             if(nrk > 1) then

               do irk = 2,nrk
           
                 evbm = max(e_of_k(1,irk),evbm)

               enddo

             endif

             ecbm = evbm
             
           else

             allocate(indx(2*neig/ispin))

             call sort(2*neig/ispin,e_of_k(1,1),indx)
           
             evbm = e_of_k(indx(n),1)
             ecbm = e_of_k(indx(n+1),1)

             if(nrk > 1) then

               do irk = 2,nrk

                 call sort(2*neig/ispin,e_of_k(1,irk),indx)
           
                 evbm = max(e_of_k(indx(n),irk),evbm)
                 ecbm = min(e_of_k(indx(n+1),irk),ecbm)

               enddo

             endif

             deallocate(indx)

           endif

           if(ecbm > evbm) then

             if(ivc == 1) then
               eref = evbm
             elseif(ivc == 2) then
               eref = ecbm
             else
               eref = (evbm+ecbm)/2
             endif
             nocc = n

             write(6,*)
             write(6,'(i8,"   Occupied bands in apparent ",              &
     &         "semiconductor/insulator",a16)') nocc,cso
             write(6,'(f12.6,"   shift applied to bands (eV) ",a16)')    &
     &             -eref*EV,cso
             write(6,'(2f12.6,"   valence band maximum and conduction",  &
     &         " band minimum (eV) ",a16)') evbm*EV,ecbm*EV,cso
             write(6,*)

           else
           
             lmetal = .TRUE.

           endif
         endif
         
         if(lmetal) then

           n = min(nint((nrk*ztot)/ispin + 0.01),2*neig*nrk/ispin)

           allocate(indx(2*neig*nrk/ispin))

           call sort(2*neig*nrk/ispin,e_of_k,indx)

           k = (indx(n)-1)/(2*neig/ispin) + 1
           j = mod((indx(n)-1),(2*neig/ispin)) + 1
           eref = e_of_k(j,k)
           nocc = 0

           write(6,*)
           write(6,'(f12.6,"   E_F, estimate of Fermi energy (eV) ",     &
     &         a16)') eref*EV,cso
           write(6,'(f12.6,"   shift applied to bands (eV) ",a16)')      &
     &             -eref*EV,cso
           write(6,*)

           deallocate(indx)

         endif

       endif

       return

       end subroutine out_band_eref
