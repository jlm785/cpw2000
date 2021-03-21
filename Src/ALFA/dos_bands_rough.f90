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

!>     Gives a rough estimate of valence band minimum and maximum and
!>     conduction band minimum

       subroutine dos_bands_rough(el,nrk,nband,ztot,ispin,               &
     &       evbb,evbm,ecbm,linsul,                                      &
     &       mxdbnd)

!      Written November 19, 2013. JLM
!      Modified 12 June, 2014. JLM
!      Modified, documentation, 19 September 2020. JLM
!      copyright  J.L.Martins, INESC-MN.

!      version 4.53 of cpw

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)             ::  nrk                           !<  number of k-points
       integer, intent(in)             ::  mxdbnd                        !<  size of number of bands
       real(REAL64), intent(in)        ::  el(mxdbnd,nrk)                !<  eigenvalues in Hartree
       integer, intent(in)             ::  nband(nrk)                    !<  number of bands for each k-points
       real(REAL64), intent(in)        ::  ztot                          !<  total charge density (electrons/cell)
       integer, intent(in)             ::  ispin                         !<  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)

!      output

       real(REAL64), intent(out)       ::  evbb                          !<  energy of bottom of valence band
       real(REAL64), intent(out)       ::  evbm                          !<  energy of maximum of valence band
       real(REAL64), intent(out)       ::  ecbm                          !<  energy of minimum of conduction band
       logical, intent(out)            ::  linsul                        !<  true if it appears to be an insulator
       
!      local variables

       integer  ::  nz2, neig

!      counters

       integer   ::  i

!      constants

       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64


       evbb = el(1,1)
       do i = 1,nrk
         if( evbb > el(1,i) )  evbb = el(1,i) 
       enddo

       neig = nband(1)
       do i = 1,nrk
         if( neig > nband(i) )  neig = nband(i) 
       enddo

!      CHECK ALL THE POSSIBILITIES (METALS!, ZTOT = 1)

       nz2 = nint(ztot/ispin + 0.01)

       if(nz2 >= neig) then
         write(6,*)
         write(6,'("  stopped    number of bands ",i6," smaller than",   &
     &    " half the number of electrons",i6," / 2")') neig,nint(ztot)

         stop

       endif

       evbm = el(nz2,1)
       do i = 1,nrk
         if( evbm < el(nz2,i) )  evbm = el(nz2,i) 
       enddo

       ecbm = el(nz2+1,1)
       do i = 1,nrk
         if( ecbm > el(nz2+1,i) )  ecbm = el(nz2+1,i) 
       enddo

       write(6,*)
       write(6,*) '  Estimates from the k-point sampling '
       write(6,*) '  Energy zero is crystal average potential '
       write(6,*)

       if(ecbm <= evbm) then

         linsul = .FALSE.
         write(6,*)
         write(6,'(" The system is metallic ")')
         write(6,*)
         write(6,'(" An estimate of the bottom of the valence band",     &
     &       " is ", f12.3," eV")') evbb*HARTREE
         write(6,'(" An estimate of the Fermi level is",f12.3,           &
     &        " eV")') 0.5*(evbm+ecbm)*HARTREE
         write(6,*)

       else

         linsul = .TRUE.
         write(6,'(" The system could be an insulator ")')
         write(6,*)
         write(6,'(" An estimate of band gap is ", f12.3," eV")')        &
     &       (ecbm-evbm)*HARTREE
         write(6,*)
         write(6,'(" An estimate of the bottom of the valence band",     &
     &       " is ", f12.3," eV")') evbb*HARTREE
         write(6,'(" An estimate of the valence band maximum is",        &
     &       f12.3," eV")') evbm*HARTREE
         write(6,'(" An estimate of the conduction band minimum",        &
     &       " is  ",f12.3," eV")') ecbm*HARTREE
         write(6,'(" An estimate of the mid-gap level is",               &
     &        f12.3," eV")') 0.5*(evbm+ecbm)*HARTREE
         write(6,*)

       endif

       return
       end subroutine dos_bands_rough
