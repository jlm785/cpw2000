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

!>     prints type of calculation

       subroutine tpage_calc(flgcal)


!      input:
!      flgcal     type of calculation

!      written may 29 1987. jlm
!      modified may 16 1990. jlm
!      version 4.1 17 january 1996. jlm
!      version 4.41 19 september 2002. jlm
!      modified 17 august 2004. jlm
!      modified 18 february 2008. jlm
!      Modified titles 14 December 2016,  JLM
!      Removed title, 16 June 2017.  JLM
!      Split in 2, 21 February 2019. JLM
!      Add EPILNG, 21 September 2019.JLM

!      Copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

!      input

       character(len=6), intent(in)      ::  flgcal                      !<  type of calculation

       write(6,*)
       write(6,*)

       if(flgcal == 'ONE   ') then
         write(6,'(//,5x,"Single geometry calculation",//)')
       elseif(flgcal == 'ONEVRD') then
         write(6,'(//,5x,"Single geometry, potential from file",//)')
       elseif(flgcal == 'MICRO ') then
         write(6,'(//,5x,"Micro-canonic Molecular Dynamics ",            &
     &             "simulation",//)')
       elseif(flgcal == 'LANG  ') then
         write(6,'(//,5x,"Langevin Molecular Dynamics simulation",//)')
       elseif(flgcal == 'LBFSYM') then
         write(6,'(//,5x,"Geometry optimization",//)')
       elseif(flgcal == 'VCSLNG') then
         write(6,'(//,5x,"Langevin Molecular Dynamics simulation",       &
     &        " including variational cell shape",//)')
       elseif(flgcal == 'VCSMIC') then
         write(6,'(//,5x,"Micro-canonic Molecular Dynamics simulation",  &
     &        " including variational cell shape",//)')
       elseif(flgcal == 'VCSLBF') then
         write(6,'(//,5x,"Geometry optimization",                        &
     &        " including variational cell shape",//)')
       elseif(flgcal == 'EPILBF') then
         write(6,'(//,5x,"geometry optimization",                        &
     &        " in an epitaxial situation",//)')
       elseif(flgcal == 'EPILNG') then
         write(6,'(//,5x,"Langevin Molecular Dynamics simulation",       &
     &        " including variational cell shape with epitaxial",        &
     &        " constraints",//)')
       elseif(flgcal == 'RSTRT ') then
         write(6,'(//,5x,"restarting old calculation",//)')
       else
         write(6,'(//,5x,"unknown type of calculation, ",                &
     &        "expect disaster",//)')
       endif
       write(6,*)
       write(6,*)

       return
       end subroutine tpage_calc
