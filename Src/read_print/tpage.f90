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

!>     prints version, time and date

       subroutine tpage(vdriv)

!      written may 29 1987. jlm
!      modified may 16 1990. jlm
!      version 4.1 17 january 1996. jlm
!      version 4.41 19 september 2002. jlm
!      modified 17 august 2004. jlm
!      modified 18 february 2008. jlm
!      Modified titles 14 December 2016,  JLM
!      Removed title, 16 June 2017.  JLM
!      Split in 2, 21 February 2019. JLM

!      copyright inesc-mn/Jose Luis Martins

!      version 4.94

       implicit none

!      input

       character(len=4), intent(in)      ::  vdriv                       !<  version of the calling program

!      local variables

       character(len=4)    ::  vers
       logical             ::  ldevel                                    !  development branch, minor version may be incompatible
       character(len=9)    ::  bdate
       character(len=8)    ::  btime


       call version(vers,ldevel)

       if(ldevel) then
         if(vdriv /= vers) then
           write(6,'("  *** STOPPED in tpage")')
           write(6,'("  library version ",a4,"  main program version ",  &
     &               a4)') vers,vdriv
 
           stop

         endif
       else
         if(vdriv(1:3) /= vers(1:3)) then
           write(6,'("  *** STOPPED in tpage")')
           write(6,'("  library version ",a4,"  main program version ",  &
     &               a4)') vers,vdriv
 
           stop

         endif

         if(vdriv(4:4) /= vers(4:4)) then
           write(6,*)
           write(6,'("     WARNING    WARNING    WARNING   in tpage")') 
           write(6,*)
           write(6,'("  library version: ",a4," cpw driver version: ",   &
     &                a4)') vdriv,vers
           write(6,*)
         endif

       endif

       call zedate(bdate)
       call zetime(btime)

       write(6,*)
       write(6,'(5x,"density-functional pseudopotential plane-wave",     &
     &     " program version ",a4)') vers
       write(6,'(5x,"run on the ",a9," at ",a8)') bdate,btime
       write(6,*)
       write(6,*)

       return
       end subroutine tpage
