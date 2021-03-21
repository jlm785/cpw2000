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

!>     random number generator adapted from numerical recipes

       function ran2(idum)

!      version 4.94 of cpw
!      version 1.5 of md

       implicit none
       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input and output
       integer, intent(inout)      :: idum                               !< random number generator seed

!      output
       real(REAL64)                :: ran2                               !< pseudo-random number

!      local varaibles

       integer        :: j
       integer, save  :: ir(97)
       integer, save  :: iy
       integer, save  :: iff = 0

!      pseudorandom number variables

       real(REAL64), parameter   :: rm = 1.4005111865831d-6
       integer, parameter        :: m = 714025
       integer, parameter        :: ia = 1366
       integer, parameter        :: ic = 150889

       if(idum < 0 .or. iff == 0)then
         iff=1
         idum=mod(ic-idum,m)
         do j=1,97
           idum=mod(ia*idum+ic,m)
           ir(j)=idum
         enddo
         idum=mod(ia*idum+ic,m)
         iy=idum
       endif
       j=1+(97*iy)/m

       if(j > 97 .or. j < 1) then
         write(6,'("  ran2:  problems with random number generator",     &
     &                 2i8)') j,idum

         stop

       endif

       iy=ir(j)
       ran2=iy*rm
       idum=mod(ia*idum+ic,m)
       ir(j)=idum

       return
       end function ran2

