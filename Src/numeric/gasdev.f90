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

!>     Pseudo-random number generator adapted from numerical recipes
!>     Returns a normally distrubuted deviate with zero mean and
!>     unit variance

       function gasdev(iseed)

!      Written 6 January 2017, from embedded code. JLM


!      version 4.94 of pw
!      version 1.5 of md

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input and output

       integer, intent(inout)      :: iseed                              !< random number generator seed

!      output

       real(REAL64)   :: gasdev                                          !< pseudo-random number

!      local variables

       real(REAL64)         ::  v1, v2, r, fac
       real(REAL64)         ::  rnd
       integer              ::  kk
       logical              ::  lexit
       real(REAL64), save   ::  gset
       integer, save        ::  iset = 0

!      constants

       real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

!      portable random number generator (ran0 from numerical recipes 2)
!      use for less than 500000 calls

       integer, parameter       ::  IA = 16807, IM = 2147483647
       integer, parameter       ::  IQ = 127773, IR = 2836

!      counters

       integer       ::  ii

       if (iset == 0) then

         do ii = 1,1000

           kk = iseed/IQ
           iseed = IA*(iseed - kk*IQ) - IR*kk
           if(iseed < 0) iseed = iseed + IM
           rnd = (iseed*UM) / (IM*UM)
           v1 = 2*rnd - UM
           kk = iseed/IQ
           iseed = IA*(iseed-kk*IQ) - IR*kk
           if(iseed < 0) iseed = iseed + IM
           rnd = (iseed*UM) / (IM*UM)
           v2 = 2*rnd - UM
           r = v1*v1 + v2*v2
           lexit = .FALSE.
           if(r < UM .and. r /= ZERO) then
             lexit = .TRUE.
             exit
           endif
         enddo

         if(.not. lexit) then
           write(6,'("   STOPPED in gasdev:   bad random numbers?")')

           stop

         endif

         fac = sqrt(-2*log(r)/r)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
       else
         gasdev = gset
         iset = 0
       endif

       return
       end function gasdev


