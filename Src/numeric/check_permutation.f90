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

!>     checks if indx is a permutation

       subroutine check_permutation(indx,ntype,lperm)

!      written June 2017. JLM
!      Copyright INESC-MN/Jose Luis Martins

!      version 4.94 of pw
!      version 1.5 of md

       implicit none

!      input

       integer, intent(in)                ::  ntype                      !> number of types of atoms
       integer, intent(in)                ::  indx(ntype)                !> data to be checked

!      output

       logical, intent(out)               ::  lperm                      !> indicates if data is a permutation.

!      local variables

       logical, allocatable               ::  lpres(:)

!      counters

       integer     ::  nt

       lperm = .TRUE.
       
!      checks range

       do nt = 1,ntype
         if(indx(nt) < 1 .OR. indx(nt) > ntype) then
           lperm = .FALSE.

           exit

         endif
       enddo

       if(lperm) then

         allocate(lpres(ntype))

         do nt = 1,ntype
           lpres(nt) = .FALSE.
         enddo

         do nt = 1,ntype
           lpres(indx(nt)) = .TRUE.
         enddo

         do nt = 1,ntype

           if(.NOT. lpres(nt)) then
             lperm = .FALSE.

             exit

           endif
         enddo

         deallocate(lpres)

       endif

       return
       end subroutine check_permutation

