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

!>     This program finds shift that brings atom to shift_origin

       subroutine shift_origin(adot,ntype,natom,rat,shift,ratio,         &
     &     mxdtyp,mxdatm)

!      writen January 10, 2017. JLM
!      Modified, documentation, August 2019. JLM
!      copyright  Jose Luis Martins/INESC-MN

!      Version 4.94 of cpw


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)


!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of types of atoms

       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  k-th component (in lattice coordinates) of the position of the n-th atom of type i

!      output

       real(REAL64), intent(out)          ::  shift(3)                   !<  shift that brings atom to shift_origin
       real(REAL64), intent(out)          ::  ratio                      !<  ratio of nearest to second nearest atom

!      local variables

       real(REAL64)       ::  dist, distmin, distmin2
       integer            ::  ntmin, jmin
       real(REAL64)       ::  r0(3),r1(3)
 
!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!      counters

       integer i, j, nt


       do i = 1,3
         r0(i) = ZERO
         r1(i) = rat(i,1,1)
       enddo
       
       call near_dist(distmin,adot,r0,r1)
       distmin2 = distmin

       
       if(ntype == 1 .and. natom(1) == 1) then

         do i = 1,3
           shift(i) = -rat(i,1,1)
         enddo
         ratio = UM
       
       else

         ntmin = 1
         jmin = 1
         
         do nt = 1,ntype
         do j = 1,natom(nt)

           do i = 1,3
             r1(i) = rat(i,j,nt)
           enddo

           call near_dist(dist,adot,r0,r1)

           if(dist < distmin) then
             ntmin = nt
             jmin = j
             distmin = dist
           endif
            
         enddo
         enddo

         distmin2 = 10*distmin + 10000000.0

         do nt = 1,ntype
         do j = 1,natom(nt)

           if(nt /= ntmin .or. j /= jmin) then
             do i = 1,3
               r1(i) = rat(i,j,nt)
             enddo

             call near_dist(dist,adot,r0,r1)

             if(dist < distmin2) then
               distmin2 = dist
             endif
 
           endif

         enddo
         enddo

         do i = 1,3
           shift(i) = -rat(i,jmin,ntmin)
         enddo
         ratio = distmin/distmin2

       endif

       return
       end subroutine shift_origin
