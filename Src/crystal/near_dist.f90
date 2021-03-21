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

!>     closest distance between atoms taking into account periodicity

       subroutine near_dist(distmin,adot,r0,r1)

!      Written 30 November 2016.  JLM
!      Modified, documentation, August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       real(REAL64), intent(in)           ::  r0(3), r1(3)               !<  positions

!      output

       real(REAL64), intent(out)          ::  distmin                    !<  closest distance between images

!      local variables

       real(REAL64)    ::  xk(3), zk(3), dist

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64

!      counters

       integer i, j, k1, k2, k3

       do i = 1,3
         xk(i) = r1(i) - r0(i)
         xk(i) = xk(i) - nint(xk(i))
       enddo
       
       distmin = ZERO
       do i = 1,3
       do j = 1,3
         distmin = distmin + xk(i)*adot(i,j)*xk(j)
       enddo
       enddo

       do k1 = -1,1
       do k2 = -1,1
       do k3 = -1,1

         zk(1) = xk(1) + k1
         zk(2) = xk(2) + k2
         zk(3) = xk(3) + k3

         dist = ZERO
         do i = 1,3
         do j = 1,3
           dist = dist + zk(i)*adot(i,j)*zk(j)
         enddo
         enddo

         if(dist < distmin) distmin = dist

       enddo
       enddo
       enddo
       
       return
       end subroutine near_dist

