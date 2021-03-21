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

!>     multiplies integer matrices 3x3

       subroutine metric_dotM(Prod,M1,M2)

!      Prod = M1*M2

!      input:
!      M1         first matrix
!      M2         second matrix

!      output:
!      Prod       product matrix

!      written by Alvaro Ladeira, 1999
!      modified, f90, 3 June 2014. JLM
!      Modified documentation, August 2019.
!      copyright Alvaro Ladeira/Jose Luis Martins/INESC-MN

!      send comments/bug reports to jlmartins@inesc-mn.pt

!      version 4.94 

       implicit none

!      input

       integer, intent(in)     ::  M1(3,3)                               !<  first matrix
       integer, intent(in)     ::  M2(3,3)                               !<  second matrix

!      output

       integer, intent(out)    ::  Prod(3,3)                             !<  Prod = M1 M2

!      counters

       integer i,j,k

       do i=1,3
          do j=1,3
             Prod(i,j)=0
             do k=1,3
                Prod(i,j) = Prod(i,j) + M1(i,k)*M2(k,j)
             enddo
          enddo
       enddo

       return
       end subroutine metric_dotM
