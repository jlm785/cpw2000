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

!>     performs elementary transformations on a matrix (Real G)

!>     iop = 0 : initializes Mtotal = identity
!>     iop = 2 : returns in M the value of Mtotal
!>     iop = 1 : transforms G <- M G M^T
!>               updates Mtotal <- M Mtotal

      subroutine metric_Tm(iop,G,M)


!     G is real*8, M integer

!     input:
!     iop         indicates if we are initializing, calculating
!                 or recovering the final result
!     G           metric matrix
!     M           transform matrix

!     output:
!     G           transformed metric matrix (iop=1)
!     M           total transformation (iop=2)

!     SAVED: Mtotal

!     written by Alvaro Ladeira, 1999
!     modified, f90, 3 June 2014. JLM
!     copyright Alvaro Ladeira/Jose Luis Martins/INESC-MN

!     send comments/bug reports to jlmartins@inesc-mn.pt

      implicit none

      integer, parameter          :: REAL64 = selected_real_kind(12)

!     input

      integer, intent(in)            ::  iop                             !<  type of operation see above.

!     input and output

      real(REAL64), intent(inout)    ::  G(3,3)                          !<  metric matrix 
      integer, intent(inout)         ::  M(3,3)                          !<  transform matrix

!     local variables

      real(REAL64)       ::  AdotB(3,3)
      integer            ::  MdotM(3,3)

      integer, save      ::  Mtotal(3,3)

!     constants

      real(REAL64), parameter :: ZERO = 0.0_REAL64

!     counters

      integer i,j,k


      If( iop == 0 ) then

        do i=1,3
          do j=1,3
            Mtotal(j,i) = 0
          enddo
          Mtotal(i,i) = 1
        enddo

      ElseIf( iop == 1 ) then 

        do i=1,3
           do j=1,3
              MdotM(i,j) = 0
              do k=1,3
                 MdotM(i,j) = MdotM(i,j) + M(i,k)*Mtotal(k,j)
              enddo
           enddo
        enddo
        do i=1,3
           do j=1,3
              Mtotal(i,j) = MdotM(i,j)
           enddo
        enddo

        do i=1,3
           do j=1,3
              AdotB(i,j) = ZERO
              do k=1,3
                 AdotB(i,j) = AdotB(i,j) + M(i,k)*G(k,j)
              enddo
           enddo
        enddo
        do i=1,3
           do j=1,3
              G(i,j) = ZERO
              do k=1,3
                 G(i,j) = G(i,j) + AdotB(i,k)*M(j,k)
              enddo
           enddo
        enddo

      ElseIf( iop == 2 ) then 

        do i=1,3
          do j=1,3
            M(j,i) = Mtotal(j,i)
          enddo
        enddo

      EndIf

      return

      end subroutine metric_Tm

