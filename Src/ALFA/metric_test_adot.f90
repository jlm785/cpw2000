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

!>    given a 3x3 matrix g tests if it corresponds to a metric

      subroutine metric_test_g(g,lmetric)

!     Written 8 December 2016 from metric_ident. JLM
!     Modified, documentation, June 2019. JLM


!     version 4.94


      implicit none

      integer, parameter          :: REAL64 = selected_real_kind(12)

!     input and output

      real(REAL64), intent(inout)        ::  g(3,3)                      !<  input metric, symmetrize don output

!     output

      logical, intent(out)               ::  lmetric                     !<  True if g is a metric

!     local variables

      real(REAL64)         ::  scal, detg

!     constants

      real(REAL64), parameter :: ZERO = 0.0_REAL64
      real(REAL64), parameter :: EPS = 1.0d-12

!     counters

      integer         ::  i, j


      lmetric = .TRUE.

!        TEST: Is G symmetric?

      scal = ZERO
      do i=1,3
      do j=1,3
        scal = max(scal,abs(g(i,j)))
      enddo
      enddo

      if(scal == ZERO) then

        write(6,*) ' metric_ident'
        write(6,*) '  Metric is zero'

        lmetric = .FALSE.

      endif


      If( (abs(g(1,2)-g(2,1)) > EPS*scal) .or.                           &
     &    (abs(g(1,3)-g(3,1)) > EPS*scal) .or.                           &
     &    (abs(g(2,3)-g(3,2)) > EPS*scal) ) then

         write(6,*) ' metric_ident'
         write(6,*) ' g_ij-g_ji ',g(1,2)-g(2,1),g(1,3)-g(3,1),           &
     &                          g(2,3)-g(3,2)
         write(6,*) ' metric_ident'
         write(6,*) 'Metric non simmetric!   exceeds tolerance', EPS

         lmetric = .FALSE.

      EndIf

      g(1,2) = (g(1,2)+g(2,1)) / 2
      g(2,1) = g(1,2)
      g(1,3) = (g(1,3)+g(3,1)) / 2
      g(3,1) = g(1,3)
      g(2,3) = (g(2,3)+g(3,2)) / 2
      g(3,2) = g(2,3)


!        TEST:  Is there any non positive diagonal element in  G ?

      If( (g(1,1) <= 0.).or.(g(2,2) <= 0.).or.(g(3,3) <= 0.) ) then

         write(6,*) ' metric_ident'
         write(6,*) ' g_ii ',g(1,1),g(2,2),g(3,3)
         write(6,*) 'Metric diagonal elements must be positive!'

         lmetric = .FALSE.

      EndIf

!        TEST:  Det(G) = 0 ?

      detg = g(1,1)*g(2,2)*g(3,3) - g(1,1)*g(2,3)*g(3,2)                 &
     &     + g(1,2)*g(2,3)*g(3,1) - g(1,2)*g(2,1)*g(3,3)                 &
     &     + g(1,3)*g(2,1)*g(3,2) - g(1,3)*g(2,2)*g(3,1)

      If( detg <= ZERO) then

         write(6,*) ' metric_ident'
         write(6,*) ' det(g) ',detg
         write(6,*) ' Det(g) less or equal to 0 !'

         lmetric = .FALSE.

      EndIf

      return

      end subroutine metric_test_g
 
 
