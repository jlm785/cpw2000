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

!>     given a unit cube and a point y inside the cube chooses 
!>     4 tetrahedral corners iq and weights z such that
!>     y(i) = sum_j=0,3  z(j) iq(j,i)
!>     Uses only 5 tetrhedra.

       subroutine cube2tetra(y,z,iq)

!      written  2 April 2019.

!      copyright  J.L.Martins, CL Reis, INESC-MN.

       implicit none
       
       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)        ::  y(3)                          !<  coordinates in the cube  0<= y <= 1

!      output

       real(REAL64), intent(out)       ::  z(0:3)                        !<  coordinates in the tetrahedra  0<= z <= 1,  sum z = 1
       integer, intent(out)            ::  iq(0:3,3)                     !<  identification of tetrahedra corner iq(i,.) = 0,1

!      local

       real(REAL64)          ::  t(3), tc

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64 , UM = 1.0_REAL64
       real(REAL64), parameter :: EPS = 100*EPSILON(UM)

!      counters

       integer            ::  i, j

!      checks input

       do j = 1,3
         if(y(j) < ZERO - EPS .or. y(j) > UM + EPS) then

           write(6,'(4x,"stopped in cube2tetra, y(",i1,") = ",e25.17)')  &
     &             j, y(j)

           stop

         endif
       enddo

       do i = 1,3
         t(i) = UM - y(i)
       enddo
         
!      nearest corner t coordinate, octante

       if(y(1) + y(2) + y(3) < UM) then

         iq(0,1) = 0
         iq(0,2) = 0
         iq(0,3) = 0

         z(1) = y(1)
         z(2) = y(2)
         z(3) = y(3)

         do j = 1,3
           do i = 1,3
             iq(j,i) = iq(0,i)
           enddo
           iq(j,j) = 1 - iq(0,j)
         enddo

       elseif(y(1) + t(2) + t(3) < UM) then

         iq(0,1) = 0
         iq(0,2) = 1
         iq(0,3) = 1

         z(1) = y(1)
         z(2) = t(2)
         z(3) = t(3)

         do j = 1,3
           do i = 1,3
             iq(j,i) = iq(0,i)
           enddo
           iq(j,j) = 1 - iq(0,j)
         enddo

       elseif(t(1) + y(2) + t(3) < UM) then

         iq(0,1) = 1
         iq(0,2) = 0
         iq(0,3) = 1

         
         z(1) = t(1)
         z(2) = y(2)
         z(3) = t(3)

         do j = 1,3
           do i = 1,3
             iq(j,i) = iq(0,i)
           enddo
           iq(j,j) = 1 - iq(0,j)
         enddo

       elseif(t(1) + t(2) + y(3) < UM) then

         iq(0,1) = 1
         iq(0,2) = 1
         iq(0,3) = 0

         
         z(1) = t(1)
         z(2) = t(2)
         z(3) = y(3)

         do j = 1,3
           do i = 1,3
             iq(j,i) = iq(0,i)
           enddo
           iq(j,j) = 1 - iq(0,j)
         enddo

       else

         iq(0,1) = 1
         iq(0,2) = 1
         iq(0,3) = 1

         do j = 1,3
           do i = 1,3
             iq(j,i) = 0
           enddo
           iq(j,j) = 1
         enddo

         tc = UM - (y(1) + y(2) + y(3))
         do i = 1,3
           z(i) = tc/2 + y(i)
         enddo

       endif

       z(0) = UM - (z(1) + z(2) + z(3))

!      checks roundoff of the output

       do j = 0,3
         if(z(j) < ZERO - 10*EPS .or. z(j) > UM + 10*EPS) then

           write(6,'(4x,"stopped in cube2tetra, z(",i1,") = ",e25.17)')  &
     &             j, z(j)

           stop

         endif

         if(z(j) < ZERO) z(j) = ZERO
         if(z(j) > UM) z(j) = UM

       enddo

       end subroutine cube2tetra
