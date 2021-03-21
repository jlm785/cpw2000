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

!>     generates the cartesian rotation matrices
!>     for the bravais lattices oriented according
!>     to the metric_print subroutine convention

       subroutine sym_rot(ibravais,nrot,rot)

!      Written October 28, 2003. JLM
!      Modified 7 June 2014, f90. JLM
!      Modified 11 June 2014, monoclinic orientation. JLM
!      Modified, documentation, December 2019. JLM
!      copyright INESC-MN/Jose Luis Martins


!      version 4.94


       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  ibravais                   !<  index of bravais lattice

!      output

       real(REAL64), intent(out)          ::  rot(3,3,48)                !<  n-th cartesian rotation matrix
       integer, intent(out)               ::  nrot                       !<  number of symetry operations

!      other variables

       integer    ::  nhalf

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  :: UMEIO = 0.5_REAL64
       real(REAL64), parameter  :: HEX = 0.866025403784438647_REAL64

!      counters

       integer    ::  i, j, n
       

       nhalf = 0

       if(ibravais == 1 .or. ibravais == 2 .or. ibravais == 3) then

!        cubic lattices

         nhalf = 24

         rot(1,1, 1) =     UM
         rot(1,2, 1) =     ZERO
         rot(1,3, 1) =     ZERO
         rot(2,1, 1) =     ZERO
         rot(2,2, 1) =     UM
         rot(2,3, 1) =     ZERO
         rot(3,1, 1) =     ZERO
         rot(3,2, 1) =     ZERO
         rot(3,3, 1) =     UM

         rot(1,1, 2) =     UM
         rot(1,2, 2) =     ZERO
         rot(1,3, 2) =     ZERO
         rot(2,1, 2) =     ZERO
         rot(2,2, 2) =    -UM
         rot(2,3, 2) =     ZERO
         rot(3,1, 2) =     ZERO
         rot(3,2, 2) =     ZERO
         rot(3,3, 2) =    -UM

         rot(1,1, 3) =    -UM
         rot(1,2, 3) =     ZERO
         rot(1,3, 3) =     ZERO
         rot(2,1, 3) =     ZERO
         rot(2,2, 3) =     UM
         rot(2,3, 3) =     ZERO
         rot(3,1, 3) =     ZERO
         rot(3,2, 3) =     ZERO
         rot(3,3, 3) =    -UM

         rot(1,1, 4) =    -UM
         rot(1,2, 4) =     ZERO
         rot(1,3, 4) =     ZERO
         rot(2,1, 4) =     ZERO
         rot(2,2, 4) =    -UM
         rot(2,3, 4) =     ZERO
         rot(3,1, 4) =     ZERO
         rot(3,2, 4) =     ZERO
         rot(3,3, 4) =     UM

         rot(1,1, 5) =     ZERO
         rot(1,2, 5) =     UM
         rot(1,3, 5) =     ZERO
         rot(2,1, 5) =     ZERO
         rot(2,2, 5) =     ZERO
         rot(2,3, 5) =     UM
         rot(3,1, 5) =     UM
         rot(3,2, 5) =     ZERO
         rot(3,3, 5) =     ZERO

         rot(1,1, 6) =     ZERO
         rot(1,2, 6) =     UM
         rot(1,3, 6) =     ZERO
         rot(2,1, 6) =     ZERO
         rot(2,2, 6) =     ZERO
         rot(2,3, 6) =    -UM
         rot(3,1, 6) =    -UM
         rot(3,2, 6) =     ZERO
         rot(3,3, 6) =     ZERO

         rot(1,1, 7) =     ZERO
         rot(1,2, 7) =    -UM
         rot(1,3, 7) =     ZERO
         rot(2,1, 7) =     ZERO
         rot(2,2, 7) =     ZERO
         rot(2,3, 7) =     UM
         rot(3,1, 7) =    -UM
         rot(3,2, 7) =     ZERO
         rot(3,3, 7) =     ZERO

         rot(1,1, 8) =     ZERO
         rot(1,2, 8) =    -UM
         rot(1,3, 8) =     ZERO
         rot(2,1, 8) =     ZERO
         rot(2,2, 8) =     ZERO
         rot(2,3, 8) =    -UM
         rot(3,1, 8) =     UM
         rot(3,2, 8) =     ZERO
         rot(3,3, 8) =     ZERO

         rot(1,1, 9) =     ZERO
         rot(1,2, 9) =     ZERO
         rot(1,3, 9) =     UM
         rot(2,1, 9) =     UM
         rot(2,2, 9) =     ZERO
         rot(2,3, 9) =     ZERO
         rot(3,1, 9) =     ZERO
         rot(3,2, 9) =     UM
         rot(3,3, 9) =     ZERO

         rot(1,1,10) =     ZERO
         rot(1,2,10) =     ZERO
         rot(1,3,10) =     UM
         rot(2,1,10) =    -UM
         rot(2,2,10) =     ZERO
         rot(2,3,10) =     ZERO
         rot(3,1,10) =     ZERO
         rot(3,2,10) =    -UM
         rot(3,3,10) =     ZERO

         rot(1,1,11) =     ZERO
         rot(1,2,11) =     ZERO
         rot(1,3,11) =    -UM
         rot(2,1,11) =     UM
         rot(2,2,11) =     ZERO
         rot(2,3,11) =     ZERO
         rot(3,1,11) =     ZERO
         rot(3,2,11) =    -UM
         rot(3,3,11) =     ZERO

         rot(1,1,12) =     ZERO
         rot(1,2,12) =     ZERO
         rot(1,3,12) =    -UM
         rot(2,1,12) =    -UM
         rot(2,2,12) =     ZERO
         rot(2,3,12) =     ZERO
         rot(3,1,12) =     ZERO
         rot(3,2,12) =     UM
         rot(3,3,12) =     ZERO

         rot(1,1,13) =     ZERO
         rot(1,2,13) =    -UM
         rot(1,3,13) =     ZERO
         rot(2,1,13) =    -UM
         rot(2,2,13) =     ZERO
         rot(2,3,13) =     ZERO
         rot(3,1,13) =     ZERO
         rot(3,2,13) =     ZERO
         rot(3,3,13) =    -UM

         rot(1,1,14) =     ZERO
         rot(1,2,14) =    -UM
         rot(1,3,14) =     ZERO
         rot(2,1,14) =     UM
         rot(2,2,14) =     ZERO
         rot(2,3,14) =     ZERO
         rot(3,1,14) =     ZERO
         rot(3,2,14) =     ZERO
         rot(3,3,14) =     UM

         rot(1,1,15) =     ZERO
         rot(1,2,15) =     UM
         rot(1,3,15) =     ZERO
         rot(2,1,15) =    -UM
         rot(2,2,15) =     ZERO
         rot(2,3,15) =     ZERO
         rot(3,1,15) =     ZERO
         rot(3,2,15) =     ZERO
         rot(3,3,15) =     UM

         rot(1,1,16) =     ZERO
         rot(1,2,16) =     UM
         rot(1,3,16) =     ZERO
         rot(2,1,16) =     UM
         rot(2,2,16) =     ZERO
         rot(2,3,16) =     ZERO
         rot(3,1,16) =     ZERO
         rot(3,2,16) =     ZERO
         rot(3,3,16) =    -UM

         rot(1,1,17) =    -UM
         rot(1,2,17) =     ZERO
         rot(1,3,17) =     ZERO
         rot(2,1,17) =     ZERO
         rot(2,2,17) =     ZERO
         rot(2,3,17) =    -UM
         rot(3,1,17) =     ZERO
         rot(3,2,17) =    -UM
         rot(3,3,17) =     ZERO

         rot(1,1,18) =    -UM
         rot(1,2,18) =     ZERO
         rot(1,3,18) =     ZERO
         rot(2,1,18) =     ZERO
         rot(2,2,18) =     ZERO
         rot(2,3,18) =     UM
         rot(3,1,18) =     ZERO
         rot(3,2,18) =     UM
         rot(3,3,18) =     ZERO

         rot(1,1,19) =     UM
         rot(1,2,19) =     ZERO
         rot(1,3,19) =     ZERO
         rot(2,1,19) =     ZERO
         rot(2,2,19) =     ZERO
         rot(2,3,19) =    -UM
         rot(3,1,19) =     ZERO
         rot(3,2,19) =     UM
         rot(3,3,19) =     ZERO

         rot(1,1,20) =     UM
         rot(1,2,20) =     ZERO
         rot(1,3,20) =     ZERO
         rot(2,1,20) =     ZERO
         rot(2,2,20) =     ZERO
         rot(2,3,20) =     UM
         rot(3,1,20) =     ZERO
         rot(3,2,20) =    -UM
         rot(3,3,20) =     ZERO

         rot(1,1,21) =     ZERO
         rot(1,2,21) =     ZERO
         rot(1,3,21) =    -UM
         rot(2,1,21) =     ZERO
         rot(2,2,21) =    -UM
         rot(2,3,21) =     ZERO
         rot(3,1,21) =    -UM
         rot(3,2,21) =     ZERO
         rot(3,3,21) =     ZERO

         rot(1,1,22) =     ZERO
         rot(1,2,22) =     ZERO
         rot(1,3,22) =    -UM
         rot(2,1,22) =     ZERO
         rot(2,2,22) =     UM
         rot(2,3,22) =     ZERO
         rot(3,1,22) =     UM
         rot(3,2,22) =     ZERO
         rot(3,3,22) =     ZERO

         rot(1,1,23) =     ZERO
         rot(1,2,23) =     ZERO
         rot(1,3,23) =     UM
         rot(2,1,23) =     ZERO
         rot(2,2,23) =    -UM
         rot(2,3,23) =     ZERO
         rot(3,1,23) =     UM
         rot(3,2,23) =     ZERO
         rot(3,3,23) =     ZERO

         rot(1,1,24) =     ZERO
         rot(1,2,24) =     ZERO
         rot(1,3,24) =     UM
         rot(2,1,24) =     ZERO
         rot(2,2,24) =     UM
         rot(2,3,24) =     ZERO
         rot(3,1,24) =    -UM
         rot(3,2,24) =     ZERO
         rot(3,3,24) =     ZERO

       elseif(ibravais == 4 .or. ibravais == 5) then

!        tetragonal lattices

         nhalf = 8

         rot(1,1, 1) =     UM
         rot(1,2, 1) =     ZERO
         rot(1,3, 1) =     ZERO
         rot(2,1, 1) =     ZERO
         rot(2,2, 1) =     UM
         rot(2,3, 1) =     ZERO
         rot(3,1, 1) =     ZERO
         rot(3,2, 1) =     ZERO
         rot(3,3, 1) =     UM

         rot(1,1, 2) =     UM
         rot(1,2, 2) =     ZERO
         rot(1,3, 2) =     ZERO
         rot(2,1, 2) =     ZERO
         rot(2,2, 2) =    -UM
         rot(2,3, 2) =     ZERO
         rot(3,1, 2) =     ZERO
         rot(3,2, 2) =     ZERO
         rot(3,3, 2) =    -UM

         rot(1,1, 3) =    -UM
         rot(1,2, 3) =     ZERO
         rot(1,3, 3) =     ZERO
         rot(2,1, 3) =     ZERO
         rot(2,2, 3) =     UM
         rot(2,3, 3) =     ZERO
         rot(3,1, 3) =     ZERO
         rot(3,2, 3) =     ZERO
         rot(3,3, 3) =    -UM

         rot(1,1, 4) =    -UM
         rot(1,2, 4) =     ZERO
         rot(1,3, 4) =     ZERO
         rot(2,1, 4) =     ZERO
         rot(2,2, 4) =    -UM
         rot(2,3, 4) =     ZERO
         rot(3,1, 4) =     ZERO
         rot(3,2, 4) =     ZERO
         rot(3,3, 4) =     UM

         rot(1,1, 5) =     ZERO
         rot(1,2, 5) =    -UM
         rot(1,3, 5) =     ZERO
         rot(2,1, 5) =    -UM
         rot(2,2, 5) =     ZERO
         rot(2,3, 5) =     ZERO
         rot(3,1, 5) =     ZERO
         rot(3,2, 5) =     ZERO
         rot(3,3, 5) =    -UM

         rot(1,1, 6) =     ZERO
         rot(1,2, 6) =    -UM
         rot(1,3, 6) =     ZERO
         rot(2,1, 6) =     UM
         rot(2,2, 6) =     ZERO
         rot(2,3, 6) =     ZERO
         rot(3,1, 6) =     ZERO
         rot(3,2, 6) =     ZERO
         rot(3,3, 6) =     UM

         rot(1,1, 7) =     ZERO
         rot(1,2, 7) =     UM
         rot(1,3, 7) =     ZERO
         rot(2,1, 7) =    -UM
         rot(2,2, 7) =     ZERO
         rot(2,3, 7) =     ZERO
         rot(3,1, 7) =     ZERO
         rot(3,2, 7) =     ZERO
         rot(3,3, 7) =     UM

         rot(1,1, 8) =     ZERO
         rot(1,2, 8) =     UM
         rot(1,3, 8) =     ZERO
         rot(2,1, 8) =     UM
         rot(2,2, 8) =     ZERO
         rot(2,3, 8) =     ZERO
         rot(3,1, 8) =     ZERO
         rot(3,2, 8) =     ZERO
         rot(3,3, 8) =    -UM

       elseif(ibravais >= 6 .and. ibravais <= 9) then

!        ortorhombic lattices

         nhalf = 4

         rot(1,1, 1) =     UM
         rot(1,2, 1) =     ZERO
         rot(1,3, 1) =     ZERO
         rot(2,1, 1) =     ZERO
         rot(2,2, 1) =     UM
         rot(2,3, 1) =     ZERO
         rot(3,1, 1) =     ZERO
         rot(3,2, 1) =     ZERO
         rot(3,3, 1) =     UM

         rot(1,1, 2) =     UM
         rot(1,2, 2) =     ZERO
         rot(1,3, 2) =     ZERO
         rot(2,1, 2) =     ZERO
         rot(2,2, 2) =    -UM
         rot(2,3, 2) =     ZERO
         rot(3,1, 2) =     ZERO
         rot(3,2, 2) =     ZERO
         rot(3,3, 2) =    -UM

         rot(1,1, 3) =    -UM
         rot(1,2, 3) =     ZERO
         rot(1,3, 3) =     ZERO
         rot(2,1, 3) =     ZERO
         rot(2,2, 3) =     UM
         rot(2,3, 3) =     ZERO
         rot(3,1, 3) =     ZERO
         rot(3,2, 3) =     ZERO
         rot(3,3, 3) =    -UM

         rot(1,1, 4) =    -UM
         rot(1,2, 4) =     ZERO
         rot(1,3, 4) =     ZERO
         rot(2,1, 4) =     ZERO
         rot(2,2, 4) =    -UM
         rot(2,3, 4) =     ZERO
         rot(3,1, 4) =     ZERO
         rot(3,2, 4) =     ZERO
         rot(3,3, 4) =     UM

       elseif(ibravais == 10) then

!        hexagonal lattice

         nhalf = 12

         rot(1,1, 1) =  UM
         rot(1,2, 1) =  ZERO
         rot(1,3, 1) =  ZERO
         rot(2,1, 1) =  ZERO
         rot(2,2, 1) =  UM
         rot(2,3, 1) =  ZERO
         rot(3,1, 1) =  ZERO
         rot(3,2, 1) =  ZERO
         rot(3,3, 1) =  UM

         rot(1,1, 2) =  UMEIO
         rot(1,2, 2) = -HEX
         rot(1,3, 2) =  ZERO
         rot(2,1, 2) =  HEX
         rot(2,2, 2) =  UMEIO
         rot(2,3, 2) =  ZERO
         rot(3,1, 2) =  ZERO
         rot(3,2, 2) =  ZERO
         rot(3,3, 2) =  UM

         rot(1,1, 3) = -UMEIO
         rot(1,2, 3) = -HEX
         rot(1,3, 3) =  ZERO
         rot(2,1, 3) =  HEX
         rot(2,2, 3) = -UMEIO
         rot(2,3, 3) =  ZERO
         rot(3,1, 3) =  ZERO
         rot(3,2, 3) =  ZERO
         rot(3,3, 3) =  UM

         rot(1,1, 4) = -UM
         rot(1,2, 4) =  ZERO
         rot(1,3, 4) =  ZERO
         rot(2,1, 4) =  ZERO
         rot(2,2, 4) = -UM
         rot(2,3, 4) =  ZERO
         rot(3,1, 4) =  ZERO
         rot(3,2, 4) =  ZERO
         rot(3,3, 4) =  UM

         rot(1,1, 5) = -UMEIO
         rot(1,2, 5) =  HEX
         rot(1,3, 5) =  ZERO
         rot(2,1, 5) = -HEX
         rot(2,2, 5) = -UMEIO
         rot(2,3, 5) =  ZERO
         rot(3,1, 5) =  ZERO
         rot(3,2, 5) =  ZERO
         rot(3,3, 5) =  UM

         rot(1,1, 6) =  UMEIO
         rot(1,2, 6) =  HEX
         rot(1,3, 6) =  ZERO
         rot(2,1, 6) = -HEX
         rot(2,2, 6) =  UMEIO
         rot(2,3, 6) =  ZERO
         rot(3,1, 6) =  ZERO
         rot(3,2, 6) =  ZERO
         rot(3,3, 6) =  UM

         rot(1,1, 7) = -UMEIO
         rot(1,2, 7) = -HEX
         rot(1,3, 7) =  ZERO
         rot(2,1, 7) = -HEX
         rot(2,2, 7) =  UMEIO
         rot(2,3, 7) =  ZERO
         rot(3,1, 7) =  ZERO
         rot(3,2, 7) =  ZERO
         rot(3,3, 7) = -UM

         rot(1,1, 8) =  UMEIO
         rot(1,2, 8) = -HEX
         rot(1,3, 8) =  ZERO
         rot(2,1, 8) = -HEX
         rot(2,2, 8) = -UMEIO
         rot(2,3, 8) =  ZERO
         rot(3,1, 8) =  ZERO
         rot(3,2, 8) =  ZERO
         rot(3,3, 8) = -UM

         rot(1,1, 9) =  UM
         rot(1,2, 9) =  ZERO
         rot(1,3, 9) =  ZERO
         rot(2,1, 9) =  ZERO
         rot(2,2, 9) = -UM
         rot(2,3, 9) =  ZERO
         rot(3,1, 9) =  ZERO
         rot(3,2, 9) =  ZERO
         rot(3,3, 9) = -UM

         rot(1,1,10) =  UMEIO
         rot(1,2,10) =  HEX
         rot(1,3,10) =  ZERO
         rot(2,1,10) =  HEX
         rot(2,2,10) = -UMEIO
         rot(2,3,10) =  ZERO
         rot(3,1,10) =  ZERO
         rot(3,2,10) =  ZERO
         rot(3,3,10) = -UM

         rot(1,1,11) = -UMEIO
         rot(1,2,11) =  HEX
         rot(1,3,11) =  ZERO
         rot(2,1,11) =  HEX
         rot(2,2,11) =  UMEIO
         rot(2,3,11) =  ZERO
         rot(3,1,11) =  ZERO
         rot(3,2,11) =  ZERO
         rot(3,3,11) = -UM

         rot(1,1,12) = -UM
         rot(1,2,12) =  ZERO
         rot(1,3,12) =  ZERO
         rot(2,1,12) =  ZERO
         rot(2,2,12) =  UM
         rot(2,3,12) =  ZERO
         rot(3,1,12) =  ZERO
         rot(3,2,12) =  ZERO
         rot(3,3,12) = -UM

       elseif(ibravais == 11) then

!        rhombohedral lattice

         nhalf = 6

         rot(1,1, 1) =  UM
         rot(1,2, 1) =  ZERO
         rot(1,3, 1) =  ZERO
         rot(2,1, 1) =  ZERO
         rot(2,2, 1) =  UM
         rot(2,3, 1) =  ZERO
         rot(3,1, 1) =  ZERO
         rot(3,2, 1) =  ZERO
         rot(3,3, 1) =  UM

         rot(1,1, 2) = -UMEIO
         rot(1,2, 2) = -HEX
         rot(1,3, 2) =  ZERO
         rot(2,1, 2) =  HEX
         rot(2,2, 2) = -UMEIO
         rot(2,3, 2) =  ZERO
         rot(3,1, 2) =  ZERO
         rot(3,2, 2) =  ZERO
         rot(3,3, 2) =  UM

         rot(1,1, 3) = -UMEIO
         rot(1,2, 3) =  HEX
         rot(1,3, 3) =  ZERO
         rot(2,1, 3) = -HEX
         rot(2,2, 3) = -UMEIO
         rot(2,3, 3) =  ZERO
         rot(3,1, 3) =  ZERO
         rot(3,2, 3) =  ZERO
         rot(3,3, 3) =  UM

         rot(1,1, 4) = -UMEIO
         rot(1,2, 4) = -HEX
         rot(1,3, 4) =  ZERO
         rot(2,1, 4) = -HEX
         rot(2,2, 4) =  UMEIO
         rot(2,3, 4) =  ZERO
         rot(3,1, 4) =  ZERO
         rot(3,2, 4) =  ZERO
         rot(3,3, 4) = -UM

         rot(1,1, 5) =  UM
         rot(1,2, 5) =  ZERO
         rot(1,3, 5) =  ZERO
         rot(2,1, 5) =  ZERO
         rot(2,2, 5) = -UM
         rot(2,3, 5) =  ZERO
         rot(3,1, 5) =  ZERO
         rot(3,2, 5) =  ZERO
         rot(3,3, 5) = -UM

         rot(1,1, 6) = -UMEIO
         rot(1,2, 6) =  HEX
         rot(1,3, 6) =  ZERO
         rot(2,1, 6) =  HEX
         rot(2,2, 6) =  UMEIO
         rot(2,3, 6) =  ZERO
         rot(3,1, 6) =  ZERO
         rot(3,2, 6) =  ZERO
         rot(3,3, 6) = -UM

       elseif(ibravais == 12 .or. ibravais == 13) then

!        monoclinic lattice

         nhalf = 2

         rot(1,1, 1) =     UM
         rot(1,2, 1) =     ZERO
         rot(1,3, 1) =     ZERO
         rot(2,1, 1) =     ZERO
         rot(2,2, 1) =     UM
         rot(2,3, 1) =     ZERO
         rot(3,1, 1) =     ZERO
         rot(3,2, 1) =     ZERO
         rot(3,3, 1) =     UM

         rot(1,1, 2) =    -UM
         rot(1,2, 2) =     ZERO
         rot(1,3, 2) =     ZERO
         rot(2,1, 2) =     ZERO
         rot(2,2, 2) =     UM
         rot(2,3, 2) =     ZERO
         rot(3,1, 2) =     ZERO
         rot(3,2, 2) =     ZERO
         rot(3,3, 2) =    -UM

       else

!        triclinic lattice

         nhalf = 1

         rot(1,1, 1) =     UM
         rot(1,2, 1) =     ZERO
         rot(1,3, 1) =     ZERO
         rot(2,1, 1) =     ZERO
         rot(2,2, 1) =     UM
         rot(2,3, 1) =     ZERO
         rot(3,1, 1) =     ZERO
         rot(3,2, 1) =     ZERO
         rot(3,3, 1) =     UM

       endif

!      now add the inversion

       nrot = 2*nhalf

       do n=1,nhalf
       do j=1,3
       do i=1,3
         rot(i,j,n+nhalf) =-rot(i,j,n)
       enddo
       enddo
       enddo

       return

       end subroutine sym_rot
