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

!>    gives the canonical lattice vectors for a Niggli metric

      subroutine Metric_Print(G,rede,avec,aconv,iprint)

!     written by Alvaro Ladeira, 1999
!     modified by Jose Luis Martins and
!     Manuel Maria Alemany, February 2000
!     bcc bug squashed july 2000
!     oriented the rhombohedral axes, October 2003, jlm
!     modified 29 May 2014, f90. JLM
!     modified documentation, August 2019. JLM
!     copyright Alvaro Ladeira/Jose Luis Martins/INESC-MN

!     version 4.94


      implicit none

      integer, parameter          :: REAL64 = selected_real_kind(12)

!     input

      real(REAL64), intent(in)           ::  G(3,3)                     !<  metric of the niggli lattice
      integer, intent(in)                ::  rede                       !<  type of lattice in Alvaro's list
      integer, intent(in)                ::  iprint                     !<  prints information if iprint  /=  0

!     output

      real(REAL64), intent(out)          ::  avec(3,3)                  !< "canonical" primitive lattice vectors for a Niggli metric
      real(REAL64), intent(out)          ::  aconv(3,3)                 !< conventional lattice vectors corresponding to a Niggli metric

!     local variables

      real(REAL64)       ::  a, b, c, d, e, f
      real(REAL64)       ::  alfa, cosalfa


!     constants

      real(REAL64), parameter :: ZERO = 0.0_REAL64, UM = 1.0_REAL64

!     counters

      integer         ::  i, j


      if(iprint  /=  0) then
      write(6,*) 
      write(6,*) ' REDE: ', rede
      write(6,*) 
      write(6,*) 'Final Metric:' 
      write(6,10) ( ( G(i,j),j=1,3 ), i=1,3 )
      write(6,*) 
      endif
 
      If( rede == 1 ) then
        a = sqrt(2*G(1,1))
        avec(1,1) = ZERO
        avec(2,1) = a/2
        avec(3,1) = a/2
        avec(1,2) = a/2
        avec(2,2) = ZERO
        avec(3,2) = a/2
        avec(1,3) = a/2
        avec(2,3) = a/2
        avec(3,3) = ZERO
         
        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a

        If( iprint /= 0 ) then 
        write(6,*) 'LATTICE: face centered cubic' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v1 = --- { 0, 1, 1 }' 
        write(6,*) '            2' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v2 = --- { 1, 0, 1 }' 
        write(6,*) '            2' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v3 = --- { 1, 1, 0 }' 
        write(6,*) '            2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,14) a 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 16 ) then 
        a = sqrt(G(1,1))
        avec(1,1) = a
        avec(2,1) = ZERO 
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = a
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = ZERO
        avec(3,3) = a

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple cubic' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = a { 0, 0, 1 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,14) a 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 11 ) then 
        a = 2*sqrt(G(1,1)/3)
        avec(1,1) =-a/2
        avec(2,1) = a/2
        avec(3,1) = a/2 
        avec(1,2) = a/2
        avec(2,2) =-a/2 
        avec(3,2) = a/2
        avec(1,3) = a/2
        avec(2,3) = a/2
        avec(3,3) =-a/2

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: body centered cubic' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v1 = --- { -1,  1,  1 }' 
        write(6,*) '            2' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v2 = --- {  1, -1,  1 }' 
        write(6,*) '            2' 
        write(6,*) 
 
        write(6,*) '            a' 
        write(6,*) '      v3 = --- {  1,  1, -1 }' 
        write(6,*) '            2' 

        write(6,*) 
        write(6,*) 'with' 
 
        write(6,14) a 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 7 ) then 
        a = sqrt(G(1,1))
        alfa = acos(G(1,2)/G(1,1))
        cosalfa = G(1,2)/G(1,1)
        avec(1,1) = a*sqrt(UM - cosalfa)/sqrt(2*UM)
        avec(2,1) = a*sqrt(UM - cosalfa)/sqrt(6*UM)
        avec(3,1) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)
        avec(1,2) =-a*sqrt(UM - cosalfa)/sqrt(2*UM)
        avec(2,2) = a*sqrt(UM - cosalfa)/sqrt(6*UM)
        avec(3,2) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)
        avec(1,3) = ZERO
        avec(2,3) =-a*sqrt(UM - cosalfa)*sqrt(2*UM/3*UM)
        avec(3,3) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)

        aconv(1,1) = a*sqrt(2*UM)*sqrt(UM - cosalfa)
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a*sqrt(2*UM)*sqrt(UM - cosalfa)/2 
        aconv(2,2) = a*sqrt(2*UM)*sqrt(UM - cosalfa)*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a*sqrt(UM + 2*cosalfa)*sqrt(3*UM)

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: rhombohedral-1' 
        write(6,*) 
        write(6,20) 
        write(6,*)
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by'
        write(6,*)

        write(6,'(22x,a13,13x,a13,10x,a15)') '1 - cos(alfa)',           &
     &         '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a30,a52)') 'v1 = {   a sqrt(-------------)',       &
     &  ',   a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(26x,a1,26x,a1,23x,a1)') '2','6','3'
        write(6,*)
        write(6,'(22x,a13,13x,a13,10x,a15)') '1 - cos(alfa)',            &
     &         '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a30,a52)') 'v2 = { - a sqrt(-------------)',        &
     &  ',   a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(26x,a1,26x,a1,23x,a1)') '2','6','3'
        write(6,*)
        write(6,'(46x,a15,10x,a15)') '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a6,23x,a3,a50)') 'v3 = {', '0 ,',                   &       
     &  ' -2a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(52x,a1,24x,a1)') '6','3 '

        write(6,*)
        write(6,*) 'with'

        write(6,18) a, alfa*180/3.141592654 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 24 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(G(3,3) - a*a/3) 
        avec(1,1) = a 
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = a/2
        avec(2,2) = a*sqrt(3*UM)/2
        avec(3,2) = ZERO
        avec(1,3) = a/2
        avec(2,3) = a*sqrt(3*UM)/6
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a/2 
        aconv(2,2) = a*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = 3*c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: rhombohedral-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '                1    sqrt(3) ' 
        write(6,*) '      v2 = a { ---, --------, 0 }' 
        write(6,*) '                2       2 ' 
        write(6,*) 
 
        write(6,*) '                1    sqrt(3)   c ' 
        write(6,*) '      v3 = a { ---, --------, --- }' 
        write(6,*) '                2       6      a' 
        write(6,*) 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 17 ) then 
        a = sqrt(G(1,1))
        alfa = acos(G(1,2)/G(1,1))
        cosalfa = G(1,2)/G(1,1)
        avec(1,1) = a*sqrt(UM - cosalfa)/sqrt(2*UM)
        avec(2,1) = a*sqrt(UM - cosalfa)/sqrt(6*UM)
        avec(3,1) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)
        avec(1,2) =-a*sqrt(UM - cosalfa)/sqrt(2*UM)
        avec(2,2) = a*sqrt(UM - cosalfa)/sqrt(6*UM)
        avec(3,2) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)
        avec(1,3) = ZERO
        avec(2,3) =-a*sqrt(UM - cosalfa)*sqrt(2*UM/3*UM)
        avec(3,3) = a*sqrt(UM + 2*cosalfa)/sqrt(3*UM)

        aconv(1,1) = a*sqrt(2*UM)*sqrt(UM - cosalfa)
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a*sqrt(2*UM)*sqrt(UM - cosalfa)/2 
        aconv(2,2) = a*sqrt(2*UM)*sqrt(UM - cosalfa)*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a*sqrt(UM + 2*cosalfa)*sqrt(3*UM)

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: rhombohedral-3' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by'
        write(6,*)

        write(6,'(22x,a13,13x,a13,10x,a15)') '1 - cos(alfa)',           &
     &         '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a30,a52)') 'v1 = {   a sqrt(-------------)',       &
     &  ',   a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(26x,a1,26x,a1,23x,a1)') '2','6','3'
        write(6,*)
        write(6,'(22x,a13,13x,a13,10x,a15)') '1 - cos(alfa)',            &
     &         '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a30,a52)') 'v2 = { - a sqrt(-------------)',        &
     &  ',   a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(26x,a1,26x,a1,23x,a1)') '2','6','3'
        write(6,*)
        write(6,'(46x,a15,10x,a15)') '1 - cos(alfa)','1 + 2 cos(alfa)'
        write(6,'(6x,a6,23x,a3,a50)') 'v3 = {', '0 ,',                   &       
     &  ' -2a sqrt(-------------), a sqrt(---------------)}'
        write(6,'(52x,a1,24x,a1)') '6','3 '

        write(6,*)
        write(6,*) 'with'
 
        write(6,18) a, alfa*180/3.141592654 
        write(6,*) 
        EndIf

      ElseIf( rede == 71 ) then 
        c = sqrt(G(1,1)) 
        a = sqrt(3*G(3,3) - c*c/3) 
        avec(1,1) = ZERO 
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a/2
        avec(2,2) = a*sqrt(3*UM)/6
        avec(3,2) =-c/3
        avec(1,3) =-a/2
        avec(2,3) = a*sqrt(3*UM)/6
        avec(3,3) =-c/3

        aconv(1,1) = a
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a/2 
        aconv(2,2) = a*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: rhombohedral-4' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
        write(6,*) '      v1 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '                1    sqrt(3)     c' 
        write(6,*) '      v2 = a { ---, --------, - ---- }' 
        write(6,*) '                2       6        3a'
        write(6,*) 
 
        write(6,*) '                  1    sqrt(3)     c ' 
        write(6,*) '      v3 = a { - ---, --------, - ---- }' 
        write(6,*) '                  2       6        3a' 
        write(6,*) 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 38 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(G(3,3)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) =-a/2 
        avec(2,2) = a*sqrt(3*UM)/2
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = ZERO
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a/2 
        aconv(2,2) = a*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple hexagonal-1' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '                 1    sqrt(3) ' 
        write(6,*) '      v2 = a {- ---, --------, 0 }' 
        write(6,*) '                 2       2 ' 
        write(6,*) 
 
        write(6,*) '      v3 = c { 0, 0, 1 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 64 ) then 
        c = sqrt(G(1,1)) 
        a = sqrt(G(2,2)) 
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a
        avec(2,2) = ZERO
        avec(3,2) = ZERO
        avec(1,3) =-a/2 
        avec(2,3) = a*sqrt(3*UM)/2
        avec(3,3) = ZERO

        aconv(1,1) = a
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a/2 
        aconv(2,2) = a*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple hexagonal-3' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '                  1    sqrt(3) ' 
        write(6,*) '      v3 = a { - ---, -------- , 0 }' 
        write(6,*) '                  2       2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 15 ) then 
        a = sqrt(G(1,1))
        c = sqrt(G(2,2))
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = ZERO
        avec(3,2) = c
        avec(1,3) =-a/2 
        avec(2,3) = a*sqrt(3*UM)/2
        avec(3,3) = ZERO
 
        aconv(1,1) = a
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) =-a/2 
        aconv(2,2) = a*sqrt(3*UM)/2
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple hexagonal-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = c { 0, 0, 1 }' 
 
        write(6,*) '                 1    sqrt(3) ' 
        write(6,*) '      v3 = a {- ---, --------, 0 }' 
        write(6,*) '                 2       2 ' 
        write(6,*) 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c
        write(6,*) 
        EndIf

      ElseIf( rede == 41 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(G(3,3)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = a
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = ZERO
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple tetragonal-1' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = c { 0, 0, 1 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 65 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt(G(2,2)) 
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = a
        avec(1,2) = b
        avec(2,2) = ZERO
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = b
        avec(3,3) = ZERO

        aconv(1,1) = b
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = a

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple tetragonal-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '      v2 = b { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = b { 0, 1, 0 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,15) a , b 
        write(6,*) 
        EndIf

      ElseIf( rede == 2 ) then 
        c = sqrt(G(1,1))
        a = sqrt(2*G(2,2)-c*c/2)
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a/2
        avec(2,2) = a/2
        avec(3,2) = c/2
        avec(1,3) =-a/2
        avec(2,3) = a/2
        avec(3,3) = c/2
 
        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: centered tetragonal-4' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {  0 , 0 , c }' 
        write(6,*) 
 
        write(6,*) '               a    a    c' 
        write(6,*) '      v2 = {  ---, ---, --- }' 
        write(6,*) '               2    2    2 '
        write(6,*) 
 
        write(6,*) '               a    a    c '
        write(6,*) '      v2 = {- ---, ---, --- }' 
        write(6,*) '               2    2    2 '
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c
        write(6,*) 
        EndIf

      ElseIf( rede == 12 ) then 
        a = sqrt(2*(G(1,1)+G(1,2))) 
        c = sqrt(-4*G(1,2)) 
        avec(1,1) = a/2
        avec(2,1) =-a/2
        avec(3,1) =-c/2 
        avec(1,2) = a/2
        avec(2,2) = a/2 
        avec(3,2) = c/2
        avec(1,3) =-a/2 
        avec(2,3) =-a/2
        avec(3,3) = c/2

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: centered tetragonal-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '              a     a     c ' 
        write(6,*) '      v1 = { ---,- ---,- --- }' 
        write(6,*) '              2     2     2 ' 
        write(6,*) 
 
        write(6,*) '              a    a    c ' 
        write(6,*) '      v2 = { ---, ---, --- }' 
        write(6,*) '              2    2    2 ' 
        write(6,*) 
 
        write(6,*) '               a     a    c ' 
        write(6,*) '      v3 = {- ---,- ---, --- }' 
        write(6,*) '               2     2    2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 13 ) then 
        c = sqrt(-4*G(1,3)) 
        a = sqrt(2*(G(1,1)+G(1,3))) 
        avec(1,1) =-a/2 
        avec(2,1) = a/2
        avec(3,1) = c/2
        avec(1,2) = a/2
        avec(2,2) =-a/2 
        avec(3,2) = c/2
        avec(1,3) = a/2
        avec(2,3) = a/2
        avec(3,3) =-c/2 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: centered tetragonal-3' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 

        write(6,*) '                a    a    c ' 
        write(6,*) '      v1 = { - ---, ---, --- }' 
        write(6,*) '                2    2    2 ' 
        write(6,*) 

        write(6,*) '              a     a    c ' 
        write(6,*) '      v2 = { ---,- ---, --- }'   
        write(6,*) '              2     2    2 ' 
        write(6,*) 

        write(6,*) '              a    a      c ' 
        write(6,*) '      v3 = { ---, ---, - --- }' 
        write(6,*) '              2    2      2 ' 

        write(6,*) 
        write(6,*) 'with' 

        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 35 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(4*G(3,3)-2*G(1,1)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = a
        avec(3,2) = ZERO
        avec(1,3) =-a/2 
        avec(2,3) =-a/2 
        avec(3,3) = c/2

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: centered tetragonal-1' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '               a     a    c ' 
        write(6,*) '      v3 = {- ---,- ---, --- }' 
        write(6,*) '               2     2    2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 49 ) then 
        c = sqrt(G(1,1)) 
        a = sqrt( 2*G(2,2) - G(1,1)/2) 
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a/2
        avec(2,2) = a/2
        avec(3,2) = c/2
        avec(1,3) =-a/2
        avec(2,3) = a/2
        avec(3,3) = c/2 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = a
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: centered tetragonal-5' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '              a    a    c ' 
        write(6,*) '      v2 = { ---, ---, --- }' 
        write(6,*) '              2    2    2 ' 
        write(6,*) 
 
        write(6,*) '                a    a    c ' 
        write(6,*) '      v3 = { - ---, ---, --- }' 
        write(6,*) '                2    2    2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,16) a , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 95 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt(G(2,2)) 
        c = sqrt(G(3,3)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = b
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = ZERO
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple orthorhombic' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = b { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = c { 0, 0, 1 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 20 ) then 
        a = sqrt( 2*(G(1,1) + G(1,2)) ) 
        b = sqrt( 2*(G(1,1) - G(1,2)) )
        c = sqrt( G(3,3)) 
        avec(1,1) = a/2
        avec(2,1) =-b/2
        avec(3,1) = ZERO
        avec(1,2) = a/2
        avec(2,2) = b/2 
        avec(3,2) = ZERO
        avec(1,3) = ZERO 
        avec(2,3) = ZERO
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: base centered orthorhombic-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '              a     b' 
        write(6,*) '      v1 = { ---,- ---, 0 }' 
        write(6,*) '              2     2' 
        write(6,*) 
 
        write(6,*) '              a    b' 
        write(6,*) '      v2 = { ---, ---, 0 }' 
        write(6,*) '              2    2' 
        write(6,*) 
 
        write(6,*) '      v3 = { 0, 0, c }' 
        write(6,*)  
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b, c 
        write(6,*) 
        EndIf

      ElseIf( rede == 62 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt( 4*G(2,2) - G(1,1) ) 
        c = sqrt(G(3,3))
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) =-a/2 
        avec(2,2) = b/2
        avec(3,2) = ZERO
        avec(1,3) = ZERO
        avec(2,3) = ZERO
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: base centered orthorhombic-4'  
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '               a    b' 
        write(6,*) '      v2 = {- ---, ---, 0 }' 
        write(6,*) '               2    2' 
        write(6,*) 
 
        write(6,*) '      v3 = { 0, 0, c }' 
        write(6,*) 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c
        write(6,*) 
        EndIf

      ElseIf( rede == 69 ) then 
        c = sqrt(G(1,1)) 
        a = sqrt( 2*(G(2,2) + G(2,3)) ) 
        b = sqrt( 2*(G(2,2) - G(2,3)) ) 
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a/2
        avec(2,2) =-b/2
        avec(3,2) = ZERO
        avec(1,3) = a/2
        avec(2,3) = b/2
        avec(3,3) = ZERO 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: base centered orthorhombic-3' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '              a      b' 
        write(6,*) '      v2 = { ---, - --- , 0 }' 
        write(6,*) '              2      2' 
        write(6,*) 
 
        write(6,*) '              a    b' 
        write(6,*) '      v3 = { ---, --- , 0 }' 
        write(6,*) '              2    2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 91 ) then 
        c = sqrt(G(1,1)) 
        a = sqrt(G(2,2)) 
        b = sqrt(4*G(3,3) - G(2,2)) 
        avec(1,1) = ZERO
        avec(2,1) = ZERO
        avec(3,1) = c
        avec(1,2) = a
        avec(2,2) = ZERO
        avec(3,2) = ZERO
        avec(1,3) =-a/2
        avec(2,3) = b/2 
        avec(3,3) = ZERO

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: base centered orthorhombic-8' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '                a     b ' 
        write(6,*) '      v3 = { - --- , --- , 0 }' 
        write(6,*) '                2     2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 93 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(G(2,2)) 
        b = sqrt(4*G(3,3) - G(1,1)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = ZERO
        avec(3,2) = c
        avec(1,3) =-a/2 
        avec(2,3) =-b/2
        avec(3,3) = ZERO

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: base centered orthorhombic-7' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = c { 0, 0, 1 }' 
        write(6,*) 
 
        write(6,*) '                a      b ' 
        write(6,*) '      v3 = { - ---, - ---, 0 }' 
        write(6,*) '                2      2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 25 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt(4*G(2,2) - G(1,1)) 
        c = sqrt(4*G(3,3) - G(1,1)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = a/2
        avec(2,2) = b/2
        avec(3,2) = ZERO
        avec(1,3) = a/2
        avec(2,3) = ZERO
        avec(3,3) = c/2

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: face centered orthorhombic-1' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '              a    b' 
        write(6,*) '      v2 = { ---, ---, 0 }' 
        write(6,*) '              2    2' 
        write(6,*) 
 
        write(6,*) '              a       c' 
        write(6,*) '      v3 = { ---, 0, --- }' 
        write(6,*) '              2       2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b, c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 36 ) then 
        a = sqrt(2*(G(1,1) + G(1,2))) 
        b = sqrt(2*(G(1,1) - G(1,2))) 
        c = sqrt(4*(G(3,3) + G(1,3))) 
        avec(1,1) = a/2
        avec(2,1) = b/2
        avec(3,1) = ZERO
        avec(1,2) = a/2
        avec(2,2) =-b/2 
        avec(3,2) = ZERO
        avec(1,3) =-a/2 
        avec(2,3) = ZERO
        avec(3,3) =-c/2 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: face centered orthorhombic-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '              a     b' 
        write(6,*) '      v1 = { ---,  ---, 0 }' 
        write(6,*) '              2     2' 
        write(6,*) 
 
        write(6,*) '              a     b' 
        write(6,*) '      v2 = { ---,- ---, 0 }' 
        write(6,*) '              2     2' 
        write(6,*) 
 
        write(6,*) '               a        c' 
        write(6,*) '      v3 = {- ---, 0,- --- }' 
        write(6,*) '               2        2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 3 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt(2*(G(3,3)-G(2,3)))
        c = sqrt(2*(G(3,3)+G(2,3))-G(1,1)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = a/2
        avec(2,2) = b/2
        avec(3,2) = c/2
        avec(1,3) = a/2
        avec(2,3) =-b/2 
        avec(3,3) = c/2

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: body centered orthorhombic-4' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = { a , 0 , 0 }' 
        write(6,*) 
 
        write(6,*) '              a     b    c' 
        write(6,*) '      v2 = { ---,  ---, ---}' 
        write(6,*) '              2     2    2' 
        write(6,*) 
 
        write(6,*) '              a     b    c' 
        write(6,*) '      v2 = { ---,- ---, --- }' 
        write(6,*) '              2     2    2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a, b, c
        write(6,*) 
        EndIf

      ElseIf( rede == 14 ) then 
        a = sqrt(2*( G(1,1) + G(2,3) )) 
        b = sqrt(2*( G(1,1) + G(1,3) )) 
        c = sqrt(2*( G(1,1) + G(1,2) )) 
        avec(1,1) = a/2
        avec(2,1) =-b/2 
        avec(3,1) = c/2
        avec(1,2) =-a/2 
        avec(2,2) = b/2
        avec(3,2) = c/2
        avec(1,3) =-a/2 
        avec(2,3) =-b/2 
        avec(3,3) =-c/2 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: body centered orthorhombic-5' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '              a     b     c ' 
        write(6,*) '      v1 = { ---,- --- , --- }' 
        write(6,*) '              2     2     2 ' 
        write(6,*) 
 
        write(6,*) '               a    b    c ' 
        write(6,*) '      v2 = {- ---, ---, --- }' 
        write(6,*) '               2    2    2 ' 
        write(6,*) 
 
        write(6,*) '               a     b     c ' 
        write(6,*) '      v3 = {- ---,- ---,- --- }' 
        write(6,*) '               2     2     2 ' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf

      ElseIf( rede == 87 ) then 
        a = sqrt(G(1,1)) 
        b = sqrt(G(2,2)) 
        c = sqrt(4*G(3,3) - G(1,1) - G(2,2)) 
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = b
        avec(3,2) = ZERO
        avec(1,3) =-a/2 
        avec(2,3) =-b/2
        avec(3,3) = c/2 

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = ZERO
        aconv(2,3) = ZERO
        aconv(3,3) = c

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: body centered orthorhombic-6' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = b { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '               a     b     c' 
        write(6,*) '      v3 = {- ---,- --- , --- }' 
        write(6,*) '               2     2     2' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,17) a , b , c 
        write(6,*) 
        EndIf
 
      ElseIf( rede == 42 ) then 
        a = sqrt(G(1,1)) 
        c = sqrt(G(2,2)) 
        d = G(1,3)/a
        b = sqrt(G(3,3) - d*d)
        avec(1,1) = a 
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO 
        avec(2,2) = c
        avec(3,2) = ZERO
        avec(1,3) = d
        avec(2,3) = ZERO
        avec(3,3) = b

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple monoclinic-1' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = c { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { d, 0, b}' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf
 
      ElseIf( rede == 66 ) then 
        a = sqrt(G(1,1)) 
        d = G(1,2)/a
        b = sqrt(G(2,2) - d*d)
        c = sqrt(G(3,3)) 
        avec(1,1) = a 
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = d
        avec(2,2) = ZERO
        avec(3,2) =-b
        avec(1,3) = ZERO
        avec(2,3) = c
        avec(3,3) = ZERO

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple monoclinic-2' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { d, 0, -b }' 
        write(6,*) 
 
        write(6,*) '      v3 = c { 0, 1, 0}' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf
 
      ElseIf( rede == 96 ) then 
        a = sqrt(G(2,2)) 
        d = G(2,3)/a
        b = sqrt(G(3,3) - d*d)
        c = sqrt(G(1,1)) 
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = a
        avec(2,2) = ZERO
        avec(3,2) = ZERO
        avec(1,3) = d
        avec(2,3) = ZERO
        avec(3,3) =-b

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
        write(6,*) 'LATTICE: simple monoclinic-3' 
        write(6,*) 
        write(6,20) 
        write(6,*) 
 
        write(6,2000) avec(1,1),avec(2,1),avec(3,1)
        write(6,2001) avec(1,2),avec(2,2),avec(3,2)
        write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = c { 0, 1, 0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = a { 1, 0, 0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { d, 0,-b}' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 4  ) then
        a = sqrt(G(3,3))
        d = 2*G(1,3)/a - a
        b = sqrt(2*(G(1,1)+G(1,2)) - (a+d)*(a+d))
        c = sqrt(2*(G(1,1)-G(1,2)))
        avec(1,1) = (a+d)/2
        avec(2,1) =-c/2 
        avec(3,1) =-b/2
        avec(1,2) = (a+d)/2
        avec(2,2) = c/2
        avec(3,2) =-b/2
        avec(1,3) = a
        avec(2,3) = ZERO
        avec(3,3) = ZERO

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-1'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = { (a+d)/2, -c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2,  c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = {       a,     0,   0 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 5  ) then
        a = sqrt(G(1,1))
        d = 2*G(1,2)/a - a
        b = sqrt(2*(G(2,2)+G(2,3)) - (a+d)*(a+d))
        c = sqrt(2*(G(2,2)-G(2,3)))
        avec(1,1) = a
        avec(2,1) = ZERO 
        avec(3,1) = ZERO
        avec(1,2) = (a+d)/2
        avec(2,2) = c/2
        avec(3,2) = b/2
        avec(1,3) = (a+d)/2
        avec(2,3) =-c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-2'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2,  c/2,  b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2, -c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 8  ) then
        c = sqrt(G(1,1))
        a = sqrt(G(2,2) + G(3,3) - 2*G(2,3))
        d = (G(2,2) - G(3,3))/a
        b = sqrt(4*G(2,3) + a*a -d*d -c*c)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = (a+d)/2
        avec(2,2) =-c/2
        avec(3,2) =-b/2
        avec(1,3) = (d-a)/2
        avec(2,3) =-c/2
        avec(3,3) =-b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-3'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2, -c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (d-a)/2, -c/2, -b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 9  ) then
        c = sqrt(G(1,1))
        a = sqrt(G(2,2) + G(3,3) - 2*G(2,3))
        d = (G(2,2) - G(3,3))/a
        b = sqrt(4*G(2,3) + a*a - d*d - c*c)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = (a+d)/2
        avec(2,2) = c/2
        avec(3,2) =-b/2
        avec(1,3) = (d-a)/2
        avec(2,3) = c/2
        avec(3,3) =-b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-4'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2,  c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (d-a)/2,  c/2, -b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 18 ) then
        c = sqrt(-2*(G(1,3) + G(2,3)))
        a = sqrt(G(1,1) + G(3,3) + 2*G(2,3))
        d = (G(1,1) - G(3,3))/a
        b = sqrt(4*G(1,3) + a*a + c*c - d*d)
        avec(1,1) = (a+d)/2
        avec(2,1) = c/2 
        avec(3,1) = b/2
        avec(1,2) =-(a+d)/2
        avec(2,2) = c/2
        avec(3,2) =-b/2
        avec(1,3) = (d-a)/2
        avec(2,3) =-c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-5'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = { (a+d)/2,  c/2,  b/2 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {-(a+d)/2,  c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (d-a)/2, -c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 19 ) then
        c = sqrt(G(1,1))
        a = sqrt(4*(G(2,2) - G(2,3)) + G(3,3) - c*c)
        d = (2*G(2,3) - G(3,3))/a
        b = sqrt(G(3,3) - d*d)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = (a+d)/2
        avec(2,2) = c/2
        avec(3,2) =-b/2
        avec(1,3) = d
        avec(2,3) = ZERO
        avec(3,3) =-b

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-6'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2,  c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = {       d,    0,   -b }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 21 ) then
        c = sqrt(G(1,1))
        a = sqrt(4*(G(3,3) - G(2,3)) + G(2,2) - c*c)
        d = (2*G(2,3) - G(2,2))/a
        b = sqrt(G(2,2) - d*d)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = d
        avec(2,2) = ZERO
        avec(3,2) = b
        avec(1,3) = (a+d)/2
        avec(2,3) = c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-7'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
        
        write(6,*) '      v2 = {       d,    0,    b }' 
        write(6,*) 
        
        write(6,*) '      v3 = { (a+d)/2,  c/2,  b/2 }' 
 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 26 ) then
        c = sqrt(G(1,1))
        a = sqrt(4*(G(2,2) - G(2,3)) + G(3,3) - c*c)
        d = (2*G(2,3) - G(3,3))/a
        b = sqrt(G(3,3) - d*d)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = (a+d)/2
        avec(2,2) =-c/2
        avec(3,2) =-b/2
        avec(1,3) = d
        avec(2,3) = ZERO
        avec(3,3) =-b

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-8'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2, -c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = {       d,    0,   -b }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 27 ) then
        a = sqrt(G(3,3))
        d = G(1,3)/a - a
        b = sqrt(G(1,1) - (a+d)*(a+d))
        c = sqrt(4*G(2,2) - G(1,1))
        avec(1,1) = a+d
        avec(2,1) = ZERO
        avec(3,1) =-b
        avec(1,2) = (a+d)/2
        avec(2,2) = c/2
        avec(3,2) =-b/2
        avec(1,3) = a
        avec(2,3) = ZERO
        avec(3,3) = ZERO

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-9'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {     a+d,    0,   -b }' 
        write(6,*) 
 
        write(6,*) '      v2 = { (a+d)/2,  c/2, -b/2 }' 
        write(6,*) 
 
        write(6,*) '      v3 = {       a,    0,    0 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 31 ) then
        a = sqrt(G(1,1))
        d = 2*G(1,3)/a - a
        c = sqrt(G(2,2))
        b = sqrt(4*G(3,3) - (a+d)*(a+d) - c*c)
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = c
        avec(3,2) = ZERO
        avec(1,3) = (a+d)/2
        avec(2,3) =-c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-10'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2, -c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 33 ) then
        a = sqrt(G(1,1))
        d = G(1,2)/a - a
        b = sqrt(G(2,2) - (a+d)*(a+d))
        c = sqrt(4*G(3,3) - G(2,2))
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = a+d
        avec(2,2) = ZERO  
        avec(3,2) = b
        avec(1,3) = (a+d)/2
        avec(2,3) =-c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-11'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {     a+d,    0,    b }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2, -c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 37 ) then
        a = sqrt(G(2,2))
        d = G(1,2)/a - a
        b = sqrt(G(1,1) - (a+d)*(a+d))
        c = sqrt(4*G(3,3) - G(1,1))
        avec(1,1) = a+d
        avec(2,1) = ZERO
        avec(3,1) = b
        avec(1,2) = a
        avec(2,2) = ZERO
        avec(3,2) = ZERO  
        avec(1,3) = (a+d)/2
        avec(2,3) = c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-12'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {     a+d,    0,    b }' 
        write(6,*) 
 
        write(6,*) '      v2 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2,  c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 40 ) then
        a = sqrt(G(1,1))
        d = 2*G(1,3)/a - a
        c = sqrt(G(2,2))
        b = sqrt(4*G(3,3) - (a+d)*(a+d) - c*c)
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = ZERO
        avec(2,2) = c
        avec(3,2) = ZERO
        avec(1,3) = (a+d)/2
        avec(2,3) = c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-13'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2,  c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 43 ) then
        c = sqrt(G(1,1))
        a = sqrt(4*(G(3,3) - G(2,3)) + G(2,2) - c*c)
        d = (2*G(2,3) - G(2,2))/a
        b = sqrt(G(2,2) - d*d)
        avec(1,1) = ZERO
        avec(2,1) = c
        avec(3,1) = ZERO
        avec(1,2) = d
        avec(2,2) = ZERO
        avec(3,2) = b
        avec(1,3) = (a+d)/2
        avec(2,3) =-c/2
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-14'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       0,    c,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {       d,    0,    b }' 
        write(6,*) 
 
        write(6,*) '      v3 = { (a+d)/2, -c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      ElseIf( rede == 47 ) then
        a = sqrt(G(1,1))
        d = G(1,2)/a
        b = sqrt(G(2,2) - d*d)
        c = sqrt(4*G(3,3) - (a+d)*(a+d) - b*b)
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = d
        avec(2,2) = ZERO
        avec(3,2) =-b
        avec(1,3) =-(a+d)/2 
        avec(2,3) = c/2 
        avec(3,3) = b/2

        aconv(1,1) = a 
        aconv(2,1) = ZERO
        aconv(3,1) = ZERO
        aconv(1,2) = ZERO 
        aconv(2,2) = c
        aconv(3,2) = ZERO
        aconv(1,3) = d
        aconv(2,3) = ZERO
        aconv(3,3) = b

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: centered monoclinic-15'
          write(6,*)
          write(6,20)
          write(6,*)

          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        write(6,*) 
        write(6,*) 'where the  vi  are given by' 
        write(6,*) 
 
        write(6,*) '      v1 = {       a,    0,    0 }' 
        write(6,*) 
 
        write(6,*) '      v2 = {       d,    0,   -b }' 
        write(6,*) 
 
        write(6,*) '      v3 = {-(a+d)/2,  c/2,  b/2 }' 
 
        write(6,*) 
        write(6,*) 'with' 
 
        write(6,19) a , b, c, d
        write(6,*) 
        EndIf

      Else
        a = sqrt(G(1,1))
        d = G(2,1)/a
        b = sqrt(G(2,2) - d*d)
        e = G(3,1)/a
        f = (G(3,2) - d*e)/b
        c = sqrt(G(3,3)-e*e-f*f)
        avec(1,1) = a
        avec(2,1) = ZERO
        avec(3,1) = ZERO
        avec(1,2) = d
        avec(2,2) = b
        avec(3,2) = ZERO
        avec(1,3) = e
        avec(2,3) = f
        avec(3,3) = c

        aconv(1,1) = a
        aconv(2,1) = ZERO 
        aconv(3,1) = ZERO
        aconv(1,2) = d
        aconv(2,2) = b
        aconv(3,2) = ZERO
        aconv(1,3) = e
        aconv(2,3) = f
        aconv(3,3) = c

        If( iprint /= 0 ) then
          write(6,*) 'LATTICE: triclinic-1'
          write(6,*) 
          write(6,20) 
          write(6,*) 
 
          write(6,2000) avec(1,1),avec(2,1),avec(3,1)
          write(6,2001) avec(1,2),avec(2,2),avec(3,2)
          write(6,2002) avec(1,3),avec(2,3),avec(3,3)
        EndIf
      EndIf
        

 2000 format(6X,' v1 = {',F10.6,',',F10.6,',',F10.6,' }')
 2001 format(6X,' v2 = {',F10.6,',',F10.6,',',F10.6,' }')
 2002 format(6X,' v3 = {',F10.6,',',F10.6,',',F10.6,' }')

  10  format(3F10.4)
  14  format(8X,'a = ',F10.6)
  15  format(8X,'a = ',F10.6,/,8X,'b = ',F10.6)
  16  format(8X,'a = ',F10.6,/,8X,'c = ',F10.6)
  17  format(8X,'a = ',F10.6,/,8X,'b = ',F10.6,/,8X,'c = ',F10.6)
  18  format(8X,'a  = ',F10.6,/,6X,'alfa = ',F10.6,' (graus)')
  19  format(8X,'a = ',F10.6,/,8X,'b = ',F10.6,/,8X,'c = ',F10.6,/,      &
     &       8X,'d = ',F10.6)
  20  format(' Canonic Vectors (based on primitive cell parameters):')

      return
      end subroutine Metric_Print



