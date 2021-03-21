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

!>     prints the crystal structure
!>     and checks the properties of the metric

       subroutine print_crystal(ipr,adot,ntype,natom,nameat,rat,         &
     & mxdtyp,mxdatm)

!      written 15 october 1993. jlm
!      modified by ivo souza
!      modified 15 january 1999. jlm
!      modified for f90, 22 May 2014. jlm
!      included CLR modifications, 11 September 2015. JLM
!      Moved CLR modifications to another subroutine, 20 October 2015. JLM
!      Modified, documentation June 2019.
!      copyright INESC-MN/Jose Luis Martins/Carlos Loia Reis/Ivo Souza

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       integer, intent(in)                ::  ipr                        !<  should be equal to one if information is to be printed.

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i
       real(REAL64), intent(in)           ::  rat(3,mxdatm,mxdtyp)       !<  lattice coordinates of atom j of type i

!      local variables

       real(REAL64)  ::    angl12,angl23,angl13
       real(REAL64)  ::    avec(3,3),bvec(3,3)
       real(REAL64)  ::    vcell,bdot(3,3)
       real(REAL64)  ::    car(3),trace

!      parameters

       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.00000001_REAL64

!      counters

       integer       ::  i, ja, jmax, nt, ntt


!      checks symmetry and positiveness of adot

       trace = adot(1,1) + adot(2,2) + adot(3,3)

       if(abs(adot(1,2)-adot(2,1))/trace > EPS) then
         write(6,'("  Stopped in print crystal.  adot is not ",          &
     &      "symmetric: ")')
         write(6,'("  adot(",i1,",",i1,") - adot(",i1,",",i1,") = ",     &
     &      e12.4)') 1,2,2,1,adot(1,2)-adot(2,1)

         stop

       endif

       if(abs(adot(1,3)-adot(3,1))/trace > EPS) then
         write(6,'("  Stopped in print crystal.  adot is not ",          &
     &      "symmetric: ")')
         write(6,'("  adot(",i1,",",i1,") - adot(",i1,",",i1,") = ",     &
     &      e12.4)') 1,3,3,1,adot(1,3)-adot(3,1)

         stop

       endif

       if(abs(adot(2,3)-adot(3,2))/trace > EPS) then
         write(6,'("  Stopped in print crystal.  adot is not ",          &
     &      "symmetric: ")')
         write(6,'("  adot(",i1,",",i1,") - adot(",i1,",",i1,") = ",     &
     &      e12.4)') 2,3,3,2,adot(2,3)-adot(3,2)

         stop

       endif

       call adot_to_avec_sym(adot,avec,bvec)
       call adot_to_bdot(adot,vcell,bdot)     

       if(ipr /= 1) return

       write(6,*)
       write(6,*)
       write(6,*)   '      CRYSTAL DATA '
       write(6,*)
       write(6,*)
       write(6,'(3x,f16.8,5x," volume")') vcell

       write(6,*)
       write(6,*) '  real-space metric'
       write(6,*)
       write(6,'(3x ,3(1x,f16.8),5x," metric  g11,g12,g13 ")')           &
     &        adot(1,1),adot(1,2),adot(1,3)
       write(6,'(20x,2(1x,f16.8),5x," metric      g22,g23 ")')           &
     &                  adot(2,2),adot(2,3)
       write(6,'(38x,     f16.8 ,5x," metric          g33 ")')           &
     &                            adot(3,3)

!      gives lengths of the lattice vectors

       write(6,*)
       write(6,'(3(f16.8,1x),"  length 1,2,3  (a.u.)")')                 & 
     &        (sqrt(adot(i,i)),i=1,3)
       
!      gives angles between lattice vectors

       angl12 = (360.0/(2*PI))*                                          &
     &    acos(adot(1,2)/(sqrt(adot(1,1))*sqrt(adot(2,2))))
       angl13 = (360.0/(2*PI))*                                          &
     &    acos(adot(1,3)/(sqrt(adot(1,1))*sqrt(adot(3,3))))
       angl23 = (360.0/(2*PI))*                                          &
     &   acos(adot(2,3)/(sqrt(adot(2,2))*sqrt(adot(3,3))))
       write(6,*)
       write(6,'(3(f16.8,1x),"  angle 12,13,23 (degrees)")')             &
     &       angl12,angl13,angl23

!      position of atoms

       write(6,*)
       write(6,'(8x,"position (lattice coord.)",9x,                      &
     &        "position (cartesian coord. a.u.)",5x,"no. type")')
       write(6,*)

       ntt = 0
       do nt=1,ntype
         jmax = natom(nt)
         do ja=1,jmax
           ntt = ntt + 1
           car(1) = (avec(1,1)*rat(1,ja,nt) + avec(1,2)*rat(2,ja,nt) +   &
     &               avec(1,3)*rat(3,ja,nt))
           car(2) = (avec(2,1)*rat(1,ja,nt) + avec(2,2)*rat(2,ja,nt) +   &
     &               avec(2,3)*rat(3,ja,nt))
           car(3) = (avec(3,1)*rat(1,ja,nt) + avec(3,2)*rat(2,ja,nt) +   &
     &               avec(3,3)*rat(3,ja,nt))
           write(6,'(2x,3(2x,f9.5),3x,3(1x,e12.5),1x,i3,3x,a2,3x,        &
     &         " position")') (rat(i,ja,nt),i=1,3),                      &
     &              (car(i),i=1,3),ntt,nameat(nt)
         enddo
       enddo


       return

       end subroutine print_crystal
