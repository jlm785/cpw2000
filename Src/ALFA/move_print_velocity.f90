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

!>     Prints the atomic velocity (and fictitious cell velocity)

       subroutine move_print_velocity(flgcal,                            &
     & adot,ntype,natom,nameat,vat,vadot,                                &
     & mxdtyp,mxdatm)

!      Written 18 September 2002. jlm
!      Modified 6 January 2017, f90. JLM
!      Modified, documentation, Jume 2019. JLM
!      Modified, EPILNG, August 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins

!      version 4.94 cpw
!      version 1.5 of md


       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       character(len=6), intent(in)       ::  flgcal                     !<  type of md calculation

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in real space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i
       character(len=2), intent(in)       ::  nameat(mxdtyp)             !<  chemical symbol for the type i

       real(REAL64), intent(in)           ::  vat(3,mxdatm,mxdtyp)       !<  d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64), intent(in)           ::  vadot(3,3)                 !<  d adot / d t  rate of change of metric
      
!      local variables

       real(REAL64)  ::    avec(3,3),bvec(3,3)
       real(REAL64)  ::    car(3)

!      counters

       integer       ::  i, ja, nt, ntt


!      checks symmetry and positiveness of adot

       call adot_to_avec_sym(adot,avec,bvec)

       if(flgcal == 'VCSLNG' .or. flgcal == 'VCSMIC') then

         write(6,*)
         write(6,'("  Time derivative of real-space Metric")')
         write(6,*)
         write(6,'(3x ,3(1x,f16.8),5x," dmet/dt  g11,g12,g13 ")')        &
     &                vadot(1,1),vadot(1,2),vadot(1,3)
         write(6,'(20x,2(1x,f16.8),5x," dmet/dt      g22,g23 ")')        &
     &                           vadot(2,2),vadot(2,3)
         write(6,'(38x,     f16.8 ,5x," dmet/dt          g33 ")')        &
     &                                      vadot(3,3)
         write(6,*)

       endif

       if(flgcal == 'EPILNG' .or. flgcal == 'VCSMIC') then

         write(6,*)
         write(6,'("  Time derivative of real-space Metric")')
         write(6,*)
         write(6,'(38x,     f16.8 ,5x," dmet/dt          g33 ")')        &
     &                                      vadot(3,3)
         write(6,*)

       endif

       if(flgcal == 'VCSLNG' .or. flgcal == 'MICRO ' .or.                &
     &    flgcal == 'LANG  ' .or. flgcal == 'VCSMIC' .or.                &
     &    flgcal == 'EPILNG') then

         write(6,*)
         write(6,'(8x,"Velocity (Lattice coord.)",20x,                   &
     &        "Velocity (Cartesian coord. m/s.)",5x,"No. type")')
         write(6,*)

         ntt = 0
         do nt = 1,ntype
           do ja = 1,natom(nt)
             ntt = ntt + 1
             car(1) = (avec(1,1)*vat(1,ja,nt) +                          &
     &                 avec(1,2)*vat(2,ja,nt) +                          &
     &                 avec(1,3)*vat(3,ja,nt))*2.19d6
             car(2) = (avec(2,1)*vat(1,ja,nt) +                          &
     &                 avec(2,2)*vat(2,ja,nt) +                          &
     &                 avec(2,3)*vat(3,ja,nt))*2.19d6
             car(3) = (avec(3,1)*vat(1,ja,nt) +                          &
     &                 avec(3,2)*vat(2,ja,nt) +                          &
     &                 avec(3,3)*vat(3,ja,nt))*2.19d6
             write(6,'(2x,3(2x,f12.8),3x,3(2x,f11.3),1x,i3,3x,a2,3x,     &
     &            " velocity")')                                         &
     &             (vat(i,ja,nt),i=1,3),(car(i),i=1,3),ntt,nameat(nt)
           enddo
         enddo
         write(6,*)

       endif

       return
       end subroutine move_print_velocity
