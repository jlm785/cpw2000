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

!>     Generates a Maxwellian velocity distribution
!>     from ran0 random numbers of Numerical Recipes.
!>     uses simple pseudo-random numbers, do not
!>     use for complex statistics!

       subroutine move_random_v_atom(iseed,tempk,vat,vadot,              &
     & ntype,natom,atmass,adot,                                          &
     & mxdtyp,mxdatm)

!      Written November 2, 1993. jlm
!      Modified 25 July 1995. jlm
!      Modified 18 January 1999. jlm
!      Modified 18 September 2002. jlm
!      Modified 6 January 2017, f90. JLM
!      Modified, documentation, June 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins

!      version 4.94 of pw
!      version 1.5 of md

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !< array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !< array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  tempk                      !< ionic temperature (in Kelvin)

       integer, intent(in)                ::  ntype                      !< number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !< number of atoms of type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !< atomic mass of atoms of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in real space

!      input and output
 
       integer, intent(inout)             ::  iseed                      !< seed for rndom number generator

!      output

       real(REAL64), intent(out)          ::  vat(3,mxdatm,mxdtyp)       !< d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64), intent(out)          ::  vadot(3,3)                 !< d adot / d t  rate of change of metric

!      local variables

       real(REAL64)         ::  avec(3,3),bvec(3,3)
       real(REAL64)         ::  v1, v2, v3, fac
       real(REAL64)         ::  gasdev
       
       real(REAL64)         ::  tmass,ptot(3),ekint

       integer              ::  ntot

!      constants

       real(REAL64), parameter  ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
       real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  ::  EV = 27.2116_REAL64
       real(REAL64), parameter  ::  TAUTOK = 11604.9_REAL64 * EV

       integer, parameter       ::  IM = 2147483647
       integer, parameter       ::  MASK = 123459876

!      counters

       integer       ::  nt, i, j, k


       tmass = ZERO
       ptot(1) = ZERO
       ptot(2) = ZERO
       ptot(3) = ZERO

       call adot_to_avec(adot,avec,bvec)

!      you may comment the next line if you have problems with ieor


       iseed = ieor(iseed,MASK)
       if(iseed > IM) iseed = mod(iseed,IM)
       if(iseed == 0) iseed = MASK

       ntot = 0

       do nt = 1,ntype
         do j = 1,natom(nt)
 
!          random velocity

           do k = 1,3

             vat(k,j,nt) = gasdev(iseed)*                                &
     &                          sqrt(tempk/(atmass(nt)*TAUTOK))
             ptot(k) = ptot(k) + atmass(nt)*vat(k,j,nt)

           enddo

           tmass = tmass + atmass(nt)
           ntot = ntot + 1

         enddo
       enddo

!      total linear momentum is zero

       ekint = ZERO
       do nt = 1,ntype
         do j = 1,natom(nt)
           do k=1,3
             vat(k,j,nt) = vat(k,j,nt) - ptot(k)/tmass
             ekint = ekint + atmass(nt)*vat(k,j,nt)*vat(k,j,nt)/2
           enddo
         enddo
       enddo

       if(ntot > 1) then
         fac = (ekint*2) / (3*ntot-3)*UM
         fac = sqrt(tempk/(fac*TAUTOK))
       else
         fac = UM
       endif

!      changes to lattice coordinates and corrects by fac

       do nt = 1,ntype
         do j = 1,natom(nt)

           v1 = vat(1,j,nt)*bvec(1,1) + vat(2,j,nt)*bvec(2,1) +          &
     &          vat(3,j,nt)*bvec(3,1)
           v2 = vat(1,j,nt)*bvec(1,2) + vat(2,j,nt)*bvec(2,2) +          &
     &          vat(3,j,nt)*bvec(3,2)
           v3 = vat(1,j,nt)*bvec(1,3) + vat(2,j,nt)*bvec(2,3) +          &
     &          vat(3,j,nt)*bvec(3,3)

           vat(1,j,nt) = v1*fac/(2*PI)
           vat(2,j,nt) = v2*fac/(2*PI)
           vat(3,j,nt) = v3*fac/(2*PI)

         enddo
       enddo

!      vadot is initialized to zero

       do i = 1,3
       do j = 1,3
         vadot(i,j) = ZERO
       enddo
       enddo

       iseed = ieor(iseed,MASK)

       return
       end subroutine move_random_v_atom


