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

!>     calculates the errors relevant for structural relaxation
!>     convergence

       subroutine force_stress_error(energy,force,stress,press,strext,   &
     & deltentpy,errfrc,errstr,minrat,minstr,                            &
     & adot,ntype,natom,                                                 &
     & mxdtyp,mxdatm)

!      Written 6 August 2019. JLM
!      Based on move_lbfgs_call.f90
!      copyright INESC-MN/Jose Luis Martins

!      version 4.94


       implicit none

       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !<  array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !<  array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  energy                     !<  energy
       real(REAL64), intent(in)           ::  force(3,mxdatm,mxdtyp)     !<  k-th component (in contravariant lattice coordinates)  of the force of the n-th atom of type i
       real(REAL64), intent(in)           ::  stress(3,3)                !<  stress tensor (in contravariant lattice coordinates)

       real(REAL64), intent(in)           ::  press                      !<  applied pressure (hartree/au)
       real(REAL64), intent(in)           ::  strext(3,3)                !<  external (applied) stress tensor in lattice coordinates * vcell (hartree/au)

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space
       integer, intent(in)                ::  ntype                      !<  number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !<  number of atoms of type i

       integer, intent(in)                ::  minrat                     !<  if =1 minimize with respect to atomic positions
       integer, intent(in)                ::  minstr                     !<  if =1 minimize with respect to all adot variables, if =2 minimize with respect to adot(3,3)

!      output

       real(REAL64), intent(out)          ::  deltentpy                  !<  enthalpy difference from last iteration
       real(REAL64), intent(out)          ::  errfrc                     !<  maximum error in force (cartesian coordiantes) Hartree/Bohr
       real(REAL64), intent(out)          ::  errstr                     !<  maximum error in stress

!      local variables

       integer, save           ::  ifirst = 0
       real(REAL64), save      ::  oldentpy
       real(REAL64)            ::  car, strimb(3,3), gdav
       real(REAL64)            ::  entpy

       real(REAL64)            :: adinv(3,3),vcell
       real(REAL64)            :: avec(3,3),bvec(3,3)

!      parameters

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: TWOPI = 2.0_REAL64*PI       

!      counters

       integer       ::  i, j, k, l
       integer       ::  nt, ja


       call adot_to_bdot(adot,vcell,adinv)
       call adot_to_avec(adot,avec,bvec)

       do i=1,3
       do j=1,3
         adinv(i,j) = adinv(i,j) / (TWOPI*TWOPI)
       enddo
       enddo

       if (minstr == 2) then

         entpy = energy

       else

         entpy = energy + press*vcell  +                                 &
     &          (strext(1,1)*adot(1,1) +                                 &
     &           strext(1,2)*adot(2,1) +                                 &
     &           strext(1,3)*adot(3,1) +                                 &
     &           strext(2,1)*adot(1,2) +                                 &
     &           strext(2,2)*adot(2,2) +                                 &
     &           strext(2,3)*adot(3,2) +                                 &
     &           strext(3,1)*adot(1,3) +                                 &
     &           strext(3,2)*adot(2,3) +                                 &
     &           strext(3,3)*adot(3,3)) / 2

       endif

       if(ifirst == 0) then
         deltentpy = ZERO
         oldentpy = entpy
       else
         deltentpy = entpy - oldentpy
         oldentpy = entpy
       endif
       
       ifirst = 1

       errfrc = ZERO

       if(minrat == 1) then
       
         do nt=1,ntype
           do ja=1,natom(nt)
             car = force(1,ja,nt)*adot(1,1)*force(1,ja,nt) +             &
     &             force(1,ja,nt)*adot(1,2)*force(2,ja,nt) +             &
     &             force(1,ja,nt)*adot(1,3)*force(3,ja,nt) +             &
     &             force(2,ja,nt)*adot(2,1)*force(1,ja,nt) +             &
     &             force(2,ja,nt)*adot(2,2)*force(2,ja,nt) +             &
     &             force(2,ja,nt)*adot(2,3)*force(3,ja,nt) +             &
     &             force(3,ja,nt)*adot(3,1)*force(1,ja,nt) +             &
     &             force(3,ja,nt)*adot(3,2)*force(2,ja,nt) +             &
     &             force(3,ja,nt)*adot(3,3)*force(3,ja,nt)
             if(abs(car) > errfrc) errfrc = abs(car)
           enddo
         enddo
         errfrc = sqrt(errfrc)
        
       endif

       errstr = ZERO

       if(minstr == 1) then

         do i = 1,3
         do j = 1,3
           strimb(i,j) = strext(i,j)                                     &
     &                 + press*vcell*adinv(i,j)-stress(i,j)
         enddo
         enddo

         do i=1,3
         do j=1,3
           gdav = ZERO
           do k=1,3
           do l=1,3
             gdav = gdav - avec(i,k)*strimb(k,l)*avec(j,l)
           enddo
           enddo
           if(errstr < abs(gdav)) errstr = abs(gdav)
         enddo
         enddo

       elseif(minstr == 2) then

         errstr = abs(stress(3,3)*adot(3,3))

       endif
       
       return
       end subroutine force_stress_error
 
