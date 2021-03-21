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

!>     given a metric adot, returns a matrix
!>     avec={a1,a2,a3}  such that adot=transpose(avec)*avec
!>     within the tolerance specified by tol
!>     avec has the highest possible symmetry compatible
!>     with tol...
!>     It also returns the primitive and conventional lattice
!>     vectors of the Niggli cell.

       subroutine adot_to_avec_aconv(adot,avec,bvec,aconv,avecnig)

      
!      Written April 2000. jlm
!      Modified, f90, metric_ident, 4 June 2014. JLM
!      Modified, documentation, August 2019. JLM
!      copyright INESC-MN/Jose Luis Martins


!      Version 4.94

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       real(REAL64), intent(in)           ::  adot(3,3)                  !<  metric in direct space in atomic units (Bohr radius)

!      output

       real(REAL64), intent(out)          ::  avec(3,3)                  !<  primitive lattice vectors that generate adot in canonical orientation
       real(REAL64), intent(out)          ::  bvec(3,3)                  !<  reciprocal lattice vectors
       real(REAL64), intent(out)          ::  aconv(3,3)                 !<  conventional lattice vextors of the Niggli cell
       real(REAL64), intent(out)          ::  avecnig(3,3)               !<  primitive lattice vextors of the Niggli cell

!      local variables

       integer           ::  bravais,mtotal(3,3),ipr
       real(REAL64)      ::  vcell,adotnig(3,3),adotsym(3,3),tol

!      counters

       integer       ::  j

!      parameters

       real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter :: EPS = 0.0000000001_REAL64




!      identify the lattice and find niggli vectors

       tol = 1.0d-9
       ipr = 0
       call metric_ident(adot,adotnig,adotsym,bravais,mtotal,            &
     &                   avec,aconv,avecnig,tol,ipr)

!      replaces by symmetrized adot

!       do i=1,3
!       do j=1,3
!         adot(j,i) = adotsym(j,i)
!       enddo
!       enddo

!      compute the reciprocal vectors and cell volume 

       bvec(1,1) = avec(2,2)*avec(3,3) - avec(3,2)*avec(2,3)
       bvec(2,1) = avec(3,2)*avec(1,3) - avec(1,2)*avec(3,3)
       bvec(3,1) = avec(1,2)*avec(2,3) - avec(2,2)*avec(1,3)
       bvec(1,2) = avec(2,3)*avec(3,1) - avec(3,3)*avec(2,1)
       bvec(2,2) = avec(3,3)*avec(1,1) - avec(1,3)*avec(3,1)
       bvec(3,2) = avec(1,3)*avec(2,1) - avec(2,3)*avec(1,1)
       bvec(1,3) = avec(2,1)*avec(3,2) - avec(3,1)*avec(2,2)
       bvec(2,3) = avec(3,1)*avec(1,2) - avec(1,1)*avec(3,2)
       bvec(3,3) = avec(1,1)*avec(2,2) - avec(2,1)*avec(1,2)

!      cell volume

       vcell = bvec(1,1)*avec(1,1) + bvec(2,1)*avec(2,1) +               &
     &         bvec(3,1)*avec(3,1)

       if(vcell < EPS) then
         write(6,'("    STOPPED in adot_to_avec_sym.    Cell volume ",   &
     &      "squared= ",e12.4)') vcell

         stop

       endif

       do j=1,3
         bvec(1,j) = 2*PI*bvec(1,j)/vcell
         bvec(2,j) = 2*PI*bvec(2,j)/vcell
         bvec(3,j) = 2*PI*bvec(3,j)/vcell
       enddo

       return

       end subroutine adot_to_avec_aconv
