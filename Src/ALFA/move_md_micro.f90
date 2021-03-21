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

!      Performs a microcanonical molecular dynamics time step/
!      D. Beeman, J. Comput. Phys. 20, 130 (76). Eqs 10a and 10c,
!      the trajectory is the same as in the Verlet algorithm,
!      only the velocity estimate is improved.

       subroutine move_md_micro(rat,vat,force,istmd,tstep,ekin,          &
     & rat1,frc1,                                                        &
     & ntype,natom,atmass,adot,                                          &
     & mxdatm,mxdtyp)

!      Written 23 October 93. jlm
!      Modified 20 January 1999. jlm
!      Modified 6 January 2017, f90. JLM
!      Modified, documentation, June 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins/Renata Wentzcovitch

!      version 4.94 of pw
!      version 1.5 of md

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !< array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !< array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  force(3,mxdatm,mxdtyp)     !< d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)

       integer,intent(in)                 ::  istmd                      !< md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(in)           ::  tstep                      !< molecular dynamics time step (in a.u.)

       integer, intent(in)                ::  ntype                      !< number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !< number of atoms of type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !< atomic mass of atoms of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in real space

!      input and output
       
       real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)       !< lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  vat(3,mxdatm,mxdtyp)       !< d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  rat1(3,mxdatm,mxdtyp)      !< previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  frc1(3,mxdatm,mxdtyp)      !< previous value of force on atom j of type i

!      output

       real(REAL64), intent(out)          ::  ekin                       !< kinetic energy of the nuclei (atomic units)

!      local variables

       real(REAL64)         ::  rat0(3), vat1(3), frc0(3)    

!      constants

       real(REAL64), parameter  :: ZERO = 0.0_REAL64

!      counters

       integer       ::  nt, j

       ekin = ZERO

       if(istmd == 1) then

         do nt = 1,ntype
           do j = 1,natom(nt)

             rat1(1,j,nt) = rat(1,j,nt)
             rat1(2,j,nt) = rat(2,j,nt)
             rat1(3,j,nt) = rat(3,j,nt)

             vat1(1) = vat(1,j,nt)
             vat1(2) = vat(2,j,nt)
             vat1(3) = vat(3,j,nt)

             frc1(1,j,nt) = force(1,j,nt)
             frc1(2,j,nt) = force(2,j,nt)
             frc1(3,j,nt) = force(3,j,nt)

             rat(1,j,nt) = rat1(1,j,nt) + tstep*(vat1(1) +               &
     &              tstep*frc1(1,j,nt)/(2*atmass(nt)))
             rat(2,j,nt) = rat1(2,j,nt) + tstep*(vat1(2) +               &
     &              tstep*frc1(2,j,nt)/(2*atmass(nt)))
             rat(3,j,nt) = rat1(3,j,nt) + tstep*(vat1(3) +               &
     &              tstep*frc1(3,j,nt)/(2*atmass(nt)))

             ekin = ekin + (atmass(nt)/2)*                               &
     &                         (vat1(1)*adot(1,1)*vat1(1) +              &
     &                          vat1(1)*adot(1,2)*vat1(2) +              &
     &                          vat1(1)*adot(1,3)*vat1(3) +              &
     &                          vat1(2)*adot(2,1)*vat1(1) +              &
     &                          vat1(2)*adot(2,2)*vat1(2) +              &
     &                          vat1(2)*adot(2,3)*vat1(3) +              &
     &                          vat1(3)*adot(3,1)*vat1(1) +              &
     &                          vat1(3)*adot(3,2)*vat1(2) +              &
     &                          vat1(3)*adot(3,3)*vat1(3))

             vat(1,j,nt) = vat1(1) + tstep*frc1(1,j,nt)/atmass(nt)
             vat(2,j,nt) = vat1(2) + tstep*frc1(2,j,nt)/atmass(nt)
             vat(3,j,nt) = vat1(3) + tstep*frc1(3,j,nt)/atmass(nt)

           enddo
         enddo

       else

         do nt = 1,ntype
           do j = 1,natom(nt)

             rat0(1) = rat1(1,j,nt)
             rat0(2) = rat1(2,j,nt)
             rat0(3) = rat1(3,j,nt)

             rat1(1,j,nt) = rat(1,j,nt)
             rat1(2,j,nt) = rat(2,j,nt)
             rat1(3,j,nt) = rat(3,j,nt)

             frc0(1) = frc1(1,j,nt)
             frc0(2) = frc1(2,j,nt)
             frc0(3) = frc1(3,j,nt)

             frc1(1,j,nt) = force(1,j,nt)
             frc1(2,j,nt) = force(2,j,nt)
             frc1(3,j,nt) = force(3,j,nt)

             vat1(1) = (rat1(1,j,nt) - rat0(1)) / tstep +                &
     &         tstep * (2*frc1(1,j,nt) + frc0(1))/(6*atmass(nt))
             vat1(2) = (rat1(2,j,nt) - rat0(2)) / tstep +                &
     &         tstep * (2*frc1(2,j,nt) + frc0(2))/(6*atmass(nt))
             vat1(3) = (rat1(3,j,nt) - rat0(3)) / tstep +                &
     &         tstep * (2*frc1(3,j,nt) + frc0(3))/(6*atmass(nt))

             rat(1,j,nt) = 2*rat1(1,j,nt) - rat0(1) +                    &
     &           tstep*tstep*frc1(1,j,nt)/atmass(nt)
             rat(2,j,nt) = 2*rat1(2,j,nt) - rat0(2) +                    &
     &           tstep*tstep*frc1(2,j,nt)/atmass(nt)
             rat(3,j,nt) = 2*rat1(3,j,nt) - rat0(3) +                    &
     &           tstep*tstep*frc1(3,j,nt)/atmass(nt)

             ekin = ekin + (atmass(nt)/2)*                               &
     &                            (vat1(1)*adot(1,1)*vat1(1) +           &
     &                             vat1(1)*adot(1,2)*vat1(2) +           &
     &                             vat1(1)*adot(1,3)*vat1(3) +           &
     &                             vat1(2)*adot(2,1)*vat1(1) +           &
     &                             vat1(2)*adot(2,2)*vat1(2) +           &
     &                             vat1(2)*adot(2,3)*vat1(3) +           &
     &                             vat1(3)*adot(3,1)*vat1(1) +           &
     &                             vat1(3)*adot(3,2)*vat1(2) +           &
     &                             vat1(3)*adot(3,3)*vat1(3))

             vat(1,j,nt) = vat1(1) + tstep*                              &
     &             (3*frc1(1,j,nt) - frc0(1))/(2*atmass(nt))
             vat(2,j,nt) = vat1(2) + tstep*                              &
     &             (3*frc1(2,j,nt) - frc0(2))/(2*atmass(nt))
             vat(3,j,nt) = vat1(3) + tstep*                              &
     &             (3*frc1(3,j,nt) - frc0(3))/(2*atmass(nt))

           enddo
         enddo

       endif

       return
       end subroutine move_md_micro
