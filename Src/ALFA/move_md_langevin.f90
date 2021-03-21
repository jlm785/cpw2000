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

!>     Performs a canonical molecular dynamics time step
!>     using Langevin molecular dynamics
!>     R. biswas and D. R. Hamann, Phys Rev B 34, 895 (1986),
!>     but using the integration method of
!>     D. Beeman, J. Comput. Phys. 20, 130 (76). Eqs 10a and 10c,
!>     for beta=0 we have a microcanonical Verlet algorithm
!>     as the integration method is different, the statistical
!>     properties should be rechecked....
!>     The previous estimate of the velocity is used for the friction force

       subroutine move_md_langevin(rat,vat,force,istmd,tstep,ekin,       &
     & beta,tempk,iseed,                                                 &
     & rat1,frc1,                                                        &
     & ntype,natom,atmass,adot,                                          &
     & mxdatm,mxdtyp)

!      Written 23 october 93. jlm
!      Modified 25 july 1995. jlm
!      Modified 18 january 1999. jlm
!      Modified 19 november 2002 (save). jlm
!      Modified 6 January 2017, f90. JLM
!      Modified 6 June 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins / Nadia Binggeli

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

       real(REAL64), intent(in)           ::  tempk                      !< ionic temperature (in Kelvin)
       real(REAL64), intent(in)           ::  beta                       !< friction coefficient/mass (in a.u.)

       integer, intent(in)                ::  ntype                      !< number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !< number of atoms of type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !< atomic mass of atoms of type i
       real(REAL64), intent(in)           ::  adot(3,3)                  !< metric in real space

!      input and output
       
       real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)       !< lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  vat(3,mxdatm,mxdtyp)       !< d rat / d t  velocity in lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  rat1(3,mxdatm,mxdtyp)      !< previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  frc1(3,mxdatm,mxdtyp)      !< previous value of force on atom j of type i

       integer, intent(inout)             ::  iseed                      !< seed for rndom number generator

!      output

       real(REAL64), intent(out)          ::  ekin                       !< kinetic energy of the nuclei (atomic units)

!      local variables

       real(REAL64)         ::  rat0(3), vat1(3), frc0(3)
       real(REAL64)         ::  facsig
       real(REAL64)         ::  avec(3,3),bvec(3,3)
       real(REAL64)         ::  gasdev, sigma(3), frcsig(3)

!      constants

       real(REAL64), parameter  ::  ZERO = 0.0_REAL64
       real(REAL64), parameter  ::  PI = 3.14159265358979323846_REAL64
       real(REAL64), parameter  ::  EV = 27.2116_REAL64
       real(REAL64), parameter  ::  TAUTOK = 11604.9_REAL64 * EV

       integer, parameter       ::  IM = 2147483647
       integer, parameter       ::  MASK = 123459876

!      counters

       integer       ::  nt, j, k

       ekin = ZERO

       call adot_to_avec(adot,avec,bvec)

!      you may comment the next line (and the last line)
!      if you have problems with ieor

       iseed = ieor(iseed,MASK)
       if(iseed > IM) iseed = mod(iseed,IM)
       if(iseed == 0) iseed = MASK

       if(istmd == 1) then
         do nt = 1,ntype
           facsig = sqrt(2*atmass(nt)*beta*tempk / (TAUTOK*tstep))
           do j = 1,natom(nt)

             rat1(1,j,nt) = rat(1,j,nt)
             rat1(2,j,nt) = rat(2,j,nt)
             rat1(3,j,nt) = rat(3,j,nt)

             vat1(1) = vat(1,j,nt)
             vat1(2) = vat(2,j,nt)
             vat1(3) = vat(3,j,nt)

!            random force

             do k = 1,3
               sigma(k) = gasdev(iseed)
             enddo

             frcsig(1) = (sigma(1)*bvec(1,1) + sigma(2)*bvec(2,1) +      &
     &                  sigma(3)*bvec(3,1)) * facsig / (2*PI)
             frcsig(2) = (sigma(1)*bvec(1,2) + sigma(2)*bvec(2,2) +      &
     &                  sigma(3)*bvec(3,2)) * facsig / (2*PI)
             frcsig(3) = (sigma(1)*bvec(1,3) + sigma(2)*bvec(2,3) +      &
     &                  sigma(3)*bvec(3,3)) * facsig / (2*PI)

             frc1(1,j,nt) = force(1,j,nt)                                &
     &                    - beta*atmass(nt)*vat1(1) + frcsig(1)
             frc1(2,j,nt) = force(2,j,nt)                                &
     &                    - beta*atmass(nt)*vat1(2) + frcsig(2)
             frc1(3,j,nt) = force(3,j,nt)                                &
     &                    - beta*atmass(nt)*vat1(3) + frcsig(3)

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
           facsig = sqrt(2*atmass(nt)*beta*tempk / (TAUTOK*tstep))
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

!            random force

             do k = 1,3
               sigma(k) = gasdev(iseed)
             enddo

             frcsig(1) = (sigma(1)*bvec(1,1) + sigma(2)*bvec(2,1) +      &
     &                  sigma(3)*bvec(3,1)) * facsig / (2*PI)
             frcsig(2) = (sigma(1)*bvec(1,2) + sigma(2)*bvec(2,2) +      &
     &                  sigma(3)*bvec(3,2)) * facsig / (2*PI)
             frcsig(3) = (sigma(1)*bvec(1,3) + sigma(2)*bvec(2,3) +      &
     &                  sigma(3)*bvec(3,3)) * facsig / (2*PI)

             frc1(1,j,nt) = force(1,j,nt)                                &
     &                    - beta*atmass(nt)*vat(1,j,nt) + frcsig(1)
             frc1(2,j,nt) = force(2,j,nt)                                &
     &                    - beta*atmass(nt)*vat(2,j,nt) + frcsig(2)
             frc1(3,j,nt) = force(3,j,nt)                                &
     &                    - beta*atmass(nt)*vat(3,j,nt) + frcsig(3)

             vat1(1) = (rat1(1,j,nt) - rat0(1)) / tstep +                &
     &         tstep * (2*frc1(1,j,nt) + frc0(1))/(6*atmass(nt))
             vat1(2) = (rat1(2,j,nt) - rat0(2)) / tstep +                &
     &         tstep * (2*frc1(2,j,nt) + frc0(2))/(6*atmass(nt))
             vat1(3) = (rat1(3,j,nt) - rat0(3)) / tstep +                &
     &         tstep * (2*frc1(3,j,nt) + frc0(3))/(6*atmass(nt))

             rat(1,j,nt) = rat1(1,j,nt) + tstep*(vat1(1) + tstep*        &
     &         ((4*frc1(1,j,nt)-frc0(1))/3)/(2*atmass(nt)))
             rat(2,j,nt) = rat1(2,j,nt) + tstep*(vat1(2) + tstep*        &
     &         ((4*frc1(2,j,nt)-frc0(2))/3)/(2*atmass(nt)))
             rat(3,j,nt) = rat1(3,j,nt) + tstep*(vat1(3) + tstep*        &
     &         ((4*frc1(3,j,nt)-frc0(3))/3)/(2*atmass(nt)))

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

       iseed = ieor(iseed,MASK)

       return
       end subroutine move_md_langevin
