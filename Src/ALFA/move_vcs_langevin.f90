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

!>     Performs a Langevin molecular dynamics time step
!>     for constant applied pressure and/or thermodynamic tension
!>     I. Souza and J. L. Martins, Phys Rev B 55, 8733 (1997).
!>     R. Biswas and D. R. Hamann, Phys Rev B 34, 895 (1986),
!>     but using the integration method of
!>     D. Beeman, J. Comput. Phys. 20, 130 (76). Eqs 10a and 10c,
!>     the cell temperature is zero

       subroutine move_vcs_langevin(rat,vat,adot,vadot,force,stress,     &
     & istmd,tstep,ekin,ekcell,strext,press,celmas,                      &
     & beta,tempk,iseed,                                                 &
     & rat1,frc1,adot1,frcel1,                                           &
     & ntype,natom,atmass,                                               &
     & mxdatm,mxdtyp)


!      Written 20 January 99. Based on Ivo Souza subroutines. JLM
!      Modified 8 March 2000. jlm
!      Modified 19 November 2020 (save). jlm
!      Modified 7 January 2017, f90. JLM
!      Modified, documentation, June 2019. JLM
!      Copyright inesc-mn/Jose Luis Martins/Ivo Souza

!      version 4.94 of pw
!      version 1.5 of md

       implicit none

       integer, parameter :: REAL64 = selected_real_kind(12)

!      input

       integer, intent(in)                ::  mxdtyp                     !< array dimension of types of atoms
       integer, intent(in)                ::  mxdatm                     !< array dimension of number of atoms of a given type

       real(REAL64), intent(in)           ::  force(3,mxdatm,mxdtyp)     !< d energy / d rat,  force on the n-th atom of type i (contravarian components, hartree/bohr)
       real(REAL64), intent(in)           ::  stress(3,3)                !< d energy / d adot,  stress tensor (contravariant)

       integer,intent(in)                 ::  istmd                      !< md step. Equal to 1 in first step of molecular dynamics
       real(REAL64), intent(in)           ::  tstep                      !< molecular dynamics time step (in a.u.)

       real(REAL64), intent(in)           ::  tempk                      !< ionic temperature (in Kelvin)
       real(REAL64), intent(in)           ::  beta                       !< friction coefficient/mass (in a.u.)

       real(REAL64), intent(in)           ::  strext(3,3)                !< external applied stress
       real(REAL64), intent(in)           ::  press                      !< external pressure
       real(real64), intent(in)           ::  celmas                     !< fictitious cell mass       

       integer, intent(in)                ::  ntype                      !< number of types of atoms
       integer, intent(in)                ::  natom(mxdtyp)              !< number of atoms of type i
       real(REAL64), intent(in)           ::  atmass(mxdtyp)             !< atomic mass of atoms of type i

!      input and output
       
       real(REAL64), intent(inout)        ::  rat(3,mxdatm,mxdtyp)       !< lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  vat(3,mxdatm,mxdtyp)       !< d rat / d t  velocity in lattice coordinates of atom j of type i

       real(REAL64), intent(inout)        ::  adot(3,3)                  !< metric in real space
       real(REAL64), intent(inout)        ::  vadot(3,3)                 !< d adot / d t  rate of change of metric

       integer, intent(inout)             ::  iseed                      !< seed for rndom number generator

       real(REAL64), intent(inout)        ::  rat1(3,mxdatm,mxdtyp)      !< previous value of lattice coordinates of atom j of type i
       real(REAL64), intent(inout)        ::  frc1(3,mxdatm,mxdtyp)      !< previous value of force on atom j of type i
       real(REAL64), intent(inout)        ::  adot1(3,3)                 !< previous value of adot
       real(REAL64), intent(inout)        ::  frcel1(3,3)                !< previous cell "force" (covariant components)

!      output

       real(REAL64), intent(out)          ::  ekin                       !< kinetic energy of the nuclei (atomic units)
       real(REAL64), intent(out)          ::  ekcell                     !< kinetic energy in Hartree of cell (fictitious)C

!      local variables

       real(REAL64)         ::  rat0(3), vat1(3), frc0(3)
       real(REAL64)         ::  vcell,bdot(3,3)
       real(REAL64)         ::  adotm1(3,3),adm1vad(3,3)
       real(REAL64)         ::  strkin(3,3)
       real(REAL64)         ::  strcont(3,3), strmix(3,3), strcov(3,3)
       real(REAL64)         ::  adot0(3,3), vadot1(3,3)
       real(REAL64)         ::  frcel0(3,3),frcell(3,3)
       real(REAL64)         ::  gggg, tram1v

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

       integer       ::  nt, i, j, k


       ekin = ZERO

       call adot_to_bdot(adot,vcell,bdot)
       call adot_to_avec(adot,avec,bvec)

!      you may comment the next line (and the last line)
!      if you have problems with ieor

       iseed = ieor(iseed,MASK)
       if(iseed > IM) iseed = mod(iseed,IM)
       if(iseed == 0) iseed = MASK

       do j=1,3
         adotm1(j,1) = bdot(j,1)/ (4*PI*PI)
         adotm1(j,2) = bdot(j,2)/ (4*PI*PI)
         adotm1(j,3) = bdot(j,3)/ (4*PI*PI)
       enddo

       do j=1,3
       do i=1,3
         adm1vad(i,j) = adotm1(i,1)*vadot(1,j) +                         &
     &                  adotm1(i,2)*vadot(2,j) +                         &
     &                  adotm1(i,3)*vadot(3,j)
         strkin(i,j) = ZERO
       enddo
       enddo

!      motion of atoms

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
     &                    - beta*atmass(nt)*vat1(1) + frcsig(1)          &
     &                    - atmass(nt)*(adm1vad(1,1)*vat1(1)+            &
     &                                  adm1vad(1,2)*vat1(2)+            &
     &                                  adm1vad(1,3)*vat1(3))
             frc1(2,j,nt) = force(2,j,nt)                                &
     &                    - beta*atmass(nt)*vat1(2) + frcsig(2)          &
     &                    - atmass(nt)*(adm1vad(2,1)*vat1(1)+            &
     &                                  adm1vad(2,2)*vat1(2)+            &
     &                                  adm1vad(2,3)*vat1(3))
             frc1(3,j,nt) = force(3,j,nt)                                &
     &                    - beta*atmass(nt)*vat1(3) + frcsig(3)          &
     &                    - atmass(nt)*(adm1vad(3,1)*vat1(1)+            &
     &                                  adm1vad(3,2)*vat1(2)+            &
     &                                  adm1vad(3,3)*vat1(3))

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


             strkin(1,1) = strkin(1,1) + vat1(1)*vat1(1)*atmass(nt)
             strkin(2,1) = strkin(2,1) + vat1(2)*vat1(1)*atmass(nt)
             strkin(3,1) = strkin(3,1) + vat1(3)*vat1(1)*atmass(nt)
             strkin(1,2) = strkin(1,2) + vat1(1)*vat1(2)*atmass(nt)
             strkin(2,2) = strkin(2,2) + vat1(2)*vat1(2)*atmass(nt)
             strkin(3,2) = strkin(3,2) + vat1(3)*vat1(2)*atmass(nt)
             strkin(1,3) = strkin(1,3) + vat1(1)*vat1(3)*atmass(nt)
             strkin(2,3) = strkin(2,3) + vat1(2)*vat1(3)*atmass(nt)
             strkin(3,3) = strkin(3,3) + vat1(3)*vat1(3)*atmass(nt)

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
     &                    - beta*atmass(nt)*vat(1,j,nt) + frcsig(1)      &
     &                    - atmass(nt)*(adm1vad(1,1)*vat(1,j,nt)+        &
     &                                  adm1vad(1,2)*vat(2,j,nt)+        &
     &                                  adm1vad(1,3)*vat(3,j,nt))
             frc1(2,j,nt) = force(2,j,nt)                                &
     &                    - beta*atmass(nt)*vat(2,j,nt) + frcsig(2)      &
     &                    - atmass(nt)*(adm1vad(2,1)*vat(1,j,nt)+        &
     &                                  adm1vad(2,2)*vat(2,j,nt)+        &
     &                                  adm1vad(2,3)*vat(3,j,nt))
             frc1(3,j,nt) = force(3,j,nt)                                &
     &                    - beta*atmass(nt)*vat(3,j,nt) + frcsig(3)      &
     &                    - atmass(nt)*(adm1vad(3,1)*vat(1,j,nt)+        &
     &                                  adm1vad(3,2)*vat(2,j,nt)+        &
     &                                  adm1vad(3,3)*vat(3,j,nt))


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

             strkin(1,1) = strkin(1,1) + vat1(1)*vat1(1)*atmass(nt)
             strkin(2,1) = strkin(2,1) + vat1(2)*vat1(1)*atmass(nt)
             strkin(3,1) = strkin(3,1) + vat1(3)*vat1(1)*atmass(nt)
             strkin(1,2) = strkin(1,2) + vat1(1)*vat1(2)*atmass(nt)
             strkin(2,2) = strkin(2,2) + vat1(2)*vat1(2)*atmass(nt)
             strkin(3,2) = strkin(3,2) + vat1(3)*vat1(2)*atmass(nt)
             strkin(1,3) = strkin(1,3) + vat1(1)*vat1(3)*atmass(nt)
             strkin(2,3) = strkin(2,3) + vat1(2)*vat1(3)*atmass(nt)
             strkin(3,3) = strkin(3,3) + vat1(3)*vat1(3)*atmass(nt)

             vat(1,j,nt) = vat1(1) + tstep*                              &
     &             (3*frc1(1,j,nt) - frc0(1))/(2*atmass(nt))
             vat(2,j,nt) = vat1(2) + tstep*                              &
     &             (3*frc1(2,j,nt) - frc0(2))/(2*atmass(nt))
             vat(3,j,nt) = vat1(3) + tstep*                              &
     &             (3*frc1(3,j,nt) - frc0(3))/(2*atmass(nt))

           enddo
         enddo
       endif

       iseed = ieor(iseed,mask)

!      cell motion

       do j = 1,3
       do i = 1,3
         strcont(i,j) = stress(i,j) + strkin(i,j) - strext(i,j) -        &
     &                 press*vcell*adotm1(i,j)
       enddo
       enddo

       do j = 1,3
       do i = 1,3
         strmix(i,j) = adot(i,1)*strcont(1,j) +                          &
     &                 adot(i,2)*strcont(2,j) +                          &
     &                 adot(i,3)*strcont(3,j)
       enddo
       enddo

       do j = 1,3
       do i = 1,3
         strcov(i,j) = strmix(i,1)*adot(1,j) +                           &
     &                 strmix(i,2)*adot(2,j) +                           &
     &                 strmix(i,3)*adot(3,j)
       enddo
       enddo

       gggg = ZERO
       do i = 1,3
       do j = 1,3
         gggg = gggg + adm1vad(i,j)*adm1vad(j,i)
       enddo
       enddo

       tram1v = adm1vad(1,1) + adm1vad(2,2) + adm1vad(3,3)

       do j = 1,3
       do i = 1,3
         frcell(i,j) = strcov(i,j)/(2*vcell*vcell) +                     &
     &      celmas*(vadot(i,1)*adm1vad(1,j) +                            &
     &              vadot(i,2)*adm1vad(2,j) +                            &
     &              vadot(i,3)*adm1vad(3,j)) -                           &
     &       celmas*tram1v*vadot(i,j) +                                  &
     &       celmas*gggg*adot(i,j)/2
       enddo
       enddo

       if(istmd == 1) then
         do j=1,3
         do i=1,3
           adot1(i,j) = adot(i,j)
           vadot1(i,j) = vadot(i,j)
           frcel1(i,j) = frcell(i,j) - beta*celmas*vadot(i,j)
           adot(i,j) = adot1(i,j) + tstep*(vadot1(i,j) +                 &
     &                     tstep*frcel1(i,j)/(2*celmas))
           vadot(i,j) = vadot1(i,j) + tstep*frcel1(i,j)/celmas
         enddo
         enddo
         ekcell = celmas*vcell*vcell*gggg/2
       else
         do j=1,3
         do i=1,3
           adot0(i,j) = adot1(i,j)
           adot1(i,j) = adot(i,j)
           frcel0(i,j) = frcel1(i,j)        
           frcel1(i,j) = frcell(i,j) - beta*celmas*vadot(i,j)       
           vadot1(i,j) = (adot1(i,j) - adot0(i,j)) / tstep +             &
     &         tstep * (2*frcel1(i,j) + frcel0(i,j))/(6*celmas)
           adot(i,j) = 2*adot1(i,j) - adot0(i,j) +                       &
     &           tstep*tstep*frcel1(i,j)/celmas
           vadot(i,j) = vadot1(i,j) + tstep*                             &
     &              (3*frcel1(i,j) - frcel0(i,j))/(2*celmas)
         enddo
         enddo
         ekcell = celmas*vcell*vcell*gggg/2
       endif

!      avoids accumulation of round-off errors

       adot(1,2) = (adot(1,2)+adot(2,1))/2
       adot(2,1) = adot(1,2)
       adot(2,3) = (adot(2,3)+adot(3,2))/2
       adot(3,2) = adot(2,3)
       adot(3,1) = (adot(3,1)+adot(1,3))/2
       adot(1,3) = adot(3,1)

       vadot(1,2) = (vadot(1,2)+vadot(2,1))/2
       vadot(2,1) = vadot(1,2)
       vadot(2,3) = (vadot(2,3)+vadot(3,2))/2
       vadot(3,2) = vadot(2,3)
       vadot(3,1) = (vadot(3,1)+vadot(1,3))/2
       vadot(1,3) = vadot(3,1)

       return
       end subroutine move_vcs_langevin
