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

!>     Prints the information about the Fermi level, occupations
!>
!>  \author       Jose Luis Martins
!>  \version      5.0.3
!>  \date         September 15, 2015, 29 November 2021.
!>  \copyright    GNU Public License v2

       subroutine print_fermi_occup(ipr,el,teleck,                       &
     & nrk,wgk,nband,                                                    &
     & frac,efermi,eband,elects,bandwid,penngap,                         &
     & mxdnrk,mxdbnd)

!      written September 15, 2015. JLM
!      Fractional level, October 15, 2017. JLM
!      based on older fermi_level subroutine
!      Modified, documentation, August 2019. JLM
!      Modified, increase threshold for printing details. 29 November 2021. JLM
!      copyright inesc-mn/Jose Luis Martins

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       integer, intent(in)                ::  mxdnrk                     !<  size of k-points
       integer, intent(in)                ::  mxdbnd                     !<  size of number of bands

       integer, intent(in)                ::  ipr                        !<  controls printing (0,1,2)
       real(REAL64), intent(in)           ::  el(mxdbnd*mxdnrk)          !<  eigenvalues in Hartree for all the k-points,
       real(REAL64), intent(in)           ::  teleck                     !<  temperature (Kelvin) of electron system.
       integer, intent(in)                ::  nrk                        !<  number of k-points for integration in the irreducible wedge of the brillouin zone
       real(REAL64), intent(in)           ::  wgk(mxdnrk)                !<  weight in the integration of k-point
       integer, intent(in)                ::  nband(mxdnrk)              !<  number of bands for each k-points
       real(REAL64), intent(in)           ::  frac(mxdbnd*mxdnrk)        !<  fractional ocupation of level j
       real(REAL64), intent(in)           ::  efermi                     !<  Fermi energy (or highest occupied state)
       real(REAL64), intent(in)           ::  eband                      !<  band energy (Hartree)
       real(REAL64), intent(in)           ::  elects                     !<  electronic temperature*entropy (hartree)
       real(REAL64), intent(in)           ::  bandwid                    !<  occupied band width estimate (Hartree)
       real(REAL64), intent(in)           ::  penngap                    !<  Penn gap estimate (Hartree)

!      local allocatable arrays

       integer, allocatable               ::  jrk(:)
       integer, allocatable               ::  ind(:)

!      local variables

       integer           ::  iel, jmin, jmax
       real(REAL64)      ::  tempau

!      counters

       integer   ::    i, j, k

!      constants

       real(REAL64), parameter ::  UM = 1.0_REAL64
       real(REAL64), parameter ::  EPS = 0.000001_REAL64
       real(REAL64), parameter ::  SMALL = EPS*EPS
       real(REAL64), parameter ::  EV = 27.2116_REAL64
       real(REAL64), parameter ::  TAUTOK = 11604.9_REAL64 * EV

       allocate(jrk(mxdnrk*mxdbnd))
       allocate(ind(mxdnrk*mxdbnd))

       iel = 0
       do i=1,nrk
         do j=1,nband(i)
           iel = iel + 1
           jrk(iel) = i
         enddo
       enddo

!      sort the energy levels

       call sort(iel,el,ind)

       tempau = teleck / TAUTOK

       if(ipr > 0) write(6,'(/,"  the fermi level is at ",f10.4,         &
     &       " [eV] ",/)') efermi*EV

       if(ipr > 2) then

         write(6,*)

         if( tempau < EPS*abs(el(ind(iel))-el(ind(1))) ) then

!          zero temperature

           do i = 1,iel
             if(abs(frac(ind(i))*(UM-frac(ind(i)))) > SMALL) then
               write(6,'(/,"  NOTE: Fractional occupancy at k point ",   &
     &              i3," frac = ",f9.6," with energy ",f10.4," [eV] ")') &
     &              jrk(ind(i)),frac(ind(i)),el(ind(i))*EV
             endif
           enddo

         else

           iel = 0
           do i=1,nrk
             write(6,'("    band ",i5,"    weigth = ",f10.5)') i,wgk(i)
             jmin = nband(i)/10
             if(jmin > 0) then
               do j=1,jmin
                 write(6,'(10f10.4)') (el(iel+k)*EV,k=1,10)
                 write(6,'(10f10.4)') (frac(iel+k),k=1,10)
                 iel = iel + 10
               enddo
             endif
             jmax = nband(i) - 10*jmin
             if(jmax > 0) then
               write(6,'(10f10.4)') (el(iel+k)*EV,k=1,jmax)
               write(6,'(10f10.4)') (frac(iel+k),k=1,jmax)
               iel = iel + jmax
             endif
           enddo

         endif

         write(6,*)

       endif

       if(ipr > 2) then
         write(6,*)
         write(6,'(5x,f14.6,"  band energy (Ha)")') eband
         write(6,'(5x,f14.6,"  TS for electrons (Ha)")') elects
         write(6,'(5x,f14.6,"  ESTIMATE of band width (eV)")')           &
     &             bandwid*EV
         write(6,'(5x,f14.6,"  ESTIMATE of Penn gap (eV)")')             &
     &             penngap*EV
         write(6,*)
       endif

       return
       end subroutine print_fermi_occup
