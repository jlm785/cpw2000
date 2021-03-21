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

!>     writes a rough ASCII plot of the density of states
!>     and the integrated density of states to the default output.
!>     It is commented for gnuplot or xmgrace!

       subroutine dos_print_ascii(emin,deltae,nhist,dhist,chist,lidos,   &
     &                            iskip,nvert)

!      written November 9 1987. jlm 
!      modified November 17, 2013. jlm
!      Modified, documentation, 19 September 2020. JLM
!      copyright  J.L.Martins, INESC-MN.

!      version 4.53 of cpw

       implicit none

       integer, parameter  :: REAL64 = selected_real_kind(12)

!      input:

       real(REAL64), intent(in)        ::  emin                          !<  minimum energy for plot
       real(REAL64), intent(in)        ::  deltae                        !<  energy step in plot
       integer, intent(in)             ::  nhist                         !<  number of points in histograms
       real(REAL64), intent(in)        ::  dhist(nhist)                  !<  density of states
       real(REAL64), intent(in)        ::  chist(nhist)                  !<  integrated density of states
       logical, intent(in)             ::  lidos                         !<  true if integrated density of states is to be computed.
       integer, intent(in)             ::  iskip                         !<  only plots every iskip point
       integer, intent(in)             ::  nvert                         !<  number of points in vertical scale

!      local variables

       integer    ::     jd, ji, na, nb 
       real(REAL64)    ::  ddos, didos, b
       real(REAL64)    ::  cmax,dmax
 
!      counters

       integer     :: i, j

!      constants

       real(REAL64), parameter :: ZERO = 0.0_REAL64
       real(REAL64), parameter :: HARTREE = 27.21138386_REAL64
       character(len=1), parameter   ::  iblnk = ' ', isyd = 'D'
       character(len=1), parameter   ::  isyi = 'I'

!      Vertical range of plot

       dmax = ZERO
       cmax = ZERO
       do i=1,nhist
         if (dhist(i) > dmax) dmax = dhist(i)
         if (chist(i) > cmax) cmax = chist(i)
       enddo

!      find scale for the dos

       ddos = dmax/(nvert*HARTREE)

       call plot_step(ddos,b)

       ddos = b * HARTREE

       didos = cmax/nvert

       call plot_step(didos,b)

       didos = b
  
       write(6,'("#")')
       write(6,'("#")')
       write(6,'("#   Density of States  (DOS)")')
       write(6,'("#")')
       write(6,'("#   Energy scale:   E_min=",f10.3," eV,     step=",    &
     &       f10.3," eV")') emin*HARTREE,deltae*iskip*HARTREE
       write(6,'("#   DOS scale:  ",f10.3,"  el/cell/eV")')              &
     &              ddos/HARTREE
       if(lidos) then
         write(6,'("#   Integrated DOS scale: ",f10.3,"  el/cell")')     &
     &              didos
       endif
       write(6,'("#")')

       if(lidos) then

         do i=1,nhist,iskip

           jd = int(dhist(i)/ddos)
           ji = int(chist(i)/didos)

           if(jd == ji) then

             if(jd == 0) then
               write(6,'(120a1)') '#',' ',' ',' ','H'
             else if(jd == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H', isyd
             else
               na = jd-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyd
             endif

           else if(ji == 0) then

             if(jd == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H', isyd
             else
               na = jd-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyd
             endif

           else if(jd == 0) then
             if(ji == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H', isyi
             else
               na = ji-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                    (iblnk,j=1,na),isyi
             endif

           else if(ji == jd+1) then

             if(jd == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H', isyd,isyi
             else
               na = jd-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyd,isyi
             endif

           else if(jd == ji+1) then

             if(ji == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H', isyi,isyd
             else
               na = ji-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyi,isyd
             endif

           else if(ji > jd+1) then

             nb = ji - jd -1
             if(jd == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  isyd,(iblnk,j=1,nb),isyi
             else
               na = jd-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyd,(iblnk,j=1,nb),isyi
             endif

           else if(jd > ji+1) then

             nb = jd - ji -1
             if(ji == 1) then
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  isyi,(iblnk,j=1,nb),isyd
             else
               na = ji-1
               write(6,'(120a1)') '#',' ',' ',' ','H',                   &
     &                  (iblnk,j=1,na),isyi,(iblnk,j=1,nb),isyd
             endif

           endif

         enddo


       else


         do i=1,nhist,iskip

           jd = int(dhist(i)/ddos)

           if(jd == 1) then
             write(6,'(120a1)') '#',' ',' ',' ','H', isyd
           else
             na = jd-1
             write(6,'(120a1)') '#',' ',' ',' ','H',                     &
     &                  (iblnk,j=1,na),isyd
           endif

         enddo

       endif

       return
       end subroutine dos_print_ascii
