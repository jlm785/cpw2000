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

       subroutine out_band_fold_xmgrace(filename,io,                     &
     &        title,subtitle,nstyle,                                     &
     &        pkn,neig,nrk,xk,e_of_k,eref,nocc,                          &
     &        nvert,xcvert,nlines,ljump,nkstep,label,xklab)

!      writes the file for later use with xmgrace
!      xmgrace "filename"

!      version 4.53. 19 October 2013. jlm
!      modified (eref) 5 February 2014. jlm
!      modified (unfold) July 2014. CLR
!      modified title, style, 6 August 2014. JLM
!      copyright  Jose Luis Martins, Carlos Loia Reis/INESC-MN

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len=*), intent(in)       ::  filename                   !  file to be written
       integer, intent(in)                ::  io                         !  tape number 

       integer, intent(in)                ::  neig                       !  number of bands
       integer, intent(in)                ::  nrk                        !  number of k-vectors

       character(len=50), intent(in)      ::  title                      !  title for plots
       character(len=140), intent(in)     ::  subtitle                   !  subtitle for plots
       integer, intent(in)                ::  nstyle                     !  choice of plot style

       real(REAL64), intent(in)           ::  pkn(nrk,neig) 

       real(REAL64), intent(in)           ::  xk(nrk)                    !  x coordinate of k-point in plot
       real(REAL64), intent(in)           ::  e_of_k(neig,nrk)           !  band energies of k-point in plot

       real(REAL64), intent(in)           ::  eref                       !  reference energy for plot
       integer, intent(in)                ::  nocc                       !  number of lines in reciprocal space

       integer, intent(in)                ::  nvert                      !  number of vertical lines in plot
       real(REAL64), intent(in)           ::  xcvert(nvert)              !  x coordinate of vertical line
       integer, intent(in)                ::  nlines                     !  number of lines in reciprocal space
       logical, intent(in)                ::  ljump(nlines)              !  indicates if the new line contains a jump from the preceeding
       integer, intent(in)                ::  nkstep(nlines)             !  number of steps in line
       character(len=6), intent(in)       ::  label(nvert+nlines)        !  label of symmetry k-points
       real(REAL64), intent(in)           ::  xklab(nvert+nlines)        !  x coordinate of label

!      local

       real(REAL64)      ::  ymax, ymin, ymtmp

!      counters

       integer      ::   irk, jrk, n, i, j, l, istart

!      constants

       real(REAL64), parameter ::  EV = 27.21138505_REAL64
       real(REAL64), parameter  :: UM = 1.0_REAL64

!      finds the energy range for the bands

       ymin = e_of_k(1,1) - eref
       do irk=1,nrk
         ymtmp = e_of_k(neig,irk) - eref
         do n=1,neig
           if(ymin > e_of_k(n,irk)-eref) ymin = e_of_k(n,irk)-eref
           if(ymtmp < e_of_k(n,irk)-eref) ymtmp = e_of_k(n,irk)-eref
         enddo
         if(irk == 1) then
           ymax = ymtmp
         else
           if(ymax > ymtmp) ymax = ymtmp
         endif
       enddo
       ymin = real(nint(ymin*27.212)) - UM
       ymax = real(nint(ymax*27.212)) + UM

       open(unit=io,file=filename,form='formatted')

       write(io,'("@    autoscale onread none ")')
       write(io,'("@    znorm 3 ")')

       write(io,'("@    world  ",f18.8,",",f18.8,",",f18.8,",",f18.8)') &
     &  xcvert(1),ymin,xcvert(nvert),ymax
       write(io,*)
       
       if(nstyle == 1) then

         write(io,'("@    frame linewidth 3.0 ")')
         write(io,*)

         write(io,'("@    title font 6")')
         write(io,'("@    title size 1.5")')
         write(io,'("@    title color 1")')
         write(io,*)

         write(io,'("@    subtitle font 6")')
         write(io,'("@    subtitle size 1.0")')
         write(io,*)

         write(io,*)
         write(io,'("@    xaxis  ticklabel font 6 ")')

       endif

       l = len(adjustl(trim(title)))
       if(l > 0) then
         write(io,*)
         write(io,"('@    title ""',a,'""')") adjustl(trim(title))
         write(io,*)
       endif

       l = len(adjustl(trim(subtitle)))
       if(l > 0) then
         write(io,*)
         write(io,"('@    subtitle ""Unfolded ',a,'""')")                &
     &            adjustl(trim(subtitle))
         write(io,*)
       endif
      
       write(io,'("@    xaxis  tick off ")')
       write(io,'("@    xaxis  tick spec type both ")')
       write(io,'("@    xaxis  tick spec ",i6)') nvert+nlines

       do n = 1,nvert+nlines
         write(io,'("@    xaxis  tick major ",i5," , ",f18.8)')         &
     &                      n-1,xklab(n)
         if(adjustl(trim(label(n))) == 'Gamma' .or.                     &
     &      adjustl(trim(label(n))) == 'gamma' .or.                     &
     &      adjustl(trim(label(n))) == 'GAMMA') then
             write(io,"('@    xaxis  ticklabel ',i5,                    &
     &              ' , ""\f{Symbol} G\f{}"" ')")  n-1
         elseif(adjustl(trim(label(n))) == 'Lambda' .or.                &
     &      adjustl(trim(label(n))) == 'lambda' .or.                    &
     &      adjustl(trim(label(n))) == 'LAMBDA') then
             write(io,"('@    xaxis  ticklabel ',i5,                    &
     &              ' , ""\f{Symbol} L\f{}"" ')")  n-1
         elseif(adjustl(trim(label(n))) == 'Sigma' .or.                 &
     &      adjustl(trim(label(n))) == 'sigma' .or.                     &
     &      adjustl(trim(label(n))) == 'SIGMA') then
             write(io,"('@    xaxis  ticklabel ',i5,                    &
     &              ' , ""\f{Symbol} S\f{}"" ')")  n-1
         elseif(adjustl(trim(label(n))) == 'Delta' .or.                 &
     &      adjustl(trim(label(n))) == 'delta' .or.                     &
     &      adjustl(trim(label(n))) == 'Delta') then
             write(io,"('@    xaxis  ticklabel ',i5,                    &
     &              ' , ""\f{Symbol} D\f{}"" ')")  n-1
         elseif(adjustl(trim(label(n))) == '-' .or.                     &
     &      adjustl(trim(label(n))) == '--' .or.                        &
     &      adjustl(trim(label(n))) == '---') then
             write(io,"('@    xaxis  ticklabel ',i5,                    &
     &              ' , ""  "" ')")  n-1
         else
           l = len(adjustl(trim(label(n))))
           if(l == 1) then
             write(io,"('@    xaxis  ticklabel ',i5,' , ""',a1,'"" ')") &
     &                      n-1,adjustl(trim(label(n)))
           elseif(l == 2) then
             write(io,"('@    xaxis  ticklabel ',i5,' , ""',a2,'"" ')") &
     &                      n-1,adjustl(trim(label(n)))
           elseif(l == 3) then
             write(io,"('@    xaxis  ticklabel ',i5,' , ""',a3,'"" ')") &
     &                      n-1,adjustl(trim(label(n)))
           else
             write(io,"('@    xaxis  ticklabel ',i5,' , ""',a6,'"" ')") &
     &                      n-1,adjustl(trim(label(n)))
           endif
         endif
       enddo


       if(nstyle == 1) then

         write(io,*)
         write(io,"('@    yaxis  label ""Energy (eV)"" ')")
         write(io,'("@    yaxis  label char size 1.5 ")')
         write(io,'("@    yaxis  label font 6 ")')
         write(io,'("@    yaxis  tick major 5 ")')
         write(io,'("@    yaxis  tick minor ticks 0 ")')
         write(io,'("@    yaxis  tick major linewidth 2.0 ")')
         write(io,'("@    yaxis  ticklabel font 6 ")')

       else

         write(io,*)
         write(io,"('@    yaxis  label ""Energy [eV]"" ')")
         write(io,'("@    yaxis  tick major 5 ")')
         write(io,'("@    yaxis  tick minor ticks 0 ")')

       endif
       

!      // put colors and stuff here
      do j=0,neig-1
         write(io,'("@    s",i3.3," symbol 1")') j
         write(io,'("@    s",i3.3," symbol size 0.330000")') j
         write(io,'("@    s",i3.3," symbol pattern 1")') j
         write(io,'("@    s",i3.3," symbol fill pattern 1")') j
         write(io,'("@    s",i3.3," line type 0")') j
      enddo
            
      do j=0,nocc-1      
         write(io,'("@    s",i3.3," symbol color 4")') j
         write(io,'("@    s",i3.3," symbol fill color 4")') j
      enddo

      do j=nocc,neig-1
         write(io,'("@    s",i3.3," symbol color 2")') j
         write(io,'("@    s",i3.3," symbol fill color 2")') j
      enddo
      
       do jrk = 1,nvert
         write(io,'("@    s",i3.3," line color 1" )') neig-1+jrk
         write(io,'("@    s",i3.3," line type 1" )') neig-1+jrk
         if(nstyle == 1) then
           write(io,'("@    s",i3.3," line linestyle 1" )') neig-1+jrk
           write(io,'("@    s",i3.3," line linewidth 1.5" )') neig-1+jrk
         else
           write(io,'("@    s",i3.3," line linestyle 3" )') neig-1+jrk
           write(io,'("@    s",i3.3," line linewidth 1.0" )') neig-1+jrk
         endif
       enddo

         write(io,'("@    s",i3.3," line color 1" )') neig+nvert
         write(io,'("@    s",i3.3," line type 1" )') neig+nvert
         write(io,'("@    s",i3.3," line linestyle 1" )') neig+nvert
         write(io,'("@    s",i3.3," line linewidth 2.0" )') neig+nvert

       
       do j=1,neig
         write(io,*)
         irk = 0
           write(io,'("@target g0.s",i3.3)') j-1
           write(io,'("@type xysize ")')
         do n=1,nlines
           if(ljump(n)) then
!             write(io,*)
             istart = 0
           else
             istart = 1
           endif
           do i=istart,nkstep(n)
             irk = irk + 1
             write(io,'(3f18.8)') xk(irk),(e_of_k(j,irk)-eref)*EV,pkn(irk,j)
           enddo
         enddo
       enddo

       do jrk = 1,nvert
         write(io,'("@target g0.s",i3.3)') neig-1+jrk
         write(io,'("@type xy ")')
         write(io,'(2f18.8)') xcvert(jrk),ymin
         write(io,'(2f18.8)') xcvert(jrk),ymax
       enddo

         write(io,'("@target g0.s",i3.3)') neig+nvert
         write(io,'("@type xy ")')
         write(io,'(2f18.8)') xk(1),0.0D0
         write(io,'(2f18.8)') xk(nrk),0.0D0



       close(unit=io)

       return
       end subroutine out_band_fold_xmgrace
