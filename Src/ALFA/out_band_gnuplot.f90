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

!>     writes the file for later use with gnuplot.
!>     Usage: gnuplot "filename"

       subroutine out_band_gnuplot(filename,io,neig,nrk,xk,e_of_k,       &
     &        eref,nvert,xcvert,nlines,ljump,nkstep,label,xklab)

!      version 4.53. 19 October 2013. jlm
!      modified (eref) 5 February 2014. jlm
!      Modified, documentation, 4 February 2020. JLM
!      copyright  Jose Luis Martins/INESC-MN

       implicit none
       integer, parameter          :: REAL64 = selected_real_kind(12)

!      input

       character(len=*), intent(in)       ::  filename                   !<  file to be written
       integer, intent(in)                ::  io                         !<  tape number 
       integer, intent(in)                ::  neig                       !<  number of bands
       integer, intent(in)                ::  nrk                        !<  number of k-vectors
       real(REAL64), intent(in)           ::  xk(nrk)                    !<  x coordinate of k-point in plot
       real(REAL64), intent(in)           ::  e_of_k(neig,nrk)           !<  band energies of k-point in plot
       real(REAL64), intent(in)           ::  eref                       !<  reference energy for plot

       integer, intent(in)                ::  nvert                      !<  number of vertical lines in plot
       real(REAL64), intent(in)           ::  xcvert(nvert)              !<  x coordinate of vertical line
       integer, intent(in)                ::  nlines                     !<  number of lines in reciprocal space
       logical, intent(in)                ::  ljump(nlines)              !<  indicates if the new line contains a jump from the preceeding
       integer, intent(in)                ::  nkstep(nlines)             !<  number of steps in line
       character(len=6), intent(in)       ::  label(nvert+nlines)        !<  label of symmetry k-points
       real(REAL64), intent(in)           ::  xklab(nvert+nlines)        !<  x coordinate of label

!      local

       real(REAL64)      ::  ymax, ymin, ymtmp
       character(len=3)  ::  cadd

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

       write(io,'("set terminal wxt enhanced")')
       write(io,*)
       write(io,'("set ylabel ''Energy (eV)'' ",                        &
     &       " font ''Helvetica-Bold , 16'' ")')
       write(io,'("set xrange [ ",f12.6," : ",f12.6," ]")')             &
     &                                  xcvert(1),xcvert(nvert)
       write(io,'("set yrange [ ",f12.6," : ",f12.6," ]")')             &
     &                                  ymin,ymax
       write(io,'("set border lw 3 ")')
       write(io,'("set xtics scale 0")')
       write(io,'("set xtics font ''Helvetica-Bold'' ")')

       do n = 1,nvert+nlines
         if(n == 1) then
           cadd = '   '
         else
           cadd = 'add'
         endif
         if(adjustl(trim(label(n))) == 'Gamma' .or.                     &
     &      adjustl(trim(label(n))) == 'gamma' .or.                     &
     &      adjustl(trim(label(n))) == 'GAMMA') then
           write(io,'("set xtics ",a3," (''{/Symbol G}'' ",             &
     &       f12.6," ) ")') cadd,xklab(n)
         elseif(adjustl(trim(label(n))) == 'Lambda' .or.                &
     &      adjustl(trim(label(n))) == 'lambda' .or.                    &
     &      adjustl(trim(label(n))) == 'LAMBDA') then
           write(io,'("set xtics ",a3," (''{/Symbol L}'' ",             &
     &       f12.6," ) ")') cadd,xklab(n)
         elseif(adjustl(trim(label(n))) == 'Sigma' .or.                 &
     &      adjustl(trim(label(n))) == 'sigma' .or.                     &
     &      adjustl(trim(label(n))) == 'SIGMA') then
           write(io,'("set xtics ",a3," (''{/Symbol S}'' ",             &
     &       f12.6," ) ")') cadd,xklab(n)
         elseif(adjustl(trim(label(n))) == 'Delta' .or.                 &
     &      adjustl(trim(label(n))) == 'delta' .or.                     &
     &      adjustl(trim(label(n))) == 'Delta') then
           write(io,'("set xtics ",a3," (''{/Symbol D}'' ",             &
     &       f12.6," ) ")') cadd,xklab(n)
         else
           l = len(adjustl(trim(label(n))))
           if(l == 1) then
             write(io,'("set xtics ",a3," (''",a1,"'' ",f12.6," ) ")')   &
     &                      cadd,adjustl(trim(label(n))),xklab(n)
           elseif(l == 2) then
             write(io,'("set xtics ",a3," (''",a2,"'' ",f12.6," ) ")')   &
     &                      cadd,adjustl(trim(label(n))),xklab(n)
           elseif(l == 3) then
             write(io,'("set xtics ",a3," (''",a3,"'' ",f12.6," ) ")')   &
     &                      cadd,adjustl(trim(label(n))),xklab(n)
           else
             write(io,'("set xtics ",a3," (''",a6,"'' ",f12.6," ) ")')   &
     &                      cadd,adjustl(trim(label(n))),xklab(n)
           endif
         endif
       enddo

       write(io,'("set ytics font ''Helvetica-Bold'' ")')
       write(io,'("plot ''-'' w lines notitle ; pause -1")')
       
       do j=1,neig
         write(io,*)
         irk = 0
         do n=1,nlines
           if(ljump(n)) then
             write(io,*)
             istart = 0
           else
             istart = 1
           endif
           do i=istart,nkstep(n)
             irk = irk + 1
             write(io,'(2f16.6)') xk(irk),(e_of_k(j,irk)-eref)*EV
           enddo
         enddo
       enddo

       do jrk = 1,nvert
         write(io,*)
         write(io,'(2f16.6)') xcvert(jrk),ymin
         write(io,'(2f16.6)') xcvert(jrk),ymax
       enddo

       close(unit=io)

       return
       end subroutine out_band_gnuplot
