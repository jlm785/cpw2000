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

!>  Plots the quantities ave (and dave if nwidth /= 0)
!>  There are a few assumptions about what is being plotted....
!>  Calls the gnuplot package  www.gnuplot.info

subroutine plot_z1D_gnuplot(ioreplay, io, ave, dave, nw, n3, height,     &
                         filename, title, ylabel, linter)


! Writen June 4, 2014. JLM
! Modified 28 June 2017.  Backslash hack.  JLM
! Modified, documentation, June 11, 2020. JLM
! Modified, name, API, 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations
  integer, intent(in)                ::  io                              !<  "tape" number

  integer, intent(in)                ::  n3                              !<  number of points

  real(REAL64), intent(in)           ::  ave(n3)                         !<  quantity to be plotted
  real(REAL64), intent(in)           ::  dave(n3)                        !<  dave is also plotted if nw =/= 0
  integer, intent(in)                ::  nw                              !<  if =0 only ave is plotted
  real(REAL64), intent(in)           ::  height                          !<  range on horizontal axis

  character(len=*), intent(in)       ::  filename                        !<  filename for plot 
  character(len=*), intent(in)       ::  title                           !<  title of plot
  character(len=*), intent(in)       ::  ylabel                          !<  label of y-axis
  logical, intent(in)                ::  linter                          !<  if .TRUE. asks if the figure should be shown 

! other variables

  character(len=1)  ::  yesno
  character(len=2)  ::  backsl = '\\'                               !  Hack to avoid the special meaning of backslash

! counters

  integer      ::  k

! constants

  real(REAL64), parameter :: BOHR = 0.5291772109_REAL64


! writes the file

  open(unit=io, file=trim(filename), status='UNKNOWN', form='FORMATTED')

  write(io,*) "set terminal wxt enhanced"
  write(io,*) "set encoding iso_8859_1" 
  write(io,*)
  if(nw == 0) then
    write(io,*) "set title  '",title," ' font 'Helvetica-Bold' "
  else
    write(io,*) "set title  '",title," average and ",backsl(1:1)
    write(io,'(i4," z-averages '' font ''Helvetica-Bold'' ")') nw
  endif
  write(io,*) "set ylabel ' ",ylabel," ' font  'Helvetica-Bold' "
  write(io,*) "set xlabel 'z ({/E \305})' font 'Helvetica-Bold' "
  write(io,*) "set xrange [ ",0," : ",height*BOHR," ] " 
  write(io,*) "set border lw 3"
  write(io,*) "set xtics font 'Helvetica-Bold' "
  write(io,*) "set ytics font 'Helvetica-Bold' "
  write(io,*) "plot '-' w lines notitle ; pause -1 "
  write(io,*)
  
  do k=1,n3
    write(io,'(2g14.6)') (k-1)*height*BOHR/real(n3), ave(k)
  enddo
  
  if(nw /= 0) then
    write(io,*)
  do k=1,n3
    write(io,'(2g14.6)') (k-1)*height*BOHR/real(n3), dave(k)
  enddo
  endif

  close(unit=io)

! asks whether a figure should be displayed on screen

  write(6,*)
  if(nw == 0) then
    write(6,*) '  The figure of the average ',title,                &
      ' was written to ',filename
  else
    write(6,*) '  The figure of the ',title,                        &
      ' average and double average was written to ',filename
  endif

  if(linter) then

    write(6,*)
    write(6,*) '  Do you want to see the plot now? (y/n)'
    write(6,*)
  
    read(5,*) yesno
    write(ioreplay,'(2x,a1,"   see interactive plot")') yesno

    if(yesno == 'y' .or. yesno == 'Y') then

      write(6,*)
      write(6,*) '  Hit return to continue '
      write(6,*)

!      call system("gnuplot " // filename // "  2> /dev/null ")
      call execute_command_line("gnuplot " // filename // "  2> /dev/null ")

    endif

  endif

  return

end subroutine plot_z1D_gnuplot
