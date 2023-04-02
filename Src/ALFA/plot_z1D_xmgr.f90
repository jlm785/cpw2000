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
!>  Writes a file to be plotted by xmgrace
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      5.06
!>  \date         21 June 1014, 2 March 2023.
!>  \copyright    GNU Public License v2

subroutine plot_z1D_xmgr(io, ave, dave, nw, n3, height,                  &
                         filename, title, ylabel)


! Modified, Documentation, name, API, 4 February 2021. JLM
! Modified, first line for KDE recognition. 20 January 2022. JLM
! Modified, several double averages, 2 March 2023. JLM

  implicit none


  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  io                              !< "tape" number

  integer, intent(in)                ::  n3                              !< number of points
  integer, intent(in)                ::  nw                              !< if = 0 only ave is plotted

  real(REAL64), intent(in)           ::  ave(n3)                         !< quantity to be plotted
  real(REAL64), intent(in)           ::  dave(n3,max(1,nw))              !< dave is also plotted if nw > 0
  real(REAL64), intent(in)           ::  height                          !< range on horizontal axis

  character(len=*), intent(in)       ::  filename                        !< filename for plot
  character(len=*), intent(in)       ::  title                           !< title of plot
  character(len=*), intent(in)       ::  ylabel                          !< label of y-axis

! other variables

  real(REAL64)      ::  ymin,ymax
  real(REAL64)      ::  xstep,ystep

! counters

  integer      ::  k, n

! constants

  real(REAL64), parameter :: BOHR = 0.5291772109_REAL64


  ymin = ave(1)
  ymax = ave(1)
  do k = 1,n3
    if(ymin > ave(k)) ymin = ave(k)
    if(ymax < ave(k)) ymax = ave(k)
  enddo
  if(nw /= 0) then
    do k = 1,n3
      if(ymin > ave(k)) ymin = ave(k)
      if(ymax < ave(k)) ymax = ave(k)
    enddo
  endif

  call plot_step(height*BOHR/5,xstep)
  call plot_step((ymax-ymin)/5,ystep)

  ymax = ystep*(int(ymax/ystep)+1)
  ymin = ystep*(int(ymin/ystep)-1)

! writes the file

  open(unit=io, file=trim(filename), status='UNKNOWN', form='FORMATTED')

  write(io,'("# Grace project file ")')
  write(io,*)

  write(io,'("@    autoscale onread none ")')
  write(io,'("@    world  ",f12.6,",",f12.6,",",f12.6,",",f12.6)')  &
   0,ymin,height*BOHR,ymax
  write(io,*)
  write(io,'("@    frame linewidth 3.0 ")')
  write(io,*)
  write(io,"('@    xaxis  label ""z (\cE\C)"" ')")
  write(io,'("@    xaxis  label char size 1.5 ")')
  write(io,'("@    xaxis  label font 6 ")')
  write(io,'("@    xaxis  tick major ",f12.4)') xstep
  write(io,'("@    xaxis  tick minor ticks 0 ")')
  write(io,'("@    xaxis  tick major linewidth 2.0 ")')
  write(io,'("@    xaxis  ticklabel font 6 ")')
  write(io,*)
  write(io,*) '@    yaxis  label " ',ylabel,' " '
  write(io,'("@    yaxis  label char size 1.5 ")')
  write(io,'("@    yaxis  label font 6 ")')
  write(io,'("@    yaxis  tick major ",f12.4)') ystep
  write(io,'("@    yaxis  tick minor ticks 0 ")')
  write(io,'("@    yaxis  tick major linewidth 2.0 ")')
  write(io,'("@    yaxis  ticklabel font 6 ")')
  write(io,*)
  write(io,*) '@   title  " ',title,' " '
  write(io,'("@    title font 6 ")')
  if(nw == 0) then
  else
    write(io,"('@   subtitle  ""Planar average and',i3,' z-averages"" ')") nw
    write(io,'("@   subtitle font 6 ")')
  endif

  write(io,*)

  do k=1,n3
    write(io,'(2g14.6)') (k-1)*height*BOHR/real(n3), ave(k)
  enddo

  if(nw > 0) then
    do n = 1,nw
      write(io,*)
      do k = 1,n3
        write(io,'(2g14.6)') (k-1)*height*BOHR/real(n3), dave(k,n)
      enddo
    enddo
  endif

  close(unit=io)

  write(6,*)
  if(nw == 0) then
    write(6,*) '  The figure of the average ',title,                &
      ' was written to ',filename
  else
    write(6,*) '  The figure of the ',title,                        &
      ' average and double average was written to ',filename
  endif

  return

end subroutine plot_z1D_xmgr
