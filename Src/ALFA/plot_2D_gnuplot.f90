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

!>  prepares a file to be used  with gnuplot
!>  with a two dimensional contour plot

subroutine plot_2D_gnuplot(ioreplay,filename,io,                  &
       ro,nx,ny,xscale,yscale)

! modified, f90, subroutine, 27 May 2014. JLM
! Modified, documentation, close correct unit, 11 June 2020. JLM
! Name, 4 february 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.99


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  character(len=*), intent(in)       ::  filename                        !<  file to be written
  integer, intent(in)                ::  io                              !<  tape number 

  integer, intent(in)                ::  nx, ny                          !<  Dimensions of grid in plane

  real(REAL64), intent(in)           ::  ro(nx,ny)                       !<  charge density interpolated on the planar grid
  real(REAL64), intent(in)           ::  xscale,yscale                   !<  aspect ratio of plot

! other variables

  character(len=1)  ::  yesno
  real(REAL64)      ::  delta
  real(REAL64)      ::  zmin, zmax, zincr
  real(REAL64)      ::  xscl, yscl
  integer           ::  ncont, ixe

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer i, j




  zmin = ro(1,1)
  zmax = ro(1,1)
  do j=1,ny
  do i=1,nx
    if(ro(i,j) < zmin) zmin = ro(i,j)
    if(ro(i,j) > zmax) zmax = ro(i,j)
  enddo
  enddo
  
! approximate number of contours

  ncont = 15

  delta = (zmax-zmin)/ncont
  
  call plot_step(delta,zincr)

! find nice limits

  ixe = int(zmin/zincr)
  zmin = zincr*(ixe-1)
  ixe = int(zmax/zincr)
  zmax = zincr*(ixe+1)
  ncont = nint((zmax-zmin)/zincr)

  write(6,*)
  write(6,'("  Values for plot:  zmin,zmax,zincr = ",3g10.3)')      &
                        zmin,zmax,zincr
  write(6,*)

! in principle one of xscale or yscale is 1

  xscl = min(UM,xscale)
  if(xscl < ZERO) xscl = UM
  yscl = min(UM,yscale)
  if(yscl < ZERO) yscl = UM

  write(6,*)
  write(6,*) '  Do you want to see a nice plot? (y/n)'
  write(6,*)
  
  read(5,*) yesno
  write(ioreplay,'(2x,a1,"   nice plot")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    open(unit=io,file=filename,form='FORMATTED')

!   commands

    write(io,'("set terminal wxt size ",g12.3,",",g12.3,            &
        " enhanced")') 800*xscl,800*yscl
    write(io,*)
    write(io,*) 'set noparam'
    write(io,*) 'set contour base'
    write(io,*) 'unset xtics'
    write(io,*) 'unset ytics'
    write(io,*) 'unset ztics'
    write(io,*) '#  comment next line if you want contour values'
    write(io,*) 'unset key'
    write(io,*) '#  comment next two lines if you want 3D view'
    write(io,*) 'set nosurface'
    write(io,*) 'set view 0,0'
    write(io,'("set cntrparam levels incremental",e10.3," ,  ",     &
         e10.3," , ",e10.3)') zmin,zincr,zmax
    write(io,*)  "splot '-' w lines notitle; pause -1 "

    write(io,'("#  mesh size is:",2i5)') nx,ny
    write(io,*)

    do j=1,ny
      do i=1,nx
        write(io,'(2x,f16.8)') ro(i,j)
      enddo
      write(io,'("#new line  ")')
      write(io,*)
    enddo

    close(unit=io)
  
    write(6,*)
    write(6,*) '  Hit return to continue '
    write(6,*)

!    syscom = 'gnuplot ' // trim(filename) // ' 2> /dev/null  '

!    call system(trim(syscom))
    call execute_command_line('gnuplot ' // trim(filename) // ' 2> /dev/null  ')
    
    write(6,*)
    write(6,*) '  You can edit the file ',trim(filename)
    write(6,*) '  if you want to change the plot style.'
    write(6,*)
    
  endif


  write(6,*)
  write(6,*) '  Do you want a .png file with the plot? (y/n)'
  write(6,*)
  
  read(5,*) yesno
  write(ioreplay,'(2x,a1,"   png file")') yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    open(unit=io,file=('png'//filename),form='FORMATTED')

!   commands

    write(io,*) '#  edit this file to change style'
    write(io,'("set terminal png size ",g12.3,",",g12.3)')          &
         800*xscl,800*yscl
    write(io,*) 'set output "density2D.png" '
    write(io,*)
    write(io,*) 'set noparam'
    write(io,*) 'set contour base'
    write(io,*) 'unset xtics'
    write(io,*) 'unset ytics'
    write(io,*) 'unset ztics'
    write(io,*) '#  comment next line if you want contour values'
    write(io,*) 'unset key'
    write(io,*) '#  comment next two lines if you want 3D view'
    write(io,*) 'set nosurface'
    write(io,*) 'set view 0,0'
    write(io,'("set cntrparam levels incremental",e10.3," ,  ",     &
         e10.3," , ",e10.3)') zmin,zincr,zmax
    write(io,*)  "splot '-' w lines notitle "

    write(io,'("#  mesh size is:",2i5)') nx,ny
    write(io,*)

    do j=1,ny
      do i=1,nx
        write(io,'(2x,f10.5)') ro(i,j)
      enddo
      write(io,'("#new line  ")')
      write(io,*)
    enddo

    close(unit=io)
 

    call execute_command_line('gnuplot ' // trim('png'//filename) // ' 2> /dev/null  ')

    write(6,*)
    write(6,*) '  The plot is in file density2D.png'
    write(6,*)
    write(6,*) '  You can edit the file ','png'//trim(filename)
    write(6,*) '  if you want to change the plot style.'
    write(6,*)
    write(6,*) '  Those files can be seen with gnuplot'
    write(6,*)
    
  endif


  return

end subroutine plot_2D_gnuplot
