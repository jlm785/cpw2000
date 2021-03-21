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

!>  quick and dirty ascii contour plot

subroutine plot_alf_plt(ro,nx,ny)

! Modified, f90, 27 May 2014. JLM
! Modified, documentation name, 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99


  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nx, ny                          !<  Dimensions of grid in plane

  real(REAL64), intent(in)           ::  ro(nx,ny)                       !<  function interpolated on the planar grid

! other variables

  real(REAL64)          ::  xmin, xmax
  integer               ::  itype
  real(REAL64)          ::  step
  integer               ::  nxx, ii, nl

  character(len=1)      ::  crow(80)

! counters

  integer      ::  i, j

! parameters

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64
  character(len=1), parameter, dimension(7)  :: posalf = (/' ','-','+','x','h','#','@'/)
  character(len=1), parameter, dimension(7)  :: calf =   (/'#','x','+',' ','-','=','@'/)


  xmin = ro(1,1)
  xmax = ro(1,1)
  do j = 1,ny
  do i = 1,nx
    if(ro(i,j) < xmin) xmin = ro(i,j)
    if(ro(i,j) > xmax) xmax = ro(i,j)
  enddo
  enddo

  if(xmin < ZERO .and. xmax < ZERO) then
    itype = -1
    step = (xmin-xmax)/5.9999
  elseif(xmin > ZERO .and. xmax > ZERO) then
    itype =  1
    step = (xmax-xmin)/5.9999
  else
    step = max(abs(xmin),abs(xmax))
    xmin = -step
    xmax = step
    step = (xmax-xmin)/5.9999
    itype = 0
  endif

  write(6,*)
  write(6,'("   xmin,xmax,step =",3e10.3)') xmin,xmax,step
  write(6,*)

  do j=1,ny
    do i=1,80
      crow(i) = ' '
    enddo
    nxx = (nx-1)/80 + 1
    ii = 0
    do i=nxx,nx,nxx
      ii = ii + 1
      if(itype == 1) then
        nl = int((ro(i,j)-xmin)/step)+1
        crow(ii) = posalf(nl)
      elseif(itype == -1) then
        nl = int((ro(i,j)-xmax)/step)+1
        crow(ii) = posalf(nl)
      else
        nl = int((ro(i,j)-xmin)/step)+1
        crow(ii) = calf(nl)
      endif
    enddo
    write(6,'(80a1)') (crow(i),i=1,80)
  enddo

  return

  end subroutine plot_alf_plt
