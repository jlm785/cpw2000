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

!>  Determines the plane for contour plot
!>  c0(i) is the corner and two vectors
!>  vx(i) and vy(i) give two sides. these are all given in units
!>  of the lattice basis vectors.

subroutine plot_get_plane(ioreplay, adot,                                &
             nx, ny, c0, dx,dy, xscale,yscale)

! Writen 5 June 2014 fom older code of February 11,2008. JLM
! Documentation, name, 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99


  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  integer, intent(in)                ::  nx, ny                          !<  Dimensions of grid in plane

! output

  real(REAL64), intent(out)         ::  c0(3)                            !<  corner (origin) of plane
  real(REAL64), intent(out)         ::  dx(3),dy(3)                      !<  step vectors that define the plane (lattice coordinates)
  real(REAL64), intent(out)         ::  xscale,yscale                    !<  aspect ratio of plot (either xscale or yscale = 1)

! other variables

  real(REAL64)      ::  vx(3),vy(3)
  real(REAL64)      ::  xlat(3)
  real(REAL64)      ::  xymax
  real(REAL64)      ::  xy,xl,yl
  character(len=1)  ::  yesno
  integer           ::  icoor               !  coordinate system

  real(REAL64)      ::  avec(3,3)           !  primitive lattice vectors
  real(REAL64)      ::  bvec(3,3)           !  reciprocal primitive lattice vectors
  real(REAL64)      ::  aconv(3,3)          !  conventional lattice vectors of the Niggli cell
  real(REAL64)      ::  avecnig(3,3)        !  primitive lattice vectors of the Niggli cell

! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: ZERO = 0.0_REAL64

! counters

  integer i, j


! gets the lattice vectors

  call adot_to_avec_aconv(adot,avec,bvec,aconv,avecnig)


  write(6,*)
  write(6,*) '   Which coordinate system do you want to use?'
  write(6,*) '   1) Primitive lattice coordinates.'
  write(6,*) '   2) Cartesian coordinates.'
  write(6,*) '   3) Conventional lattice coordinates.'
  write(6,*)
  write(6,*) '   Enter your choice (1-3).'
  write(6,*)
  
  read(5,*) icoor
  write(ioreplay,'(2x,i8,"   coordinate choice")') icoor
  
  if(icoor <1 .or. icoor > 3) then
    write(6,*) '   Wrong choice, using primitive lattice'
    icoor = 1
  endif

! prints information about the coordinate system

  if(icoor == 1) then

    write(6,*)
    write(6,*) '   using PRIMITIVE LATTICE COORDINATES '
    write(6,*)
    write(6,*) '   The primitive lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",       &
        f8.4,")")') avec(1,1),avec(1,2),avec(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",      &
        f8.4,")")') avec(2,1),avec(2,2),avec(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",      &
        f8.4,")")') avec(3,1),avec(3,2),avec(3,3)
     write(6,*)

  elseif(icoor == 2) then

    write(6,*)
    write(6,*) '   using CARTESIAN COORDINATES (atomic units)'
    write(6,*)

  elseif(icoor == 3) then

    write(6,*)
    write(6,*) '   using CONVENTIONAL LATTICE COORDINATES '
    write(6,*)
    write(6,*) '   The conventional lattice vectors are: '
    write(6,*)
    write(6,'("  a1 = (",f8.4,")   a2 = (",f8.4,")   a3 = (",       &
        f8.4,")")') aconv(1,1),aconv(1,2),aconv(1,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",      &
        f8.4,")")') aconv(2,1),aconv(2,2),aconv(2,3)
     write(6,'("       (",f8.4,")        (",f8.4,")        (",      &
        f8.4,")")') aconv(3,1),aconv(3,2),aconv(3,3)
     write(6,*)

  endif

! enters data

  write(6,*)
  write(6,*) '   enter corner of plane,  x1,x2,x3 '
  write(6,*)

  read(5,*) c0(1),c0(2),c0(3)
  write(ioreplay,'(3(2x,f20.10),"   corners")') c0(1),c0(2),c0(3)

  write(6,*)
  write(6,*) '   enter x direction of plane, v1,v2,v3 '
  write(6,*)

  read(5,*) vx(1),vx(2),vx(3)
  write(ioreplay,'(3(2x,f20.10),"   1st direction")')               &
                vx(1),vx(2),vx(3)

  write(6,*)
  write(6,*) '   enter y direction of plane,  u1,u2,u3 '
  write(6,*)

  read(5,*) vy(1),vy(2),vy(3)
  write(ioreplay,'(3(2x,f20.10),"   2nd direction")')               &
                 vy(1),vy(2),vy(3)

! transforms back to lattice coordinates

  if(icoor == 2) then

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + c0(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      c0(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + vx(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      vx(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + vy(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      vy(i) = xlat(i)
    enddo
    
  elseif(icoor == 3) then

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + aconv(i,j)*c0(j)
      enddo
    enddo
    do i=1,3
      c0(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + aconv(i,j)*vx(j)
      enddo
    enddo
    do i=1,3
      vx(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + aconv(i,j)*vy(j)
      enddo
    enddo
    do i=1,3
      vy(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + c0(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      c0(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + vx(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      vx(i) = xlat(i)
    enddo

    do i=1,3
      xlat(i) = ZERO
      do j=1,3
        xlat(i) = xlat(i) + vy(j)*bvec(j,i)
      enddo
      xlat(i) = xlat(i) / (2*PI)
    enddo
    do i=1,3
      vy(i) = xlat(i)
    enddo
    
  endif

! check that vx and vy are orthogonal, and find the ratio
! between their lengths.

  xy = ZERO
  xl = ZERO
  yl = ZERO
  do i=1,3
    do j=1,3
      xy = xy + vx(i)*adot(i,j)*vy(j)
      xl = xl + vx(i)*adot(i,j)*vx(j)
      yl = yl + vy(i)*adot(i,j)*vy(j)
    enddo
  enddo
  xl = sqrt(xl)
  yl = sqrt(yl)

  if (abs(xy)/(xl*yl) > 1.0d-6) then
    write(6,'("  WARNING  vectors in cplot are not ",               &
          "orthogonal  cos(fi) = ",f10.3)') xy/(xl*yl)
    write(6,*)
    write(6,*) 'do you want to orthogonalize second vector? (y/n)'
    write(6,*)
    
    read(5,*) yesno
    write(ioreplay,'(2x,a1,"   orthogonal axis")') yesno
    
    if(yesno == 'y' .or. yesno == 'y') then
      vy(1) = vy(1) - (xy/(xl*xl))*vx(1)
      vy(2) = vy(2) - (xy/(xl*xl))*vx(2)
      vy(3) = vy(3) - (xy/(xl*xl))*vx(3)
      yl = 0.0
      do i=1,3
        do j=1,3
          yl = yl + vy(i)*adot(i,j)*vy(j)
        enddo
      enddo
      yl = sqrt(yl)
    endif
  endif
  xymax = max(xl,yl)
  xscale = xl/xymax
  yscale = yl/xymax

! compute step length.

  do i=1,3
    dx(i) = vx(i)/(nx-1)
    dy(i) = vy(i)/(ny-1)
  enddo

  return

end subroutine plot_get_plane
