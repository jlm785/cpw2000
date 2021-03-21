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

!>  Weights for super-quadratic interpolation on a 3D grid using 32 points.
!>  It is exact for quadratic functions. It is also exact for
!>  some higher order monomials if: none of the coordinates has a
!>  power of three or higher and only one coordinate at most has a power of two.

recursive subroutine quad_3D32pt(x,fl)

! Written 30 November 2013. JLM 
! Modified documentation August 2019.  JLM
! Modified to output weights, and declared recursive 
! to make it thread-safe. 17 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)        ::  x(3)                               !<  point to be interpolated in grid. Should be 0 < d(i) < 1.

! output

  real(REAL64), intent(out)       ::  fl(-1:2,-1:2,-1:2)                 !<  function on a 4x4x4 grid

! local variables

  real(REAL64)        ::  d(3)                                      !  point to be interpolated in grid. Should be 0 < d(i) < 1.
  real(REAL64)        ::  fg(-1:1,-1:1,-1:1)                        !  function weights on the 3x3x3 grid

! counters

  integer  ::  i1,i2,i3, j1,j2,j3, k1,k2,k3


! constants

  real(REAL64), parameter ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64


! traps for out of bounds  uncomment the following lines if you want to check
! input values, the result would be a thread unsafe subroutine

!   logical    ::  lfirst = .TRUE.
!   integer    ::  nout = 0
!   real(REAL64), parameter ::  EPS = 0.001_REAL64
!   
! checks bounds
! 
!   if(lfirst) then
! 
!     if(x(1) < -EPS .or. x(1) > UM+EPS .or.                          &
! &      x(2) < -EPS .or. x(2) > UM+EPS .or.                          &
! &      x(3) < -EPS .or. x(3) > UM+EPS) THEN
! 
!       nout = nout + 1
!         write(6,*)
!         write(6,*) '     WARNING'
!         write(6,*)
!         write(6,*) '  quad_3D32pt called with bad arguments'
!         write(6,'("   0 < ",3f8.2," < 1  ")') x(1),x(2),x(3)
!         write(6,*)
!                  
!       if(nout >= 10) then
!         write(6,*)
!         write(6,*) '  TOO MANY WARNINGS.  WILL NOT WARN AGAIN'
!         write(6,*)
!         lfirst = .FALSE.
!       endif
! 
!     endif
! 
!   endif

! loop over 3x3x3 interpolation direction

  
  
  fl(:,:,:) = ZERO

  do i1 = -1,1,2
  do i2 = -1,1,2
  do i3 = -1,1,2

    d(1) = i1*x(1) + (1-i1)/2
    d(2) = i2*x(2) + (1-i2)/2
    d(3) = i3*x(3) + (1-i3)/2

!   loop over 3x3x3 grid points

!   call quad_3D10pt(d,fg)

    fg(:,:,:) = ZERO

    fg( 0, 0, 0) = UM + d(1)*d(2) + d(2)*d(3) + d(3)*d(1) - d(1)*d(1) - d(2)*d(2) - d(3)*d(3)

    fg(-1, 0, 0) = d(1) * ( d(1) - UM ) / 2
    fg( 0,-1, 0) = d(2) * ( d(2) - UM ) / 2
    fg( 0, 0,-1) = d(3) * ( d(3) - UM ) / 2

    fg( 1, 0, 0) = d(1) * ( UM + d(1) - 2*d(2) -2*d(3) ) / 2
    fg( 0, 1, 0) = d(2) * ( UM + d(2) - 2*d(3) -2*d(1) ) / 2
    fg( 0, 0, 1) = d(3) * ( UM + d(3) - 2*d(1) -2*d(2) ) / 2

    fg( 0, 1, 1) = d(2) * d(3)
    fg( 1, 0, 1) = d(3) * d(1)
    fg( 1, 1, 0) = d(1) * d(2)

    do j1 = -1,1
    do j2 = -1,1
    do j3 = -1,1

      k1 = i1*j1 + (1-i1)/2 
      k2 = i2*j2 + (1-i2)/2 
      k3 = i3*j3 + (1-i3)/2

      fl(k1,k2,k3) = fl(k1,k2,k3) + (UM-d(1))*(UM-d(2))*(UM-d(3)) * fg(j1,j2,j3)

    enddo
    enddo
    enddo

  enddo
  enddo
  enddo

  return
end subroutine quad_3D32pt
