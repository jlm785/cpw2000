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

!>  Generates a set reciprocal space directions for the alcagoita representation
!>  representation of effective masses.  Based on a octahedral
!>
!>  \author       Jose Luis Martins
!>  \version      5.10
!>  \date         20 December 2023.
!>  \copyright    GNU Public License v2

subroutine berry_mass_directions(nlin, adot, npoint, xpoint)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  nlin                            !<  number of points in the octahedron edge.
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

! output

  integer, intent(out)               ::  npoint                          !<  number of points generated
  real(REAL64), intent(out)          ::  xpoint(3,4*nlin*nlin)           !<  directions as points on unit sphere

! local variables

  integer         ::  ivert(3,6)                            !  octahedron vertices
  integer         ::  iedge(2,12)                           !  octahedron edges
  integer         ::  iface(3,8)                            !  octahedron faces

  real(REAL64)    ::  xx, yy
  real(REAL64)    ::  avec(3,3), bvec(3,3)
  real(REAL64)    ::  xk(3)

! constants

  real(REAL64)    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer         ::  i, j, k, m, icount


  if(nlin < 3) then
    write(6,*)
    write(6,*) '   Stopped in berry_mass_directions'
    write(6,*) '   nlin does not make sense ',nlin
    write(6,*)

    stop

  endif

  ivert(:,:) = 0
  do i = 1,3
    ivert(i,i) = 1
    ivert(i,i+3) = -1
  enddo

  icount = 0
  do j = 1,5
  do k = j+1,6
    if(ivert(1,j) + ivert(1,k) /=0 .or. ivert(2,j) + ivert(2,k) /=0 .or.     &
       ivert(3,j) + ivert(3,k) /=0) then
      icount = icount + 1
      iedge(1,icount) = j
      iedge(2,icount) = k
    endif
  enddo
  enddo

  icount = 0
  do i = 1,4
  do j = i+1,5
    if(ivert(1,i) + ivert(1,j) /=0 .or. ivert(2,i) + ivert(2,j) /=0 .or.       &
       ivert(3,i) + ivert(3,j) /=0) then
      do k = j+1,6
        if((ivert(1,j) + ivert(1,k) /=0 .or. ivert(2,j) + ivert(2,k) /=0 .or.  &
            ivert(3,j) + ivert(3,k) /=0)   .and.                               &
           (ivert(1,i) + ivert(1,k) /=0 .or. ivert(2,i) + ivert(2,k) /=0 .or.  &
            ivert(3,i) + ivert(3,k) /=0)) then
          icount = icount + 1
          iface(1,icount) = i
          iface(2,icount) = j
          iface(3,icount) = k
        endif
      enddo
    endif
  enddo
  enddo

! vertices

  do i = 1,6
  do j = 1,3
    xpoint(j,i) = UM*ivert(j,i)
  enddo
  enddo

  icount = 6

! edges

  do i = 1,12
    do k = 1,nlin-2
      xx = (UM*k) / (UM*(nlin-1))
      icount = icount+1
      do j = 1,3
        xpoint(j,icount) = xx*ivert(j,iedge(1,i)) + (UM-xx)*ivert(j,iedge(2,i))
      enddo
    enddo
  enddo

! faces

  if(nlin > 3) then
    do i = 1,8
      do k = 1,nlin-3
        if(nlin-2-k > 0) then
          do m = 1,nlin-2-k
            xx = (UM*k) / (UM*(nlin-1))
            yy = (UM*m) / (UM*(nlin-1))
            icount = icount+1
            do j = 1,3
              xpoint(j,icount) = xx*ivert(j,iface(1,i)) + yy*ivert(j,iface(2,i)) +   &
                                 (UM-xx-yy)*ivert(j,iface(3,i))
            enddo
          enddo
        endif
      enddo
    enddo
  endif

  npoint = icount

! points on a sphere

  do i =1,npoint
    xx = xpoint(1,i)*xpoint(1,i) + xpoint(2,i)*xpoint(2,i)+ xpoint(3,i)*xpoint(3,i)
    xx = UM / sqrt(xx)
    do j = 1,3
      xpoint(j,i) = xx*xpoint(j,i)
    enddo
  enddo

! converts to lattice coordinates

  call adot_to_avec_sym(adot,avec,bvec)

  do i = 1,npoint
    do j = 1,3
      xk(j) = ZERO
      do k = 1,3
        xk(j) = xk(j) + xpoint(k,i)*avec(k,j)
      enddo
    enddo
    do j = 1,3
      xpoint(j,i) = xk(j)
    enddo
  enddo


  return

end subroutine berry_mass_directions
