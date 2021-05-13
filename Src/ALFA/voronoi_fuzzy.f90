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

!>  calculates the fuzzy Voronoy polyhedra weight "vpw"
!>  at a point r following the scheme of Becke (JChemPhys 88, 2547, 1988)
!>  sum vpw = 1, 0 < vpw < 1, and if r is close to atom i then vpw(i) is close to 1.
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         before 1994, modified 12 october 2006, modernized 20 April 2021
!>  \copyright    GNU Public License v2

subroutine voronoi_fuzzy(r, nneigh, rcar, rbs, vpw)

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)      ::  r(3)                                 !<  point where weight is evaluated (cartesian)

  integer, intent(in)           ::  nneigh                               !<  number of candidate neighbors
  real(REAL64), intent(in)      ::  rcar(3,nneigh)                       !<  cartesian coordinate of neighboring atoms
  real(REAL64), intent(in)      ::  rbs(nneigh)                          !<  size of atom i (Bragg-Slater table for example)

! output

  real(REAL64), intent(out)     ::  vpw(nneigh)                          !<  voronoi polyedra weight of atom i for point r.

! allocatable arrays

  real(REAL64), allocatable     ::  d(:), p(:)

! local variables

  real(REAL64)       ::  dij,xmuij,x,xuij,xaij

! parameters

  real(REAL64), parameter    ::  ONE = 1.0_REAL64
  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer     ::  i, j, k


  allocate(d(nneigh),p(nneigh))

  if(nneigh == 1) then

    vpw(1) = ONE

  else

    do i = 1,nneigh
      d(i) = sqrt((r(1)-rcar(1,i))*(r(1)-rcar(1,i)) +               &
                  (r(2)-rcar(2,i))*(r(2)-rcar(2,i)) +               &
                  (r(3)-rcar(3,i))*(r(3)-rcar(3,i)))
      p(i) = ONE
    enddo

    do i = 2,nneigh
    do j = 1,i-1
      dij = sqrt((rcar(1,j)-rcar(1,i))*(rcar(1,j)-rcar(1,i)) +      &
                 (rcar(2,j)-rcar(2,i))*(rcar(2,j)-rcar(2,i)) +      &
                 (rcar(3,j)-rcar(3,i))*(rcar(3,j)-rcar(3,i)))
      xmuij = (d(i)-d(j))/dij
      xuij = (rbs(i)-rbs(j)) / (rbs(i)+rbs(j))

      if(xuij > (sqrt(2*ONE)-ONE)) then
        xaij = -ONE/2
      elseif(xuij < -(sqrt(2*ONE)-ONE)) then
        xaij =  ONE/2
      else
        xaij = xuij / (xuij*xuij - ONE)
      endif

      x = xmuij + xaij * (ONE - xmuij*xmuij)
      do k = 1,3
        x = x*(3*ONE-x*x)/2
      enddo
      p(i) = p(i)*(ONE-x)/2

      x = -xmuij - xaij * (ONE - xmuij*xmuij)
      do k = 1,3
        x = x*(3*ONE-x*x)/2
      enddo
      p(j) = p(j)*(ONE-x)/2
    enddo
    enddo

    x = ZERO
    do i=1,nneigh
      x = x+p(i)
    enddo
    do i=1,nneigh
      vpw(i) = p(i)/x
    enddo

  endif

  deallocate(d,p)

  return
end subroutine voronoi_fuzzy
