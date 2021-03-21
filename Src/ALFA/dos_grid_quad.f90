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

!>  Interpolates super-quadraticaly a 3D function defined either on a regular
!>  periodic grid or on a local grid.

recursive subroutine dos_grid_quad(eint,rk,nx,egrid,lper,nxmax)          

! Written 30 November 2013. JLM 
! Modified, documentation, 19 September 2020. JLM
! Modified to do one band at a time and declared recursive
! to make it thread safe, 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)             ::  nxmax(3)                           !<  maximum number of k-points in each direction in regular grid

  real(REAL64), intent(in)        ::  rk(3)                              !<  coordinates of point to be interpolated (unit is grid step, zero is first point)
  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in regular grid
  real(REAL64), intent(in)        ::  egrid(nxmax(1),nxmax(2),nxmax(3))  !<  eigenvalues in Hartree in regular grid
  logical, intent(in)             ::  lper                               !<  True if grid is periodic

! output

  real(REAL64), intent(out)       ::  eint                               !<  interpolated energies

! local variables

  integer             ::  kl(3)
  real(REAL64)        ::  x(3)                                           !  point to be interpolated in grid. Should be 0 < d(i) < 1.
  real(REAL64)        ::  fl(-1:2,-1:2,-1:2)                             !  weight of function on the 4x4x4 grid

! counters

  integer          ::  i, i1,i2,i3, j1,j2,j3

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64
  real(REAL64), parameter ::  EPS = 0.001_REAL64

! x should be between 0 an 1 (within error of EPS in the
! absence of periodicity. 
! kl takes the values 1,...,nx for the periodic case.
! kl takes the values 2,...,nx-2 for the non-periodic case .

! Does not test to be thread-safe.  Uncomment for general applications.

  if(lper) then

    do i=1,3
      kl(i) = floor(rk(i))
      x(i) = rk(i) - kl(i)
      kl(i) = mod(kl(i),nx(i))
      if(kl(i) < 0) kl(i) = kl(i) + nx(i)
      kl(i) = kl(i) + 1
    enddo

  else

    do i=1,3
      kl(i) = floor(rk(i))
      if(kl(i) == 0 .and. floor(rk(i)+EPS) == 1) kl(i) = 1
      if(kl(i) == (nx(i) - 2) .and. floor(rk(i)-EPS) == (nx(i) - 3) )    &
                    kl(i) = nx(i) - 3
      x(i) = rk(i) - kl(i)
      kl(i) = kl(i) + 1
!       if(kl(i) < 2)   stop   ' kl(i) < 2 '
!       if(kl(i) > nx(i) - 2)   stop   ' kl(i)> nx(i) - 2 '
    enddo

  endif

  call quad_3D32pt(x,fl)

! loop over 4x4x4 interpolation grid

  eint = ZERO

  do i1 = -1,2
    j1 = kl(1) + i1
    if(lper) then
      if(j1 <= 0) j1 = j1 + nx(1)
      if(j1 > nx(1)) j1 = j1 - nx(1)
    endif
    do i2 = -1,2
      j2 = kl(2) + i2
      if(lper) then
        if(j2 <= 0) j2 = j2 + nx(2)
        if(j2 > nx(2)) j2 = j2 - nx(2)
      endif
      do i3 = -1,2
        j3 = kl(3) + i3
        if(lper) then
          if(j3 <= 0) j3 = j3 + nx(3)
          if(j3 > nx(3)) j3 = j3 - nx(3)
        endif

        eint = eint + egrid(j1,j2,j3) * fl(i1,i2,i3)

      enddo
    enddo
  enddo


  return
end subroutine dos_grid_quad
