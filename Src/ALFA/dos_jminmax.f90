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

!>  Gets the range of bands within plotting range

subroutine dos_jminmax(nx,egrid,nhist,ehist,ezero,lper,jmin,jmax,mxdbnd)

! Written 18 October 2020 based on earlier code. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98  

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)             ::  mxdbnd                             !<  size of number of bands
  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in regular grid

  real(REAL64), intent(in)        ::  egrid(nx(1),nx(2),nx(3),mxdbnd)    !<  eigenvalues in Hartree in regular grid
  integer, intent(in)             ::  nhist                              !<  number of points in histogram
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  histogram energy
  real(REAL64), intent(in)        ::  ezero                              !<  zero of energy
  logical, intent(in)             ::  lper                               !<  true if grid is periodic

! output:

  integer, intent(out)            ::  jmin, jmax                         !<  range of bands

! local variables

  logical       ::  ljump
  integer       ::  ilast, jm

!  
! counters

  integer   ::  j
  integer   ::  i1, i2, i3


  if(lper) then
    ilast = 0
  else
    ilast = 1
  endif

  ljump = .FALSE.
  do j = 1,mxdbnd
    jmin = j
    do i3 = 1,nx(3)-ilast
      do i2 = 1,nx(2)-ilast
        do i1 = 1,nx(1)-ilast
          if(egrid(i1,i2,i3,j) - ezero > ehist(1) ) then
            ljump = .TRUE.

            exit

          endif
        enddo

        if(ljump) exit

      enddo

      if(ljump) exit

    enddo

    if(ljump) exit

  enddo

  ljump = .FALSE.
  jmax = 1
  do i3 = 1,nx(3)-ilast
  do i2 = 1,nx(2)-ilast
  do i1 = 1,nx(1)-ilast
    do j = mxdbnd,1,-1
      jm = j
      if(egrid(i1,i2,i3,j) - ezero < ehist(nhist) ) then
        exit

      endif
    enddo
    jmax = max(jm,jmax)
  enddo
  enddo
  enddo

  return
end subroutine dos_jminmax
