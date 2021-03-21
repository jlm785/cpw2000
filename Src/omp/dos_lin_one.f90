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

!>  calculates the density of states with a linear interpolation in tetrahedra
!>  for a single band
!>  openmp instructions

subroutine dos_lin_one(nx,egrid,nhist,ehist,dhist,chist,chlow,           &
                       ezero,ispin,lunif,lper,lidos,vol)

! Written November 19, 2013.
! Modified, documentation, 19 September 2020. JLM
! Modified to do one band at a time, efficient idos, 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)                   ::  nx(3)                        !<  number of k-points in each direction in regular grid
  integer, intent(in)                   ::  nhist                        !<  number of points in histogram

  real(REAL64), intent(in)              ::  egrid(nx(1),nx(2),nx(3))     !<  eigenvalues in Hartree in regular grid
  logical, intent(in)                   ::  lunif                        !<  true if grid is uniform
  logical, intent(in)                   ::  lper                         !<  true if grid is periodic
  logical, intent(in)                   ::  lidos                        !<  true if integrated density of states is to be computed.
  real(REAL64), intent(in)              ::  ehist(nhist)                 !<  histogram energy
  real(REAL64), intent(in)              ::  vol                          !<  volume (in reciprocal space) of each element
  real(REAL64), intent(in)              ::  ezero                        !<  zero of energy
  integer, intent(in)                   ::  ispin                        !<  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)

! output

  real(REAL64), intent(out)             ::  dhist(nhist)                 !<  histogram dos
  real(REAL64), intent(out)             ::  chist(nhist)                 !<  histogram intergrated dos)
  real(REAL64), intent(out)             ::  chlow(nhist)                 !<  add to chist for n >= j at the end 

! local variables

  real(REAL64)      ::  ec(8)                                            !  energies at corner of cube
  integer           ::  ilast

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64

! counters

  integer   ::  i1,i2,i3, ip1,ip2,ip3
  integer   ::  j

! loop over k points

  if(lper) then
    ilast = 0
  else
    ilast = 1
  endif

  dhist(:) = ZERO
  chist(:) = ZERO
  chlow(:) = ZERO

!$omp   parallel do default(shared) private(i1,i2,i3,ip1,ip2,ip3,j,ec)   & 
!$omp&  reduction(+:dhist,chist)
  do i3 = 1,nx(3)-ilast
  do i2 = 1,nx(2)-ilast
  do i1 = 1,nx(1)-ilast

    if(lper) then
      ip1 = mod(i1,nx(1)) + 1
      ip2 = mod(i2,nx(2)) + 1
      ip3 = mod(i3,nx(3)) + 1
    else
      ip1 = i1 + 1
      ip2 = i2 + 1
      ip3 = i3 + 1
    endif

!   bottom  4 3   top  8 7
!           1 2        5 6

    ec(1) = egrid(i1,i2,i3) - ezero
    ec(2) = egrid(ip1,i2,i3) - ezero
    ec(3) = egrid(ip1,ip2,i3) - ezero
    ec(4) = egrid(i1,ip2,i3) - ezero
    ec(5) = egrid(i1,i2,ip3) - ezero
    ec(6) = egrid(ip1,i2,ip3) - ezero
    ec(7) = egrid(ip1,ip2,ip3) - ezero
    ec(8) = egrid(i1,ip2,ip3) - ezero

    call dos_cub(vol,ec,nhist,ehist,dhist,chist,chlow,lunif,lidos)

  enddo
  enddo
  enddo
!$omp  end parallel do

  do j = 1,nhist
    dhist(j) = ispin*dhist(j)
    chist(j) = ispin*chist(j)
    chlow(j) = ispin*chlow(j)
  enddo

  return
end subroutine dos_lin_one
