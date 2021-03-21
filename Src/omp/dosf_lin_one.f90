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

!>  calculates the density of states multiplied by a function
!>  with a linear interpolation in tetrahedra
!>  uses openmp instructions
!>
!>  \author       Jose Luis Martins, Carlos Loia Reis
!>  \version      4.99
!>  \date         2020-2021
!>  \copyright    GNU Public License v2

subroutine dosf_lin_one(nx,egrid,fkgrid,nhist,ehist,dhist,ezero,ispin,  &
               lunif,lper,vol)

! adapted from dos_lin
! Written November 2018.
! Modified, documentation, 19 September 2020. JLM
! Modified, iband, 20 October 2020. JLM
! copyright Carlos Loia Reis/ J.L.Martins, INESC-MN.

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in regular grid
  integer, intent(in)             ::  nhist                              !<  number of points in histogram

  real(REAL64), intent(in)        ::  egrid(nx(1),nx(2),nx(3))           !<  eigenvalues in Hartree in regular grid
  real(REAL64), intent(in)        ::  fkgrid(nx(1),nx(2),nx(3) )         !<  matrix in regular grid

  logical, intent(in)             ::  lunif                              !<  true if grid is uniform
  logical, intent(in)             ::  lper                               !<  true if grid is periodic
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  histogram energy
  real(REAL64), intent(in)        ::  vol                                !<  volume (in reciprocal space) of each element
  real(REAL64), intent(in)        ::  ezero                              !<  zero of energy
  integer, intent(in)             ::  ispin                              !<  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)

! input and output

  real(REAL64), intent(out)       ::  dhist(nhist)                       !<  histogram "dos*fk"

! local variables

  real(REAL64), allocatable  ::  ec(:), Fkc(:)                                        !  energies at corner of cube

! constants

  real(REAL64), parameter      :: ZERO = 0.0_REAL64

! counters

  integer   ::  j, ilast
  integer   ::  i1,i2,i3, ip1,ip2,ip3

! loop over k points and bands

  allocate(ec(8), Fkc(8))

  if(lper) then
    ilast = 0
  else
    ilast = 1
  endif

  dhist(:) = ZERO

!$omp   parallel do default(shared) private(i1,i2,i3,ip1,ip2,ip3,j,ec,Fkc)   &
!$omp&  reduction(+:dhist)
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

    ec(1) = egrid(i1,i2,i3)    - ezero
    ec(2) = egrid(ip1,i2,i3)   - ezero
    ec(3) = egrid(ip1,ip2,i3)  - ezero
    ec(4) = egrid(i1,ip2,i3)   - ezero
    ec(5) = egrid(i1,i2,ip3)   - ezero
    ec(6) = egrid(ip1,i2,ip3)  - ezero
    ec(7) = egrid(ip1,ip2,ip3) - ezero
    ec(8) = egrid(i1,ip2,ip3)  - ezero

    Fkc(1) = fkgrid(i1,i2,i3)
    Fkc(2) = fkgrid(ip1,i2,i3)
    Fkc(3) = fkgrid(ip1,ip2,i3)
    Fkc(4) = fkgrid(i1,ip2,i3)
    Fkc(5) = fkgrid(i1,i2,ip3)
    Fkc(6) = fkgrid(ip1,i2,ip3)
    Fkc(7) = fkgrid(ip1,ip2,ip3)
    Fkc(8) = fkgrid(i1,ip2,ip3)

    call dosf_cub(vol,ec,Fkc,nhist,ehist,dhist,lunif)

  enddo
  enddo
  enddo
!$omp end parallel do

  do j = 1,nhist
    dhist(j) = ispin*dhist(j)
  enddo

  return
end subroutine dosf_lin_one
