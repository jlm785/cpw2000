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

!>  Calculates the density of states with a linear interpolation in tetrahedra
!>  openmp instructions

  subroutine dos_lin(nx,egrid,nhist,ehist,dhist,chist,ezero,ispin,  &
               lunif,lper,lidos,vol,mxdbnd)

! Written November 19, 2013.
! Modified, documentation, 19 September 2020. JLM
! Modified to call dos_lin_one, dos_jminmax. 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)             ::  mxdbnd                             !<  size of number of bands
  integer, intent(in)             ::  nx(3)                              !<  number of k-points in each direction in regular grid
  integer, intent(in)             ::  nhist                              !<  number of points in histogram

  real(REAL64), intent(in)        ::  egrid(nx(1),nx(2),nx(3),mxdbnd)    !<  eigenvalues in Hartree in regular grid
  logical, intent(in)             ::  lunif                              !<  true if grid is uniform
  logical, intent(in)             ::  lper                               !<  true if grid is periodic
  logical, intent(in)             ::  lidos                              !<  true if integrated density of states is to be computed.
  real(REAL64), intent(in)        ::  ehist(nhist)                       !<  histogram energy
  real(REAL64), intent(in)        ::  vol                                !<  volume (in reciprocal space) of each element
  real(REAL64), intent(in)        ::  ezero                              !<  zero of energy
  integer, intent(in)             ::  ispin                              !<  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)

! input and output

  real(REAL64), intent(out)       ::  dhist(nhist)                       !<  histogram dos
  real(REAL64), intent(out)       ::  chist(nhist)                       !<  histogram intergrated dos)

! local arrays

  real(REAL64), allocatable  ::  dhistloc(:)
  real(REAL64), allocatable  ::  chistloc(:)
  real(REAL64), allocatable  ::  chlowloc(:)
  real(REAL64), allocatable  ::  chlow(:)

! constants

  real(REAL64), parameter ::  ZERO = 0.0_REAL64

! counters

  integer   ::  j, m
  integer   ::  jmin,jmax


! loop over k points and bands

  if(lidos) then

    jmin = 1
    jmax = mxdbnd

  else

    call dos_jminmax(nx,egrid,nhist,ehist,ezero,lper,jmin,jmax,mxdbnd)

  endif

  dhist(:) = ZERO
  chist(:) = ZERO

  allocate(dhistloc(nhist))
  allocate(chistloc(nhist))
  allocate(chlowloc(nhist))
  
  allocate(chlow(nhist))

  chlow(:) = ZERO

  do j = jmin,jmax

    call dos_lin_one(nx,egrid(:,:,:,j),nhist,ehist,dhistloc,chistloc,chlowloc,       &
               ezero,ispin,lunif,lper,lidos,vol)

    do m = 1,nhist
      dhist(m) = dhist(m) + dhistloc(m)
      chist(m) = chist(m) + chistloc(m)
      chlow(m) = chlow(m) + chlowloc(m)
    enddo

  enddo

  do m = 1,nhist
    do j = m,nhist
      chist(j) = chist(j) + chlow(m)
    enddo
  enddo

  deallocate(dhistloc)
  deallocate(chistloc)
  deallocate(chlowloc)
  deallocate(chlow)

  return
end subroutine dos_lin
