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

!>  Calculates the density of states and integrated dos.
!>  First it interpolates super-quadraticaly on a finer mesh,
!>  and then does the linear tetrahedron interpolation on that finer mesh
!>  uses openmp

  subroutine dos_quad(nx,egrid,ninter,nhist,ehist,dhist,chist,      &
        ezero,ispin,lunif,lper,lidos,mxdbnd)

! Written December 2, 2013.
! Modified, documentation, 19 September 2020. JLM
! Modified, egrid, dos_lin_one, dos_jminmax, 18 October 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.98 of cpw

  implicit none

  integer, parameter  :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)             ::  mxdbnd                        !<  size of number of bands
  integer, intent(in)             ::  nx(3)                         !<  number of k-points in each direction in regular grid
  integer, intent(in)             ::  nhist                         !<  number of points in histogram

  real(REAL64), intent(in)  ::  egrid(nx(1),nx(2),nx(3),mxdbnd)     !<  eigenvalues in Hartree in regular grid
  logical, intent(in)             ::  lunif                         !<  true if grid is uniform
  logical, intent(in)             ::  lper                          !<  true if grid is periodic
  logical, intent(in)             ::  lidos                         !<  true if integrated density of states is to be computed.
  real(REAL64), intent(in)        ::  ehist(nhist)                  !<  histogram energy
  integer, intent(in)             ::  ninter                        !<  number of interpolated points for big grid (~4)
  real(REAL64), intent(in)        ::  ezero                         !<  zero of energy
  integer, intent(in)             ::  ispin                         !<  spin degeneracy (2 for non-spin-polarized, 1 for spin-orbit)

! output

  real(REAL64), intent(out)       ::  dhist(nhist)                  !<  histogram dos
  real(REAL64), intent(out)       ::  chist(nhist)                  !<  histogram integrated dos

! local arrays

  real(REAL64), allocatable  ::  egbig(:,:,:)                   !  denser, interpolated eigenvalue array

  real(REAL64), allocatable  ::  dhistloc(:)
  real(REAL64), allocatable  ::  chistloc(:)
  real(REAL64), allocatable  ::  chlowloc(:)
  real(REAL64), allocatable  ::  chlow(:)

  real(REAL64), allocatable  ::  egridloc(:,:,:)

! local variables

  real(REAL64)   ::  eint                            !  interpolated bands
  integer        ::  nxx(3)
  real(REAL64)   ::  rk(3)
  real(REAL64)   ::  vol

! counters

  integer     ::  i, j, k, n, m
  integer     ::  jmin,jmax

! constants

  real(REAL64), parameter ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64


! loop over k points and bands

  if(lidos) then

    jmin = 1
    jmax = mxdbnd

  else

    call dos_jminmax(nx,egrid,nhist,ehist,ezero,lper,jmin,jmax,mxdbnd)

  endif

  dhist(:) = ZERO
  chist(:) = ZERO

  nxx(1) = ninter*nx(1)
  nxx(2) = ninter*nx(2)
  nxx(3) = ninter*nx(3)

  allocate(egbig(nxx(1),nxx(2),nxx(3)))
  allocate(egridloc(nx(1),nx(2),nx(3)))

  allocate(dhistloc(nhist))
  allocate(chistloc(nhist))
  allocate(chlowloc(nhist))
  
  allocate(chlow(nhist))

  chlow(:) = ZERO

  vol = UM / (nxx(1)*nxx(2)*nxx(3))

  do n = jmin,jmax

!   quad3D_32pt.f90 may not be thread safe.  If you need parallelism
!   remove the test from that subroutine and make sure you compile
!   everything at the same time, and double check.

!$omp   parallel do default(shared) private(i,j,k,rk,eint)
    do k = 1,nxx(3)
    do j = 1,nxx(2)
    do i = 1,nxx(1)
      rk(1) = (UM*(i-1)) / ninter
      rk(2) = (UM*(j-1)) / ninter
      rk(3) = (UM*(k-1)) / ninter

      call dos_grid_quad(eint,rk,nx,egrid(:,:,:,n),.TRUE.,nx)

      egbig(i,j,k) = eint

    enddo
    enddo
    enddo
!$omp end parallel do

    call dos_lin_one(nxx,egbig,nhist,ehist,dhistloc,chistloc,chlowloc,   &
               ezero,ispin,lunif,lper,lidos,vol)

    do m = 1,nhist
      dhist(m) = dhist(m) + dhistloc(m)
      chist(m) = chist(m) + chistloc(m)
      chlow(m) = chlow(m) + chlowloc(m)
    enddo

  enddo

!   do m = 1,nhist
!     do j = m,nhist
!       chist(j) = chist(j) + chlow(m)
!     enddo
!   enddo

  deallocate(egbig)
  deallocate(egridloc)

  deallocate(dhistloc)
  deallocate(chistloc)
  deallocate(chlowloc)
  deallocate(chlow)

  return
end subroutine dos_quad
