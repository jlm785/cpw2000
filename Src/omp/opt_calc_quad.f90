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

!>  Performs the calculation of the dielectric tensor (independent electron approximation)
!>  and related optical functions, and calls subroutines for output file printing.


subroutine opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehist, e_re, e_im,    &
                       nx,ny,nz, nhist, nvtc)

! Written by Carlos Loia Reis. July 2020
! Modified, documentation, 20 September 2020. JLM
! Modified, output, egrid, 18 October 2020. JLM
! Modified, uses files to avoid exceeding RAM. 12 december 2020. JLM
! copyright  Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL32 = selected_real_kind(6)

! input

  integer, intent(in)                ::  nx,ny,nz                        !<  grid size
  integer, intent(in)                ::  nhist                           !<  energy values in histograms
  integer, intent(in)                ::  nvtc                            !<  number of valence bands times number of conduction bands

  integer, intent(in)                ::  ispin                             !<  ispin = 1 with spin-orbit, ispin = 2 no spin-orbit
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in real space

  character(len=16), intent(in)      ::  filegrid                        !<  file with grid data
  integer, intent(in)                ::  io_grid                         !<  tape number for grid

  integer, intent(in)                ::  ninter                          !<  number of interpolated points for big grid (~4)
  real(REAL64), intent(in)           ::  ehist(nhist)                    !<  histogram energy

! output

  real(REAL64), intent(out)          ::  e_re(nhist), e_im(nhist)        !<  real and imaginary parts of the dielectric function

! local allocatable arrays

  real(REAL64), allocatable          ::  dhist(:)

  real(REAL64), allocatable          ::  e_big(:,:,:)
  real(REAL64), allocatable          ::  f_big(:,:,:)

  real(REAL64), allocatable          ::  egrid(:,:,:)                    !  excitation values in one "joint band"
  real(REAL64), allocatable          ::  fkgrid(:,:,:)                   !  dipole matrix elements in one "joint band"

  real(REAL32), allocatable          ::  egrid_32(:,:,:)                 !  excitation values in one "joint band"
  real(REAL32), allocatable          ::  fkgrid_32(:,:,:)                !  dipole matrix elements in one "joint band"

! local variables

  real(REAL64)                 ::  eint                                  !  interpolated bands
  integer                      ::  nxx(3)
  real(REAL64)                 ::  rk(3)

  integer                      ::  irec_err, ir_size

  integer                      ::  nkgrid(3)

  real(REAL64)                 ::  vcell, bdot(3,3)
  real(REAL64)                 ::  vol
  logical                      ::  lunif

! constants

  real(REAL64), parameter      :: ZERO = 0.0_REAL64
  real(REAL64), parameter      :: UM = 1.0_REAL64
  real(REAL64), parameter      :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter      :: FOUR_PI_SQ = 4*PI*PI
  real(REAL64), parameter      :: EPS = 1.0E-08_REAL64

! counters

  integer          :: i, j, k, iband

  nkgrid(1) = nx
  nkgrid(2) = ny
  nkgrid(3) = nz

  e_im(:) = ZERO
  e_re(:) = ZERO

  allocate(dhist(nhist))

  dhist(:) = ZERO

  allocate(egrid_32(nx,ny,nz))
  allocate(fkgrid_32(nx,ny,nz))

  allocate(egrid(nx,ny,nz))
  allocate(fkgrid(nx,ny,nz))

  inquire(iolength = ir_size) egrid_32(:,:,:), fkgrid_32(:,:,:)

  open(unit = io_grid, file = trim(filegrid), access="direct", recl=ir_size)

  if(ninter == 1) then

    lunif = .TRUE.
    vol = UM / (nkgrid(1)*nkgrid(2)*nkgrid(3))

    do iband = 1,nvtc

      read(unit = io_grid, rec = iband, iostat=irec_err) egrid_32(:,:,:), fkgrid_32(:,:,:)

      egrid(:,:,:) = egrid_32(:,:,:)
      fkgrid(:,:,:) = fkgrid_32(:,:,:)

      call dosf_lin_one(nkgrid, egrid, fkgrid, nhist, ehist, dhist,      &
                        ZERO, ispin, lunif, .TRUE., vol)

      do i = 1,nhist
        if ( abs(ehist(i)) > EPS) then
           e_im(i) = e_im(i) + dhist(i)/(ehist(i)*ehist(i))
        endif
      enddo

      call progress(iband,nvtc)
      call flush(6)

    enddo

  else

    nxx(1) = ninter*nkgrid(1)
    nxx(2) = ninter*nkgrid(2)
    nxx(3) = ninter*nkgrid(3)

    allocate(e_big(nxx(1),nxx(2),nxx(3)))
    allocate(f_big(nxx(1),nxx(2),nxx(3)))

    lunif = .TRUE.
    vol = UM / (nxx(1)*nxx(2)*nxx(3))

    do iband = 1,nvtc

! interpolates energy and function

      read(unit = io_grid, rec = iband, iostat=irec_err) egrid_32(:,:,:), fkgrid_32(:,:,:)

      egrid(:,:,:) = egrid_32(:,:,:)
      fkgrid(:,:,:) = fkgrid_32(:,:,:)

!$omp   parallel do default(shared) private(i,j,k,rk,eint)
      do k = 1,nxx(3)
      do j = 1,nxx(2)
      do i = 1,nxx(1)
        rk(1) = (UM*(i-1)) / ninter
        rk(2) = (UM*(j-1)) / ninter
        rk(3) = (UM*(k-1)) / ninter

        call dos_grid_quad(eint,rk,nkgrid,egrid,.TRUE.,nkgrid)

        e_big(i,j,k) = eint

        call dos_grid_quad(eint,rk,nkgrid,fkgrid,.TRUE.,nkgrid)

        f_big(i,j,k) = eint

      enddo
      enddo
      enddo
!$omp end parallel do

      call dosf_lin_one(nxx, e_big, f_big,  nhist, ehist, dhist,         &
                        ZERO, ispin, lunif, .TRUE., vol)

      do i = 1,nhist
        if ( abs(ehist(i)) > EPS) then
           e_im(i) = e_im(i) + dhist(i)/(ehist(i)*ehist(i))
        endif
      enddo

      call progress(iband,nvtc)
      call flush(6)

    enddo
  endif

  deallocate(egrid_32)
  deallocate(fkgrid_32)

  deallocate(egrid)
  deallocate(fkgrid)

  deallocate(dhist)
  if(ninter /= 1) then
    deallocate(e_big)
    deallocate(f_big)
  endif

  close(unit = io_grid)

  write(6,*)
  write(6,*)

  call adot_to_bdot(adot, vcell, bdot)

! calculates real part and rescales

  do i = 1,nhist
    e_im(i) = e_im(i)*FOUR_PI_SQ/vcell
  enddo

  call opt_kk_im(ehist, e_im, e_re, nhist)

  return
end subroutine opt_calc_quad

