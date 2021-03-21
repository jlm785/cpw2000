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

!>  Provides a very ROUGH estimate of the optical gap and square oscillator strength

subroutine opt_rad_rough_pcvsq(filegrid, io_grid, egapopt, pcvsq, nvtc, nx,ny,nz)


! Written 16 October 2020. JLM
! Modified, use of files. 14 December 2020. JLM
! copyright  Jose Luis Martins/INESC-MN

! version 4.99

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)
  integer, parameter          ::  REAL32 = selected_real_kind(6)

! input

  integer, intent(in)                ::  nx,ny,nz                        !<  grid size
  integer, intent(in)                ::  nvtc                            !<  number of valence bands times number of conduction bands

  character(len=16), intent(in)      ::  filegrid                        !<  file with grid data
  integer, intent(in)                ::  io_grid                         !<  tape number for grid

! output

  real(REAL64), intent(out)          ::  egapopt                         !<  optical gap
  real(REAL64), intent(out)          ::  pcvsq                           !<  pcv^2

! local allocatable arrays

  real(REAL32), allocatable          ::  egrid_32(:,:,:)                 !  excitation values in one "joint band"
  real(REAL32), allocatable          ::  fkgrid_32(:,:,:)                !  dipole matrix elements in one "joint band"

  real(REAL32), allocatable          ::  xmin_32(:)                      !  minimum energy for band pair  


! local variables

  real(REAL64)           ::  emin
  real(REAL64)           ::  fkmax

  integer                ::  irec_err, ir_size

! constants

  real(REAL64), parameter      ::  ZERO = 0.0_REAL64
  real(REAL64), parameter      ::  EPS = 0.0001_REAL64

! counters

  integer                      ::  i, j, k, m


  allocate(xmin_32(nvtc))

  allocate(egrid_32(nx,ny,nz))
  allocate(fkgrid_32(nx,ny,nz))

  inquire(iolength = ir_size) egrid_32(:,:,:), fkgrid_32(:,:,:)
 
  open(unit = io_grid, file = trim(filegrid), access="direct", recl=ir_size)

! finds minimum excitation energy

  do m = 1,nvtc

    read(unit = io_grid, rec = m, iostat=irec_err) egrid_32(:,:,:)

    xmin_32(m) = egrid_32(1,1,1)
    do k = 1,nz
    do j = 1,ny
    do i = 1,nx
      if(xmin_32(m) > egrid_32(i,j,k)) xmin_32(m) = egrid_32(i,j,k)
    enddo
    enddo
    enddo
  enddo

  emin = xmin_32(1)
  do m = 1,nvtc
    if(emin > xmin_32(m)) emin = xmin_32(m)
  enddo

  egapopt = emin

  emin = emin + EPS
  fkmax = ZERO

! finds maximum dipole square among states close to the gap (quasi-degeneracy)

  do m = 1,nvtc

    read(unit = io_grid, rec = m, iostat=irec_err) egrid_32(:,:,:), fkgrid_32(:,:,:)

    do k = 1,nz
    do j = 1,ny
    do i = 1,nx

      if(egrid_32(i,j,k) < emin) then
        if(fkmax < fkgrid_32(i,j,k)) fkmax = fkgrid_32(i,j,k)
      endif

    enddo
    enddo
    enddo
  enddo

  pcvsq = fkmax
  
  close(unit = io_grid)

  deallocate(egrid_32)
  deallocate(fkgrid_32)

  return
end subroutine opt_rad_rough_pcvsq
