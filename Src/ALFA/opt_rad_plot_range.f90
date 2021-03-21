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

!>  finds a reasonable scale for the energy range of
!>  the relevant joint density of sates and the relevant band pairs

  subroutine opt_rad_plot_range(filegrid, io_grid, nhtarg,               &
    erange, nrange, emin, deltae, nhist, nx,ny,nz, nvtc)

! Written December 9, 2013.
! Modified, documentation, 19 September 2020. JLM
! Modified, egrid, 18 October 2020. JLM
! Modified, using files, 14 December 2020. JLM
! copyright  J.L.Martins, INESC-MN.

! version 4.99 of cpw

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)
  integer, parameter          ::  REAL32 = selected_real_kind(6)

! input:

  integer, intent(in)                ::  nvtc                            !<  size of number of pairs of bands
  integer, intent(in)                ::  nx,ny,nz                        !<  number of k-points in each direction in regular grid

  character(len=16), intent(in)      ::  filegrid                        !<  file with grid data
  integer, intent(in)                ::  io_grid                         !<  tape number for grid

  integer, intent(in)                ::  nhtarg                          !<  approximate number of points for plot
  real(REAL64), intent(in)           ::  erange                          !<  do not consider excitation energies above optical band gap plus this value.

! output

  real(REAL64), intent(out)          ::  emin                            !<  minimum recomended energy for histogram 
  real(REAL64), intent(out)          ::  deltae                          !<  recomended energy step for histogram 
  integer, intent(out)               ::  nhist                           !<  recomended number of points for histogram
  integer, intent(out)               ::  nrange                          !<  number of band pairs within erange

! local allocatable arrays

  real(REAL32), allocatable          ::  egrid_32(:,:,:)                 !  excitation values in one "joint band"
  real(REAL32), allocatable          ::  fkgrid_32(:,:,:)                !  dipole matrix elements in one "joint band"

  real(REAL32), allocatable          ::  xmin_32(:)                         !  minimum energy for band pair  

! local variables

  real(REAL64)        ::  b, emax
  integer             ::  irec_err, ir_size

! counters

  integer   ::  i, j ,k, m

! constants

  real(REAL64), parameter :: HARTREE = 27.21138386_REAL64


! find emin and nrange

  allocate(xmin_32(nvtc))

  allocate(egrid_32(nx,ny,nz))
  allocate(fkgrid_32(nx,ny,nz))

  inquire(iolength = ir_size) egrid_32(:,:,:), fkgrid_32(:,:,:)
 
  open(unit = io_grid, file = trim(filegrid), access="direct", recl=ir_size)

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
  
  close(unit = io_grid)

  deallocate(egrid_32)
  deallocate(fkgrid_32)

  emin = xmin_32(1)
  do m = 1,nvtc
    if(emin > xmin_32(m)) emin = xmin_32(m)
  enddo

  nrange = 1
  do m = 1,nvtc
    if(xmin_32(m) < emin + erange) nrange = m
  enddo

  deallocate(xmin_32)

! find scale for the plot

  deltae = HARTREE*(erange) / nhtarg

  call plot_step(deltae,b)

  deltae = b/HARTREE

  emin = deltae*real(int((emin-erange/5)/deltae)-1)
  emax = emin + erange + deltae
  nhist = nint((emax-emin)/deltae) + 1

  return
end subroutine opt_rad_plot_range
