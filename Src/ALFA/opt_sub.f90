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

!>  Calculates the dielectric tensor (independent electron approximation)
!>  and related optical functions.  master subroutine

subroutine opt_sub(ioreplay)

! Written by Carlos Loia Reis. July 2020
! Modified, documentation, 20 September 2020. JLM
! Modified, rearrange subroutines, October 23 2020. JLM
! copyright  Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter                 :: REAL64 = selected_real_kind(12)
  integer, parameter                 :: REAL32 = selected_real_kind(6)

! input

  integer, intent(in)                :: ioreplay                         !<  tape number for reproducing calculations

! local allocatable variables

  real(REAL64), allocatable          ::  el(:,:)                         !  eigenvalues in Hartree in irreducible wedge

  integer, allocatable               ::  kmap(:,:,:,:)                   !  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, allocatable               ::  indk(:,:)                       !  index of the six k-points neighbouring k-point i
  real(REAL64),allocatable           ::  rk(:,:)                         !  component in lattice coordinates of the k-point in the mesh
  integer, allocatable               ::  nband(:)                        !  number of bands for each k-points
  real(REAL64),allocatable           ::  wgk(:)                          !  weight in the integration of k-point (not used)

  real(REAL64), allocatable          ::  ehist(:)                        !  histogram energy
  real(REAL64), allocatable          ::  e_re(:), e_im(:)                !  real and imaginary parts of the dielectric function

  real(REAL64), allocatable          ::  e_of_k(:,:)                     !  band energies of k-point in plot
  real(REAL64), allocatable          ::  e_of_k_so(:,:)                  !  spin-orbit band energies of k-point in plot

! local variables                                                       

  integer                            ::  mxdbnd                         !  array dimension for the number of bands

  integer                            ::  nx, ny ,nz                     !  grid size

  integer                            ::  ninter                         !  number of interpolated points for big grid (~4)
  integer                            ::  nhist                          !  number of energy values in histogram

  integer                            ::  ntrans                         !  number of symmetry operations in the factor group
  integer                            ::  mtrx(3,3,48)                   !  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer                            ::  nrk                            !  number of k-points for integration in the irreducible wedge of the brillouin zone

  real(REAL64)                       ::  ztot                           !  total number of electrons
  integer                            ::  neig                           !  number of bands
  integer                            ::  nval                           !  number of valence bands
  integer                            ::  ncond                          !  number of conduction bands
  integer                            ::  nvtc                           !  nval*ncond or...

  real(REAL64)                       :: adot(3,3)                       !  metric in real space

  real(REAL64)                       ::  vcell, bdot(3,3)

  integer                            ::  ix,iy
  character(len = 1)                 ::  yesno
  integer                            ::  ispin
  integer                            ::  ix_in, iy_in
  integer                            ::  io, iotape

  character(len = 14)                ::  filename                        !  input file name

  character(len=50)                  ::  title                           !  title for plots
  character(len=140)                 ::  subtitle                        !  subtitle for plots

  logical                            ::  lscl                            !  write for scalar (no-spin-orbit)
  logical                            ::  lso                             !  also write for spin-orbit
  integer                            ::  identif                         !  identifier that is almost unique

  integer                            ::  nhtarg                          !  target number of points in  histogram
  integer                            ::  nkgrid(3)
  real(REAL64)                       ::  ezero, deltae

  integer                            ::  io_unsym
  integer                            ::  io_dhdrk
  integer                            ::  io_tmp
  integer                            ::  io_grid

  character(len=16)                  ::  filedhdrk
  character(len=16)                  ::  fileunsym
  character(len=16)                  ::  filetmp
  character(len=16)                  ::  filegrid

! constants

  real(REAL64), parameter      :: ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter      :: HARTREE = 27.21138386_REAL64

! counters

  integer            :: i, irk


! parameters for the calculation

  filename = "dos_file.dat"
  filetmp = "tmp_tmp.dat"
  fileunsym = "tmp_unsym.dat"
  filegrid = "tmp_grid.dat"

  io = 21
  iotape = 15

  io_dhdrk = 80
  io_unsym = 81
  io_tmp = 82
  io_grid = 83

! these parameters are hardcoded.  Change them if you need more precision

  nhtarg = 2000

  write(6,*)
  write(6,*) '   Dielectric tensor calculation'
  write(6,*) 
  write(6,*) '   which component pair do you want to compute?'
  write(6,*)
  write(6,*) '   (1,1) = xx, (2,2) = yy , (3,3) = zz, etc ...'
  write(6,*) '   (0,0)    computes scalar '
  write(6,*) '   (-1,-1)  computes diagonal '
  write(6,*) '   (-2,-2)  computes all '
  write(6,*) '    enter i,j pair '

  read(5,*) ix_in, iy_in 
  write(ioreplay,*) ix_in, iy_in,'     enter i,j pair'

  write(6,*)
  write(6,*) '   do you want results including Spin-Orbit (y/n) ?'
  read(5,*) yesno
  write(ioreplay,*) yesno,'     spin-orbit for optical'

  if(yesno == 'y' .or. yesno == 'Y') then
    ispin   = 1
  else
    ispin   = 2
  endif

  write(6,*)
  write(6,*) '   How many interpolation points do you want (~4)?'
  write(6,*) '   Enter 1 for no interpolation.'
  write(6,*) '   Do not use interpolation for perfect superlattices.'
  read(5,*) ninter
  write(ioreplay,*) ninter,'     number of superquadratic interpolation points'

  if(ninter < 1) then
    ninter = 1
    write(6,*)
    write(6,*) '   wrong value, ninter set to ', ninter
    write(6,*)
  endif

  if(ninter > 10) then
    ninter = 10
    write(6,*)
    write(6,*) '   too large value, ninter set to ', ninter
    write(6,*)
  endif


! reads data file


  call dos_read_size(trim(filename), io, nrk, mxdbnd, nkgrid)

  nx = nkgrid(1)
  ny = nkgrid(2)
  nz = nkgrid(3)

  allocate(e_of_k(mxdbnd,nrk))
  allocate(e_of_k_so(2*mxdbnd,nrk))
  allocate(nband(nrk))
  allocate(wgk(nrk))
  allocate(rk(3,nrk))
  allocate(indk(6,nrk))
  allocate(kmap(3,nx,ny,nz))

  call dos_read_data(trim(filename), io, title, subtitle,                &
    lscl, lso, identif,                                                  &
    nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                           &
    nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                       &
    nrk, mxdbnd)

  call adot_to_bdot(adot,vcell,bdot)

  neig = nband(1)
  do irk = 1,nrk
    if(neig > nband(irk)) neig = nband(irk)
  enddo

  write(6,*)
  write(6,*) '   The title and subtitle associated with the data file are:' 
  write(6,*)
  write(6,*) '   ',title
  write(6,*) '   ',subtitle
  write(6,*)

  if(ispin == 1) then
    if(.not. lso) then
      write(6,*) '    returning in opt_sub:  file does not contain spin-orbit info'

      return

    endif
    mxdbnd = 2*mxdbnd
    neig = 2*neig
    nval = nint(ztot)
    ncond = neig - nval
    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k_so(:,:)
    do irk = 1,nrk
      nband(irk) = 2*nband(irk)
    enddo

    filedhdrk = 'opt_dhdrk_so.dat'

    write(6,*)
    write(6,*) '     Optical functions calculated with spin-orbit'
    write(6,*)

  else
    if(.not. lscl) then
      write(6,*) '    returning in dos_sub:  file does not contain scalar info'

      return

    endif

    nval = nint(0.5*ztot)
    ncond = neig - nval

    allocate(el(mxdbnd,nrk))
    el(:,:) = e_of_k(:,:)

    filedhdrk = 'opt_dhdrk.dat'

    write(6,*)
    write(6,*) '     Optical functions calculated WITHOUT spin-orbit'
    write(6,*)

  endif


! generates unsymmetrized grid

  call opt_grid_unsym(el, filedhdrk, io_dhdrk, fileunsym, io_unsym, kmap, mtrx, &
    neig, nx,ny,nz, nrk, mxdbnd)

  deallocate(kmap)
  deallocate(indk)
  deallocate(rk)
  deallocate(wgk)


! sets histogram range

  call opt_set_in_grid(0,0, adot, fileunsym, io_unsym, filetmp, io_tmp,       &
                       filegrid, io_grid, neig, nval, ncond, nx, ny, nz)


  ezero = ZERO

  call opt_plot_range(el, neig, nval, nhtarg, deltae, nhist, nrk, mxdbnd)

  write(6,*)
  write(6,'("  The optical functions will be plotted up to ",2x,f10.3," eV")') nhist*deltae*HARTREE
  write(6,*)

  allocate(ehist(nhist))
  allocate(e_im(nhist))
  allocate(e_re(nhist))

  do i = 1,nhist
    ehist(i) = deltae*i
  enddo

! does the job

  nvtc = nval*ncond

  if(ix_in == -1 .and. iy_in == -1) then

    do ix = 1,3

      call opt_set_in_grid(ix, ix, adot, fileunsym, io_unsym, filetmp, io_tmp,   &
                           filegrid, io_grid, neig, nval, ncond, nx, ny, nz)

      call opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehist,  &
                         e_re, e_im, nx,ny,nz, nhist, nvtc)

      call opt_write(ix, iy, iotape, title, subtitle,                    &
                     nhist, ehist, e_re, e_im, vcell, ztot, ispin)

    enddo

  elseif(ix_in == 0 .and. iy_in == 0) then

    ix = 0
    iy = 0

    call opt_set_in_grid(ix, iy, adot, fileunsym, io_unsym, filetmp, io_tmp,     &
                         filegrid, io_grid, neig, nval, ncond, nx, ny, nz)

    call opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehist,    &
                      e_re, e_im, nx,ny,nz, nhist, nvtc)

    call opt_write(ix, iy, iotape, title, subtitle,                      &
                   nhist, ehist, e_re, e_im, vcell, ztot, ispin)

  elseif(ix_in == -2 .and. iy_in == -2) then

    do ix = 1,3
    do iy = 1,3

      call opt_set_in_grid(ix, iy, adot, fileunsym, io_unsym, filetmp, io_tmp,   &
                           filegrid, io_grid, neig, nval, ncond, nx, ny, nz)

      call opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehist,  &
                         e_re, e_im, nx,ny,nz, nhist, nvtc)

      call opt_write(ix, iy, iotape, title, subtitle,                    &
                     nhist, ehist, e_re, e_im, vcell, ztot, ispin)

    enddo
    enddo

  elseif(ix_in > 0 .and. ix_in < 4 .and. iy_in > 0 .and. iy_in < 4) then

    ix = ix_in
    iy = iy_in

    call opt_set_in_grid(ix, iy, adot, fileunsym, io_unsym, filetmp, io_tmp,     &
                         filegrid, io_grid, neig, nval, ncond, nx, ny, nz)

    call opt_calc_quad(ispin, adot, filegrid, io_grid, ninter, ehist,    &
                        e_re, e_im, nx,ny,nz, nhist, nvtc)

    call opt_write(ix, iy, iotape, title, subtitle,                      &
                   nhist, ehist, e_re, e_im, vcell, ztot, ispin)

  else

    write(6,*) '  incorrect values of ix,iy, run again'

  endif

  deallocate(ehist)
  deallocate(e_re)
  deallocate(e_im)

  return
end subroutine opt_sub

