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

!>  Writes the file with the band information for later calculation of 
!>  the density of states or optical response

  subroutine dos_read_data(filename, io, title, subtitle,                &
    lscl, lso, identif,                                                  &
    nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                           &
    nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                       &
    mxdnrk, mxdbnd)


! Reverse of out_dos_write, 13 December 2020. JLM
! copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdnrk                          !<  size of k-points
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  integer, intent(in)                ::  nrk                             !<  number of irreducible k-points
  integer, intent(in)                ::  nx, ny, nz                      !<  original k-point mesh

  character(len=*), intent(in)       ::  filename                        !<  file that should be written
  integer, intent(in)                ::  io                              !<  tape numbers

! output

  character(len=50), intent(out)     ::  title                           !<  title for plots
  character(len=140), intent(out)    ::  subtitle                        !<  subtitle for plots

  logical, intent(out)               ::  lscl                            !<  write for scalar (no-spin-orbit)
  logical, intent(out)               ::  lso                             !<  also write for spin-orbit
  integer, intent(out)               ::  identif                         !<  identifier that is almost unique to the calculation


  real(REAL64), intent(out)          ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(out)          ::  adot(3,3)                       !<  metric in direct space
  integer, intent(out)               ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(out)               ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(out)               ::  nband(mxdnrk)                   !<  number of bands for each k-points
  real(REAL64), intent(out)          ::  rk(3,mxdnrk)                    !<  component in lattice coordinates of the k-point in the mesh
  real(REAL64), intent(out)          ::  wgk(mxdnrk)                     !<  weight in the integration of k-point
  integer, intent(out)               ::  indk(6,mxdnrk)                  !<  index of the six k-points neighbouring k-point i
  integer, intent(out)               ::  kmap(3,nx,ny,nz)                !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation

  real(REAL64), intent(out)          ::  e_of_k(mxdbnd,nrk)              !<  band energies of k-point in plot
  real(REAL64), intent(out)          ::  e_of_k_so(2*mxdbnd,nrk)         !<  spin-orbit band energies of k-point in plot

! local variables

  integer       ::  nx2,ny2,nz2
  integer       ::  nrk2
  integer       ::  ios

! counters

  integer                            ::  irk


  open(unit=io, file= trim(filename), status='OLD', form='UNFORMATTED', iostat=ios)

  if(ios /= 0) then
    write(6,*) '  Stopped in dos_read_data: ', trim(filename), ' not found'

    stop

  endif

  read(io) nrk2

  if(nrk2 /= nrk) then
    write(6,'("   stopped in dos_read_data:  nrk, nrk2 = ",2i8)') nrk, nrk2

    stop

  endif

  if(nrk > mxdnrk) then
    write(6,'("   stopped in dos_read_data:  nrk, mxdnrk = ",2i8)') nrk, mxdnrk

    stop

  endif

  read(io) nband(1:nrk)

  do irk = 1,nrk
    if(nband(irk) > mxdbnd) then
      write(6,'("   stopped in dos_read_data:  irk, nband(irk), mxdbnd = ",3i8)')  &
             irk, nband(irk), mxdbnd

      stop

    endif
  enddo

  read(io) nx2,ny2,nz2

  if(nx2 /= nx .or. ny2 /= ny .or. nz2 /= nz) then
    write(6,'("   stopped in dos_read_data:  inconsistent nxyz ",6i8)') nx,ny,nz, nx2,ny2,nz2

    stop

  endif

  read(io) title, subtitle
  read(io) lso, lscl, identif

  read(io) ztot
  read(io) adot

  read(io) ntrans
  read(io) mtrx

  read(io) rk(:,1:nrk)
  read(io) wgk(1:nrk)
  read(io) indk(:,1:nrk)
  read(io) kmap

  if(lscl) then
    do irk = 1, nrk
      read(io) e_of_k(1:nband(irk),irk)
    enddo
  endif

  if(lso) then
    do irk = 1, nrk
      read(io) e_of_k_so(1:2*nband(irk),irk)
    enddo
  endif

  close(unit=io)

  return
  end subroutine dos_read_data
