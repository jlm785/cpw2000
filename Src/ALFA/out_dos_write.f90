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

  subroutine out_dos_write(filename, io, title, subtitle,                &
    lscl, lso, identif,                                                  &
    nrk, nx, ny, nz, ztot, adot, ntrans, mtrx,                           &
    nband, rk, wgk, indk, kmap, e_of_k, e_of_k_so,                       &
    mxdnrk, mxdbnd)


! Merge of out_dos_write with out_opt_write, 13 December 2020. JLM
! copyright  Jose Luis Martins/Carlos Loia Reis/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdnrk                          !<  size of k-points
  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  character(len=*), intent(in)       ::  filename                        !<  file that should be written
  integer, intent(in)                ::  io                              !<  tape numbers

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  logical, intent(in)                ::  lscl                            !<  write for scalar (no-spin-orbit)
  logical, intent(in)                ::  lso                             !<  write for spin-orbit
  integer, intent(in)                ::  identif                         !<  identifier that is almost unique to the calculation

  integer, intent(in)                ::  nrk                             !<  number of irreducible k-points
  integer, intent(in)                ::  nx, ny, nz                      !<  original k-point mesh

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)
  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(in)                ::  nband(mxdnrk)                   !<  number of bands for each k-points
  real(REAL64), intent(in)           ::  rk(3,mxdnrk)                    !<  component in lattice coordinates of the k-point in the mesh
  real(REAL64), intent(in)           ::  wgk(mxdnrk)                     !<  weight in the integration of k-point
  integer, intent(in)                ::  indk(6,mxdnrk)                  !<  index of the six k-points neighbouring k-point i
  integer, intent(in)                ::  kmap(3,nx,ny,nz)                !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation

  real(REAL64), intent(in)           ::  e_of_k(mxdbnd,nrk)              !<  band energies of k-point in plot
  real(REAL64), intent(in)           ::  e_of_k_so(2*mxdbnd,nrk)         !<  spin-orbit band energies of k-point in plot

! counters

  integer                            ::  irk


  open(unit=io, file= trim(filename), form='UNFORMATTED')

  write(io) nrk
  write(io) nband(1:nrk)
  write(io) nx,ny,nz

  write(io) title, subtitle
  write(io) lso, lscl, identif

  write(io) ztot
  write(io) adot

  write(io) ntrans
  write(io) mtrx

  write(io) rk(:,1:nrk)
  write(io) wgk(1:nrk)
  write(io) indk(:,1:nrk)
  write(io) kmap

  
  if(lscl) then
    do irk = 1, nrk
      write(io) e_of_k(1:nband(irk),irk)
    enddo
  endif

  if(lso) then
    do irk = 1, nrk
      write(io) e_of_k_so(1:2*nband(irk),irk)
    enddo
  endif

  close(unit=io)

  return
  end subroutine out_dos_write
