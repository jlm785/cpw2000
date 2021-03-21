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

!>  Writes the files with the oscillator strengths for later processing.

  subroutine out_opt_write(title, subtitle, identif,                     &
    neig, nval, ztot, adot, ntrans, mtrx,                                &
    nrk, rk, wght, indk, kmap,                                           &
    mxdbnd, mxdpnt, nx,ny,nz)

! Extracted from out_opt_ie, 7 December 2020. JLM

! copyright  Carlos Loia Reis/Jose Luis Martins/INESC-MN

! version 4.99

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands
  integer, intent(in)                ::  nx,ny,nz                        !<  k-point mesh size
  integer, intent(in)                ::  mxdpnt                          !<  dimensions for dos k-points

  character(len=50), intent(in)      ::  title                           !<  title for plots
  character(len=140), intent(in)     ::  subtitle                        !<  subtitle for plots

  integer, intent(in)                ::  identif                         !<  identifier that is almost unique to the calculation

  integer, intent(in)                ::  neig                            !<  total number of eigenvectors
  integer, intent(in)                ::  nval                            !<  number of valence states

  real(REAL64), intent(in)           ::  ztot                            !<  total charge density (electrons/cell)

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space
  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

  integer, intent(in)                ::  nrk                             !<  number of k-points
  real(REAL64), intent(in)           ::  rk(3,mxdpnt)                    !<  component in lattice coordinates of the k-point in the mesh
  real(REAL64), intent(in)           ::  wght(mxdpnt)                    !<  weight in the integration of k-point

  integer, intent(in)                ::  kmap(3,nx,ny,nz)                !<  kmap(1,...) corresponding k-point, kmap(2,...) = 1 additional inversion, kmap(3,...) symmetry operation
  integer, intent(in)                ::  indk(6,mxdpnt)                  !<  index of the six k-points neighbouring k-point i
  
! local variables

  integer                            ::  io21

  integer                            ::  ncond


  io21 = 21

  ncond = neig - nval

  open(unit=io21, file="opt_struct.dat", form="unformatted")

  write(io21) title, subtitle
  write(io21) identif

  write(io21) nx,ny,nz
  write(io21) ntrans
  write(io21) mtrx
  write(io21) nrk
  write(io21) rk
  write(io21) wght
  write(io21) indk
  write(io21) kmap

  write(io21) ztot
  write(io21) neig, nval, ncond
  write(io21) adot

  write(io21) mxdbnd

  close(io21)

  return
  end subroutine out_opt_write
