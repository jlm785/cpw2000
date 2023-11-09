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

!>  Performs the interpolation for the calculation the effective mass
!>  for a given k-vector and direction with k.p model
!>
!>  \author       Jose Luis Martins
!>  \version      5.08
!>  \date         8 November 2023.
!>  \copyright    GNU Public License v2

subroutine out_mass_kdotp_xk(rkpt, xk, nmodel, npt, delta,               &
       deidk, d2eidk2,                                                   &
       adot, h0, dh0drk, d2h0drk2,                                       &
       mxdbnd)


! Extracted from out_effective_mass when it was converted to out_mass_kdotp

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  mxdbnd                          !<  array dimension for the number of bands

  real(REAL64), intent(in)           ::  rkpt(3)                         !<  k-point reciprocal lattice coordinates
  real(REAL64), intent(in)           ::  xk(3)                           !<  k-direction in reciprocal lattice coordinates
  integer, intent(in)                ::  npt                             !<  2*npt+1 is the total number of interpolation points
  real(REAL64), intent(in)           ::  delta                           !<  spacing between the poins used in the interpolation

  integer, intent(in)                ::  nmodel                          !<  size of the kdotp matrix

  real(REAL64), intent(in)           ::  adot(3,3)                       !<  metric in direct space

  complex(REAL64), intent(in)        ::  h0(mxdbnd,mxdbnd)               !<  <Psi|H|Psi> without spin-orbit
  complex(REAL64), intent(in)        ::  dh0drk(mxdbnd,mxdbnd,3)         !<  d <Psi|H|Psi> d k
  complex(REAL64), intent(in)        ::  d2h0drk2  (mxdbnd,mxdbnd,3,3)   !<  d^2 <Psi|H|Psi> d k^2


! output

  real(REAL64), intent(out)          ::  deidk(mxdbnd)                   !<  d E / d xk  (lattice coordinates)
  real(REAL64), intent(out)          ::  d2eidk2(mxdbnd)                 !<  d^2 E / d xk^2  (lattice coordinates)

! allocatable arrays

  real(REAL64), allocatable          ::  ei_l(:,:)                       !  eigenvalue no. i. in the line (hartree)
  real(REAL64), allocatable          ::  rk_l(:,:)                       !  k-point on the line

  real(REAL64), allocatable          ::  xin(:), yin(:)

  integer, allocatable               ::  ipl(:)                          !  most similar levels.


! local variables

  real(REAL64)      ::  yk(3), xknorm
  real(REAL64)      ::  vcell, bdot(3,3)

  real(REAL64)      ::  y(0:2), dy(0:2)

! constants

  real(REAL64), parameter     ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter     ::  EPS = 1.0E-12_REAL64

! counters

  integer    ::  i, j, k, n


! normalizes xk

  call adot_to_bdot(adot, vcell, bdot)

  xknorm = ZERO
  do i = 1,3
  do j = 1,3
    xknorm = xknorm + xk(i)*bdot(i,j)*xk(j)
  enddo
  enddo

  if(xknorm < EPS) then
    write(6,*)
    write(6,*) '   stopped in out_mass_kdotp_xk, xknorm = ', xknorm
    write(6,*)

    stop

  endif

  do i = 1,3
    yk(i) = xk(i) / sqrt(xknorm)
  enddo

  allocate(ei_l(mxdbnd,-npt:npt))
  allocate(rk_l(3,-npt:npt))

  do n = -npt,npt

    do k = 1,3
      rk_l(k,n) = rkpt(k) - n*delta*yk(k)
    enddo

    call kdotp_diag_nopsi(rk_l(:,n), nmodel, ei_l(:,n),                  &
             rkpt, h0, dh0drk, d2h0drk2,                                 &
             mxdbnd)

  enddo

! tries to match left and right sides

  allocate(ipl(nmodel))

  call out_mass_match(nmodel, npt, ei_l, ipl,                            &
       mxdbnd)

! does the fit

  allocate(xin(-npt:npt))
  allocate(yin(-npt:npt))

  do n = -npt,npt
    xin(n) = n*delta
  enddo

  do n = 1,nmodel

    do j = -npt,-1
      yin(j) = ei_l(ipl(n),j) - ei_l(n,0)
    enddo

    do j = 0,npt
      yin(j) = ei_l(n,j) - ei_l(n,0)
    enddo

    call poly_interp(y, dy, xin, yin, 2*npt, 2)

    d2eidk2(n) = y(2)
    deidk(n) = y(1)

  enddo

  deallocate(ei_l)
  deallocate(rk_l)

  deallocate(xin,yin)
  deallocate(ipl)

  return

end subroutine out_mass_kdotp_xk
