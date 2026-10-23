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

!>  Symmetry check of the interpolated modified tight-binding
!>  hamiltonian:  largest difference of the eigenvalues in the star
!>  (with time reversal) of a general k-point.
!>
!>  \author       claude
!>  \version      5.13
!>  \date         23 October 2026.
!>  \copyright    GNU Public License v2

subroutine ao_int_sym_check(noiData, ntrans, mtrx, xsymerr)

! Written 23 October 2026 from ao_interpolation_prepare. claude

  use NonOrthoInterp

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  type(noiData_t)                    ::  noiData                         !<  see NonOrthoInterp

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group

! output

  real(REAL64), intent(out)          ::  xsymerr                         !<  largest eigenvalue difference in the star (hartree)

! local allocatable arrays

  real(REAL64), allocatable          ::  ev0(:)                          !  eigenvalues at the k-point
  real(REAL64), allocatable          ::  ev(:)                           !  eigenvalues at a symmetry related k-point

! local variables

  real(REAL64)      ::  rksym(3)                                         !  k-point for the symmetry test
  real(REAL64)      ::  rkpt(3)                                          !  symmetry related k-point
  integer           ::  itr                                              !  time reversal

! parameters

  real(REAL64), parameter :: ZERO = 0.0_REAL64
  real(REAL64), parameter :: HARTREE = 27.211386246_REAL64

! counters

  integer           ::  isym, j


  allocate(ev0(noiData%nband))
  allocate(ev(noiData%nband))

  rksym(1) = 0.1234_REAL64
  rksym(2) = 0.2345_REAL64
  rksym(3) = 0.3456_REAL64
  call NonOrthoInterpRun(noiData, rksym, ev0)

  xsymerr = ZERO
  do isym = 1,ntrans
    do itr = 1,-1,-2
      do j = 1,3
        rkpt(j) = itr*(mtrx(j,1,isym)*rksym(1) +                         &
                       mtrx(j,2,isym)*rksym(2) +                         &
                       mtrx(j,3,isym)*rksym(3))
      enddo
      call NonOrthoInterpRun(noiData, rkpt, ev)
      do j = 1,noiData%nband
        xsymerr = max(xsymerr,abs(ev(j) - ev0(j)))
      enddo
    enddo
  enddo

  write(6,*)
  write(6,'("   Largest eigenvalue difference in the star of a ",        &
      &     "k-point (symmetry check): ",e10.2," eV")') xsymerr*HARTREE
  write(6,*)

  deallocate(ev0)
  deallocate(ev)

  return

end subroutine ao_int_sym_check
