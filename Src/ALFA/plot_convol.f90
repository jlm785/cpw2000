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

!>  performs a convolution with a square function
!>
!>  \author       Jose Luis Martins
!>  \version      5.07
!>  \date         June 2014.
!>  \copyright    GNU Public License v2

subroutine plot_convol(n, x, a, c)

! Written possibly in June 2014. JLM
! Documentation, 4 February 2021. JLM
! Documentation details, 3 August 2023.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: INT64 = selected_int_kind(18)

! input

  integer, intent(in)         :: n                                       !<  number of points
  real(REAL64), intent(in)    :: x                                       !<  width of square function
  real(REAL64), intent(in)    :: a(n)                                    !<  original function

! output

  real(REAL64), intent(out)   :: c(n)                                    !<  function convoluted with square function

! local variables

!  complex(REAL64), allocatable  :: b(:,:)        !
!  complex(REAL64), allocatable  :: wrk(:)        !
!  complex(REAL64), allocatable  :: trigs(:)      !

  real(REAL64), allocatable  ::  b(:,:)
  real(REAL64), allocatable  ::  wrk(:,:)
  real(REAL64), allocatable  ::  trigs(:)

  integer  :: ifax(19)                           !
  integer  :: plan2
  integer(INT64) :: plan3

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64
  real(REAL64), parameter  :: PI=3.141592653589793_REAL64

! counters

  integer ::   i,j


  allocate(b(2,0:n-1))
  allocate(wrk(2,n))
  allocate(trigs(4*n))

  do i=1,n
!    b(i-1,1) = cmplx(a(i),ZERO,REAL64)
    b(1,i-1) = a(i)
    b(2,i-1) = ZERO
  enddo

  call cfft_prepare(b,wrk,trigs,ifax,1,n,n,1,1,plan2,plan3)

  call cfft_mlt_c16(b,wrk,trigs,ifax,1,n,n,1,1,plan2,plan3)

  call cfft_finish(plan2,plan3)

  do i=1,n-1
    j = mod(i+n/2,n) - n/2
!    b(i,1) = b(i,1)*sin(PI*x*j)/(PI*x*j)
    b(1,i) = b(1,i)*sin(PI*x*j)/(PI*x*j)
    b(2,i) = b(2,i)*sin(PI*x*j)/(PI*x*j)
  enddo

  call cfft_prepare(b,wrk,trigs,ifax,1,n,n,1,-1,plan2,plan3)

  call cfft_mlt_c16(b,wrk,trigs,ifax,1,n,n,1,-1,plan2,plan3)

  do i=1,n
!    c(i) = real(b(i-1,1))/n
    c(i) = b(1,i-1) / n
  enddo

  call cfft_finish(plan2,plan3)

  deallocate(b)
  deallocate(wrk)
  deallocate(trigs)

  return

end subroutine plot_convol

