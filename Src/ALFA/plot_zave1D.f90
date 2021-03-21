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

!>  calcultes a function in real space, averages in the "xy" direction
!>  and plots it in the "z"direction.

subroutine plot_zave1D(ave, func, nplane, id, n1,n2,n3, ng, kgv)

! Writen June 4, 2014.jlm
! Modified, documentation, name 4 February 2021. JLM
! copyright  Jose Luis Martins/INESC-MN

  implicit none

! version 4.99

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ng                              !<  total number of g-vectors with length less than gmax
  integer, intent(in)                ::  kgv(3,ng)                       !<  i-th component (reciprocal lattice coordinates) of the n-th g-vector ordered by stars of increasing length
  complex(REAL64), intent(in)        ::  func(ng)                        !<  function in reciprocal space
  integer, intent(in)                ::  nplane                          !<  number of atomic planes
  integer, intent(in)                ::  id, n1, n2, n3                  !<  dimensions for the FFT

! output

  real(REAL64), intent(out)          ::  ave(n3)                         !<  layer average of func in the 3d direction of fft grid
  

! allocatable arrays

  complex(REAL64), allocatable       :: chd(:,:,:)                  !  charge density in the fft grid
  complex(REAL64), allocatable       :: wrk(:)                      !  work array.

! other variables

  integer      ::  mxdwrk
  integer      ::  k1, k2, k3, kd

! constants

  real(REAL64), parameter  :: ZERO = 0.0_REAL64

! counters

  integer      ::  i, j, k


! allocates and initialize charge density array

  
  mxdwrk = 2*max(n1*n2,n1*n3,n2*n3)

  allocate(chd(id,n2,n3))
  allocate(wrk(mxdwrk))

  do k=1,n3
  do j=1,n2
  do i=1,id
    chd(i,j,k) = cmplx(ZERO,ZERO,REAL64)
  enddo
  enddo
  enddo

! wraparound for charge

  do i=1,ng
    k1 = kgv(1,i)
    kd = n1*(k1/n1)
    if (k1 < 0) kd = kd - n1
    k1 = k1 - kd + 1
    k2 = kgv(2,i)
    kd = n2*(k2/n2)
    if (k2 < 0) kd = kd - n2
    k2 = k2 - kd + 1
    k3 = kgv(3,i)
    kd = n3*(k3/n3)
    if (k3 < 0) kd = kd - n3
    k3 = k3 - kd + 1

    chd(k1,k2,k3) = chd(k1,k2,k3) + func(i)
  enddo

  call cfft_c16(chd, id, n1,n2,n3, -1, wrk, mxdwrk)

  do k=1,n3
    ave(k) = ZERO
    do j=1,n2
    do i=1,n1
      ave(k) = ave(k) + real(chd(i,j,k),REAL64)
    enddo
    enddo
    ave(k) = ave(k) / (n1*n2*nplane)
  enddo

  deallocate(chd)
  deallocate(wrk)


  return

end subroutine plot_zave1D
