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

!>  sorts an array by the heapsort method
!>  and performs the same permutation on one real and
!>  one integer associated arrays
!>  adapted from http://rosettacode.org
!>  see also W. H. Preuss et al. Numerical Recipes
!>
!>  \author       Jose Luis Martins
!>  \version      5.01
!>  \date         20 April 2021
!>  \copyright    GNU Public License v2

subroutine heapsort_par(n, a, b, ia)

! I was too lazy to write the sort in place...

  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)           ::  n                                    !<  length of array

! input and output

  real(REAL64), intent(inout)   ::  a(n)                                 !<  array to be sorted
  real(REAL64), intent(inout)   ::  b(n)                                 !<  real array to be reindexed
  integer, intent(inout)        ::  ia(n)                                !<  integer array to be reindexed

! allocatable arrays

  integer, allocatable          ::  indx(:)
  real(REAL64), allocatable     ::  aux(:)
  integer, allocatable          ::  iaux(:)

! counter

  integer     ::  i

  allocate(indx(n))
  allocate(aux(n))
  allocate(iaux(n))

  call sort(n,a,indx)

  do i = 1,n
    aux(i) = a(i)
  enddo
  do i = 1,n
    a(i) = aux(indx(i))
  enddo

  do i = 1,n
    aux(i) = b(i)
  enddo
  do i = 1,n
    b(i) = aux(indx(i))
  enddo

  do i = 1,n
    iaux(i) = ia(i)
  enddo
  do i = 1,n
    ia(i) = iaux(indx(i))
  enddo

  deallocate(indx)
  deallocate(aux)
  deallocate(iaux)

end subroutine heapsort_par
