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

!>  Identifies the inverse symmetry operation.
!>
!>  \author       Jose Luis Martins
!>  \version      5.06
!>  \date         11 November 2022.
!>  \copyright    GNU Public License v2

subroutine sym_inverse_op(invop,                                         &
     ntrans, mtrx, tnp)

! Written 11 November 2022, loosely based on sym_test.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  ntrans                          !<  number of symmetry operations in the factor group
  integer, intent(in)                ::  mtrx(3,3,48)                    !<  rotation matrix (in reciprocal lattice coordinates) for the k-th symmetry operation of the factor group
  real(REAL64), intent(in)           ::  tnp(3,48)                       !<  2*pi* i-th component (in lattice coordinates) of the fractional translation vector associated with the k-th symmetry operation of the factor group

! output

  integer, intent(out)               ::  invop(48)                       !<  invop(j) identifies the inverse operation of j

! local variables

  integer        ::  itmp(3,3)
  logical        ::  lnull
  integer        ::  ninv(48)
  real(REAL64)   ::  frac(3)

! parameters

  real(REAL64), parameter :: PI = 3.14159265358979323846_REAL64
  real(REAL64), parameter :: EPS = 1.0E-8_REAL64

! counters

  integer i, j, k, n, m

! paranoid check

  if(ntrans < 1 .or. ntrans > 48) then
    write(6,*)
    write(6,*) '   stopped in sym_inverse_op'
    write(6,*)
    write(6,'("    ntrans = ",i10)') ntrans
    write(6,*)

    stop

  endif

! initial calculations

  if(ntrans == 1) then
    invop(1) = 1
  else
    do n = 1,ntrans
      do m = n,ntrans

        do i = 1,3
        do j = 1,3
          itmp(i,j) = 0
          do k = 1,3
            itmp(i,j) = itmp(i,j) + mtrx(i,k,n)*mtrx(k,j,m)
          enddo
        enddo
        enddo
        do i = 1,3
          itmp(i,i) = itmp(i,i) - 1
        enddo

        lnull = .TRUE.
        do i = 1,3
          do j = 1,3
            if(itmp(i,j) /= 0) then
              lnull = .FALSE.

              exit

            endif
          enddo

          if(lnull == .FALSE.) exit

        enddo

        if(lnull) then
          invop(n) = m
          invop(m) = n

          exit

        endif

      enddo
    enddo

  endif

! paranoid checks

  do n = 1,ntrans
    ninv(n) = 0
  enddo

  do n = 1,ntrans
    ninv(invop(n)) = ninv(invop(n)) + 1
  enddo

  do n = 1,ntrans
    if(ninv(n) /= 1) then
      write(6,*)
      write(6,*) '   stopped in sym_inverse_op'
      write(6,*)
      write(6,'(48i5)') (invop(i),i=1,ntrans)
      write(6,*)

      stop

    endif
  enddo

  do n = 1,ntrans

    do i = 1,3
      frac(i) = mtrx(1,i,invop(n))*tnp(1,invop(n)) +                              &
                mtrx(2,i,invop(n))*tnp(2,invop(n)) +                              &
                mtrx(3,i,invop(n))*tnp(3,invop(n)) + tnp(i,n)
      frac(i) = frac(i) / (2*PI)
    enddo
    do i = 1,3
      if(abs(frac(i)-nint(frac(i))) > EPS) then
        write(6,*)
        write(6,*) '   stopped in sym_inverse_op'
        write(6,*)
        write(6,'("  op = ",i5,"  i = ",i5,4x,f12.4)') n, i, frac(i)
        write(6,*)

        stop

      endif
    enddo
  enddo

  return

end subroutine sym_inverse_op
