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

!>  generates an ~8 digit integer from time data
!>  can be used as an identifier, although it is not foolproof.

subroutine get_identifier(identif)

! Written 8 December 2020. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.99

  implicit none
  integer, parameter          :: REAL64 = selected_real_kind(12)

! output

  integer, intent(out)          ::  identif                              !<  almost ynique identifier

! local variables

  real(REAL64)                  ::  t1, t2, x
  integer                       ::  iseed, nalg
  real(REAL64),external      :: ran2

  call zeelap(t1)

  t1 = t1 + 123456.0
  nalg = nint(log10(t1))
  t1 = t1 * 10.0**(8-nalg)

  call zesec(t2)
  nalg = nint(log10(t2))
  t2 = t2 * 10.0**(6-nalg)

  iseed = -nint(abs(t2)+1234)
  x = ran2(iseed)

  identif = nint(x*t1)

  return
end subroutine get_identifier
