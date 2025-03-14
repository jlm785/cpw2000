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

!>  Given a plot step it suggests a nicer step  b = 1/2/5 * 10^n and b~a
!>
!>  \author       Jose Luis Martins
!>  \version      5.11
!>  \date         15 November 2013. 13 March 2025.
!>  \copyright    GNU Public License v2

subroutine plot_step(a,b)

! written 15 November 2013, based on older code. jlm
! Modified documentation August 2019.  JLM
! Modified, Indentation, a <0. 13 March 2025. JLM


  implicit none
  integer, parameter  :: REAL64 = selected_real_kind(12)

! input
  real(REAL64), intent(in)   ::  a                                  !<  approximatestep desired

! output
  real(REAL64), intent(out)  ::  b                                  !<  b>a and b is a round number

! local variables
  integer          :: ixe
  real(REAL64)     :: xl, xm, xe

! constants
  real(REAL64), parameter :: ZERO = 0.0_REAL64 , UM = 1.0_REAL64

  if(a < tiny(UM)) then
    b = UM
  else
    xl = log10(a)
    ixe = int(xl)
    if(ixe < 0) ixe = ixe - 1
    xe = real(ixe)
    xm = xl - xe
    if(xm < 0.15) then
      xm = ZERO
    else if(xm < 0.5) then
      xm = log10(2*UM)
    else if(xm < 0.85) then
      xm = log10(5*UM)
    else
      xm = UM
    endif
    xl = xe + xm
    b = ((10*UM)**xl)
  endif

  return

end subroutine plot_step
